#include "CavityDetection.h"

int CavityDetection::InitVDB(int dss, double vs)
{
	openvdb::initialize();

	_proteinGrid = openvdb::DoubleGrid::create();

	_voxelSize = vs;
	_halfWidth = 3;
	deformStepSize = dss;

	return 1;
}

int CavityDetection::InitMoleculerSurface(std::string msFile)
{
	ReadMolecularSurfaceTriangleMesh(msFile);

	_xform = openvdb::math::Transform::createLinearTransform(_voxelSize);
	_proteinGrid = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(*_xform, _meshVertices, _meshFaces, _halfWidth);

	int tav; // total number of active voxels. We want it to be near 100k, otherwise we ajust the voxel size
	SelectProteinInterested(tav);
	_voxelSize /= cbrt(100000.0 / tav);
	
	_proteinGrid->clear();
	_xform = openvdb::math::Transform::createLinearTransform(_voxelSize);
	_proteinGrid = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(*_xform, _meshVertices, _meshFaces, _halfWidth);
	
	_proteinInterested.clear();
	SelectProteinInterested(tav);
	
	FillColorVector();
	GetCenterAndDist();

	return 1;
}

int CavityDetection::FirstPass()
{
	_membraneGrid = _OriMembraneGrid->deepCopy();
	_numBVM = 0;
	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
		openvdb::Coord coord = iter.getCoord();
		_blocked[coord] = false;
		++_numBVM;
	}

	InitProteinBlocked();

	Deform();

	return 1;
}

int CavityDetection::SecondPass()
{
	Clear();

	_membraneGrid->clear();
	_membraneGrid = _OriMembraneGrid->deepCopy();
	_numBVM = 0;
	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
		openvdb::Coord coord = iter.getCoord();
		_blocked[coord] = false;
		++_numBVM;
	}

	InitProteinBlocked();
	InitFinalCCLabel();
	DeformToMarkRepresentativePockets();

	return 1;
}

int CavityDetection::InitMembrane()
{
	GetConvexHull();
	BuildAndRelaxConvexHull();

	return 1;
}

int CavityDetection::setParameters(double d, double a, double v)
{
	m_moveSpeed = 0.5; // ratio of voxel size
	m_persistenceThreshold = d; // length in A
	m_areaMinimum=a; // A^2
	m_depthMinimum = d; // A
	m_volumeMinimum = v; // A^3

	return 1;
}

int CavityDetection::countActiveVoxels()
{
	_totalNumOfActiveVoxels = 0;
	for (openvdb::DoubleGrid::ValueOnIter iter = _OriMembraneGrid->beginValueOn(); iter.test(); ++iter)
		++_totalNumOfActiveVoxels;
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter)
		++_totalNumOfActiveVoxels;
	
	return 1;
}

int CavityDetection::Deform()
{
	persistenceGraph = new PGraph();

	double step = 50;
	double forceRatio = m_moveSpeed;
	int numDeform = 0;
	int forceCounter = 0;

	for (int i = 0; i < deformStepSize; i++, ++numDeform){
		std::cout << "Deforming setep: " << i << std::endl;
		if (i % 3 == 2){
			ReinitializeLevelSet();
		}
		MeanCurvatureDeformation(step, forceRatio, forceCounter);
		
		MoveInfo();
		SelectMembraneInterested();
		LabelConnectedComponents();
		AppendPersistenceGraph(i);

		if (AllProteinBlocked())
			break;
	}

	CheckPersistenceGraphCompleteness();
	AssembleChildren();
	AssemblePGraphData();

	AssembleHomologyComponent(numDeform);

	ReconstructMajorComponentTree();

	return 1;
}

int CavityDetection::DeformToMarkRepresentativePockets()
{
	InitLabel();

	std::sort(reps.begin(), reps.end(), PNICompare());

	int finalLabel = 0;
	_totalCCNum = 0;
	double step = 50;
	double forceRatio = m_moveSpeed;
	int forceCounter = 0;
	for (int i = 0; i < deformStepSize; i++){
		std::cout << "Deforming setep: " << i << std::endl;
		if (i % 3 == 2){
			ReinitializeLevelSet();
		}

		MeanCurvatureDeformation(step, forceRatio, forceCounter);

		std::vector<int> wantedLabels;
		std::map<int, double> labelToDetph;
		SelectStepLabelPair(wantedLabels, labelToDetph, i, reps);

		LabelCCAndExtractReps(wantedLabels, labelToDetph, i, finalLabel);

		if (AllProteinBlocked())
			break;
	}

	return 1;
}

int CavityDetection::PocketAtomMap(std::string xyzrFile)
{
	// read position
	std::fstream file(xyzrFile.c_str());
	std::string str;

	while (std::getline(file, str)){
		std::stringstream ss(str);
		double x, y, z;
		ss >> x >> y >> z;

		_atomLocations.push_back(openvdb::Vec3d(x, y, z));
	}

	std::cout << "Reading XYZR file complete... " << _atomLocations.size() << std::endl;

	// assemble kd tree
	PointCloud<double> cloud;
	cloud.pts.resize(_atomLocations.size());
	for (int i = 0; i < _atomLocations.size(); ++i){
		cloud.pts[i].x = _atomLocations[i].x();
		cloud.pts[i].y = _atomLocations[i].y();
		cloud.pts[i].z = _atomLocations[i].z();
	}

	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double> >, PointCloud<double>, 3> MyKDTree;
	MyKDTree kdTree(3, cloud, KDTreeSingleIndexAdaptorParams(20));
	kdTree.buildIndex();

	std::cout << "Building KDTree complete..." << std::endl;

	for (int i = 0; i < _totalCCNum; ++i) {
		std::set<size_t> indices;
		std::cout << _finalLabelToCC[i].size() << std::endl;
		for (int j = 0; j < _finalLabelToCC[i].size(); ++j) {
			openvdb::Coord coord = _finalLabelToCC[i][j];
			openvdb::Vec3d wldCoord = _xform->indexToWorld(coord);
			const double query[3] = { wldCoord.x(), wldCoord.y(), wldCoord.z() };
			int num = 1;
			std::vector<size_t> r(num);
			std::vector<double> d(num);

			num = kdTree.knnSearch(&query[0], num, &r[0], &d[0]);
			if (r[0] < _atomLocations.size())
				indices.insert(r[0]);
		}
		_CCLabelToAtoms[i] = indices;
	}

	//openvdb_viewer::Viewer viewer = openvdb_viewer::init("Viewer", false);
	//View(viewer, 0);

	return 1;
}

int CavityDetection::ReadMolecularSurfaceTriangleMesh(std::string msFile)
{
	std::ifstream in;
	in.open(msFile.c_str());

	ReadMeshData(in);

	return 1;
}

int CavityDetection::MeanCurvatureDeformation(double& tStep, double& forceRatio, int& i)
{
	//Test all blocked
	if (AllBlocked()){
		std::cout << "Warning: Membrane is all-blocked, but protein is still not all-blocked..."<< std::endl;
		exit(0);
	}

	//Estimate Laplacian
	int count = 0;
	omp_set_num_threads(8);
	#pragma omp parallel
	{
		for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
			#pragma omp single nowait
			{
				//std::cout << "Dealing with voxel: " << iter.getCoord() << std::endl;

				openvdb::Coord currentCoord;
				currentCoord = iter.getCoord();

				openvdb::Coord frontCoord, rearCoord, leftCoord, rightCoord, upCoord, downCoord;
				frontCoord = currentCoord; frontCoord.z()++;
				rearCoord = currentCoord; rearCoord.z()--;
				leftCoord = currentCoord; leftCoord.x()--;
				rightCoord = currentCoord; rightCoord.x()++;
				upCoord = currentCoord; upCoord.y()++;
				downCoord = currentCoord; downCoord.y()--;

				openvdb::DoubleGrid::Accessor accessor = _membraneGrid->getAccessor();
				double current;
				current = accessor.getValue(currentCoord);
				double front, rear, left, right, up, down;
				front = accessor.getValue(frontCoord);
				rear = accessor.getValue(rearCoord);
				left = accessor.getValue(leftCoord);
				right = accessor.getValue(rightCoord);
				up = accessor.getValue(upCoord);
				down = accessor.getValue(downCoord);

				double lap;
				lap = (front + rear - 2 * current) / pow(_voxelSize, 2)
					+ (right + left - 2 * current) / pow(_voxelSize, 2)
					+ (up + down - 2 * current) / pow(_voxelSize, 2);

				#pragma omp critical
				{
					_laplacian[currentCoord] = lap;
					count++;
				}
			}
		}
	}

	double backGround;
	backGround = _membraneGrid->background();

	//deform
	bool dispSuccess = false;

	dispSuccess = true;

	_displacement.clear();
	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter) {
		openvdb::Coord coord = iter.getCoord();
		_displacement[coord] = _laplacian[coord];
		if (_displacement[coord] < -0.5*forceRatio*_voxelSize)
			_displacement[coord] = -0.5*forceRatio*_voxelSize;
	}

	int totalActive = 0;
	int nearlyStable = 0;;
	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
		openvdb::Coord coord = iter.getCoord();
		double preVal= iter.getValue();

		if (!_blocked[coord]){
			++totalActive;

			if (_displacement[coord] > -0.0001) {
				iter.setValue(iter.getValue() + forceRatio * _voxelSize);
			}
			else {
				iter.setValue(iter.getValue() + _displacement[coord] + forceRatio * _voxelSize);
			}
		}

		openvdb::DoubleGrid::Accessor accessor= _proteinGrid->getAccessor();
		if (accessor.getValue(coord) < *iter - _voxelSize){
			iter.setValue(accessor.getValue(coord) + _voxelSize);
			_blocked[coord] = true;
		}

		if (iter.getValue() > backGround){
			iter.setValue(0.99*backGround);
		}
		if (iter.getValue() < -backGround){
			iter.setValue(-0.99*backGround);
		}
	}

	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		if (_proteinInterested[iter.getCoord()] && !_proteinBlocked[iter.getCoord()]){
			openvdb::DoubleGrid::Accessor accessor = _membraneGrid->getAccessor();

			if (accessor.getValue(iter.getCoord()) > *iter){
				_proteinBlocked[iter.getCoord()] = true;
			}
		}
	}

	return 1;
}

int CavityDetection::ReinitializeLevelSet()
{
	// level set
	openvdb::DoubleGrid::Ptr newMembraneGrid = openvdb::tools::levelSetRebuild<openvdb::DoubleGrid>(*_membraneGrid, 0, _halfWidth);

	// initialize blocked
	std::map<openvdb::Coord, bool> newBlocked;
	_numBVM = 0;
	for (openvdb::DoubleGrid::ValueOnIter newIter = newMembraneGrid->beginValueOn(); newIter.test(); ++newIter)
	{
		openvdb::Coord coord = newIter.getCoord();
		newBlocked[coord] = false;
		++_numBVM;
	}
	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter)
	{
		openvdb::Coord coord = iter.getCoord();
		openvdb::DoubleGrid::Accessor accessor = newMembraneGrid->getAccessor();

		if (accessor.isValueOn(coord) && _blocked[coord])
			newBlocked[coord] = true;
	}

	//switch
	_membraneGrid->clear();
	_membraneGrid = newMembraneGrid;
	_blocked.clear();
	_blocked = newBlocked;
	newBlocked.clear();

	return 1;
}

int CavityDetection::ReadMeshData(std::istream& in)
{
	std::string line;

	while (getline(in, line)){
		std::stringstream ss(line);
		std::string token;

		ss >> token;

		if (token == "v") ReadPosition(ss);
		if (token == "f") ReadFace(ss);
	}


	return 1;
}

int CavityDetection::ReadPosition(std::stringstream& ss)
{
	double x, y, z;
	ss >> x >> y >> z;

	_meshVertices.push_back(openvdb::Vec3s(1*x, 1*y, 1*z));

	return 1;
}

int CavityDetection::ReadFace(std::stringstream& ss)
{
	std::vector<int> faceIndices;
	std::string token;

	while (ss >> token){
		faceIndices.push_back(ParseFaceIndex(token));
	}

	_meshFaces.push_back(openvdb::Vec3I(faceIndices[0], faceIndices[1], faceIndices[2]));

	return 1;
}

int CavityDetection::ParseFaceIndex(const std::string& token)
{
	// parse indices of the form
	//
	// p/[t]/[n]
	//
	// where p is an index into positions, t is an index into
	// texcoords, n is an index into normals, and [.] indicates
	// that an index is optional

	std::stringstream in(token);
	std::string indexstring;
	int indices[3] = { -1, -1, -1 };
	int i = 0;

	while (getline(in, indexstring, '/'))
	{
		std::stringstream ss(indexstring);
		ss >> indices[i++];
	}

	// decrement since indices in OBJ files are 1-based
	return indices[0] - 1;
}

int CavityDetection::MoveToOrigin()
{
	//Calculate barycenter
	openvdb::Vec3s center(0, 0, 0);
	for (int i = 0; i < _meshVertices.size(); ++i)
		center += _meshVertices[i];
	center /= _meshVertices.size();

	for (int i = 0; i < _meshVertices.size(); ++i)
		_meshVertices[i] -= center;

	return 1;
}

int CavityDetection::GetConvexHull()
{
	std::vector<Point_3> points;
	for (int i = 0; i < _meshVertices.size(); ++i){
		Point_3 p(_meshVertices[i].x(), _meshVertices[i].y(), _meshVertices[i].z());
		points.push_back(p);
	}

	Polyhedron_3 poly;
	CGAL::convex_hull_3(points.begin(), points.end(), poly);
	CGAL::Polygon_mesh_processing::triangulate_faces(poly);
	
	std::map<Vertex_iterator, int> index;
	int idx = 0;
	for (Vertex_iterator vi = poly.vertices_begin(); vi != poly.vertices_end(); ++vi){
		_membraneVertices.push_back(openvdb::Vec3s(vi->point().x(), vi->point().y(), vi->point().z()));
		index[vi] = idx; ++idx;
	}

	for (Facet_iterator fi = poly.facets_begin(); fi != poly.facets_end(); ++fi){
		_membraneFaces.push_back(openvdb::Vec3I(index[fi->halfedge()->vertex()], 
												index[fi->halfedge()->next()->vertex()], 
												index[fi->halfedge()->next()->next()->vertex()]));
	}

	return 1;
}

int CavityDetection::BuildAndRelaxConvexHull()
{
	//build
	openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(_voxelSize);
	_OriMembraneGrid = openvdb::tools::meshToLevelSet<openvdb::DoubleGrid>(*xform, _membraneVertices, _membraneFaces, _halfWidth);

	return 1;
}

int CavityDetection::SelectProteinInterested(int& tav)
{
	_proteinInterested.clear();
	_proteinNumOfBoundaryVoxels = 0;

	int num=0;
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		num++;

		_proteinInterested[iter.getCoord()] = false;

		openvdb::Coord coord = iter.getCoord();
		for (int i = 0; i < 27; ++i){
			openvdb::Coord adjCoord = coord;
			int ii = i;

			int a, b, c;
			a = ii / 9; ii = ii - a * 9;
			b = ii / 3; ii = ii - b * 3;
			c = ii;

			if (a == 1)
				adjCoord.x()++;
			else if (a == 2)
				adjCoord.x()--;
			if (b == 1)
				adjCoord.y()++;
			else if (b == 2)
				adjCoord.y()--;
			if (c == 1)
				adjCoord.z()++;
			else if (c == 2)
				adjCoord.z()--;

			openvdb::DoubleGrid::Accessor accessor = _proteinGrid->getAccessor();
			if (accessor.getValue(adjCoord)*iter.getValue() <= 0)
				_proteinInterested[iter.getCoord()] = true;
		}

		if (_proteinInterested[iter.getCoord()])
			++_proteinNumOfBoundaryVoxels;
	}

	std::cout << "Total # active voxels of protein: " << num << std::endl;
	std::cout << "Total # boundary voxels of protein: " << _proteinNumOfBoundaryVoxels << std::endl;
	tav = num;

	return 1;
}

int CavityDetection::SelectMembraneInterested()
{
	_membraneInterested.clear();

	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
		_membraneInterested[iter.getCoord()] = false;

		openvdb::Coord coord = iter.getCoord();
		for (int i = 0; i < 27; ++i){
			openvdb::Coord adjCoord = coord;
			int ii = i;

			int a, b, c;
			a = ii / 9; ii = ii - a * 9;
			b = ii / 3; ii = ii - b * 3;
			c = ii;
		
			if (a == 1)
				adjCoord.x()++;
			else if(a==2)
				adjCoord.x()--;
			if (b == 1)
				adjCoord.y()++;
			else if (b == 2)
				adjCoord.y()--;
			if (c == 1)
				adjCoord.z()++;
			else if (c == 2)
				adjCoord.z()--;

			openvdb::DoubleGrid::Accessor accessor = _membraneGrid->getAccessor();
			if (accessor.getValue(adjCoord)*iter.getValue() <= 0)
				_membraneInterested[iter.getCoord()] = true;
		}
	}

	return 1;
}

int CavityDetection::InitProteinBlocked()
{
	int num = 0;
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		if (_proteinInterested[iter.getCoord()]){
			_proteinBlocked[iter.getCoord()] = false;
			++num;
		}
	}

	return 1;
}

bool CavityDetection::AllProteinBlocked()
{
	bool all = true;

	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		if (_proteinInterested[iter.getCoord()]){
			if (!_proteinBlocked[iter.getCoord()]){
				all = false;
				break;
			}
		}
	}

	return all;
}

bool CavityDetection::AllBlocked()
{
	bool all = true;

	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
		if (!_blocked[iter.getCoord()]){
			all = false;
			break;
		}
	}

	return all;
}

int CavityDetection::OutPutTriangleMesh(openvdb::DoubleGrid::Ptr grid)
{
	std::vector<openvdb::Vec3s> points;
	std::vector<openvdb::Vec3I> triangles;
	std::vector<openvdb::Vec4I> quads;

	openvdb::tools::volumeToMesh<openvdb::DoubleGrid>(*grid, points, triangles, quads);

	std::ofstream out("mesh.obj");

	for (int i = 0; i < points.size(); ++i)
		out << "v " << points[i].x() << " " << points[i].y() << " " << points[i].z() << std::endl;

	for (int i = 0; i < triangles.size(); ++i)
		out << "f " << triangles[i].x() + 1 << " " << triangles[i].z() + 1 << " " << triangles[i].y() + 1 << std::endl;

	for (int i = 0; i < quads.size(); ++i){
		out << "f " << quads[i].x() + 1 << " " << quads[i].z() + 1 << " " << quads[i].y() + 1 << std::endl;
		out << "f " << quads[i].x() + 1 << " " << quads[i].w() + 1 << " " << quads[i].z() + 1 << std::endl;
	}

	out.close();

	return 1;
}

//int CavityDetection::View(openvdb_viewer::Viewer& viewer, int counter)
//{
//	openvdb::io::File outFile("mygrids.vdb");
//
//	openvdb::GridPtrVec outGrids;
//	
//	//outGrids.push_back(_membraneGrid);
//	outGrids.push_back(_proteinGrid);
//
//	outFile.write(outGrids);
//	outFile.close();
//
//	openvdb::GridCPtrVec allGrids;
//	openvdb::io::File inFile("mygrids.vdb");
//	inFile.open();
//	openvdb::GridPtrVecPtr inGrids = inFile.getGrids();
//	allGrids.insert(allGrids.end(), inGrids->begin(), inGrids->end());
//
//	std::cout << "Opening viewer\n";
//	viewer.open();
//
//	std::cout << "View all grids\n";
//	viewer.view(allGrids, _proteinInterested, _finalCCLabel, _proteinBlocked, colorVector, _xform, counter, center, dist);
//	//viewer.view(allGrids, _membraneInterested, _ccLabel, _blocked, colorVector, _xform, counter, center, dist);
//
//	return 1;
//}

//int CavityDetection::ViewPocket(openvdb_viewer::Viewer& pViewer, openvdb::DoubleGrid::Ptr pocket)
//{
//	openvdb::io::File outFile("pocketGrid.vdb");
//
//	openvdb::GridPtrVec outGrids;
//	outGrids.push_back(pocket);
//	outGrids.push_back(_membraneGrid);
//
//	outFile.write(outGrids);
//	outFile.close();
//
//	openvdb::GridCPtrVec allGrids;
//
//	openvdb::io::File inFile("pocketGrid.vdb");
//	inFile.open();
//	openvdb::GridPtrVecPtr inGrids = inFile.getGrids();
//	allGrids.insert(allGrids.end(), inGrids->begin(), inGrids->end());
//
//	std::cout << "Opening viewer\n";
//	pViewer.open();
//
//	std::cout << "View all grids\n";
//	pViewer.view(allGrids, _proteinInterested, _finalCCLabel, _proteinBlocked, colorVector, _xform, 0, center, dist);
//
//	return 1;
//}

int CavityDetection::GetBoundingBox()
{
	openvdb::CoordBBox boundingBox;
	boundingBox = _proteinGrid->evalActiveVoxelBoundingBox();

	std::cout << "Bounding Box: " << boundingBox << std::endl;

	return 1;
}

int CavityDetection::LabelConnectedComponents()
{
	// for all interesting membrane cells which are not blocked.
	// they are sperated by blocked ones.
	CCPos.clear();
	InitLabel();

	openvdb::Coord seed;
	int label = 0;
	openvdb::Vec3d pos;
	openvdb::DoubleGrid::ValueOnIter cIter = _proteinGrid->beginValueOn();

	while (!AllLabeled(seed, cIter)){
		int voxelNum = 0;
		int csVoxelNum = 0;

		SpanSeed(seed, label, pos, voxelNum, csVoxelNum);

		_labelToCCSize[label] = voxelNum;
		_labelToCCSizeMembrane[label] = csVoxelNum;

		// record info
		CCPos[label] = pos / voxelNum;

		label++;
	}

	_numCC = label;

	return 1;
}

int CavityDetection::InitLabel()
{
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		if (_proteinInterested[iter.getCoord()]){
			_ccLabel[iter.getCoord()] = -1;
		}
	}

	for (openvdb::DoubleGrid::ValueOnIter iter = _membraneGrid->beginValueOn(); iter.test(); ++iter){
		if (_membraneInterested[iter.getCoord()]){
			_ccLabelMembrane[iter.getCoord()] = -1;
		}
	}

	return 1;
}

bool CavityDetection::AllLabeled(openvdb::Coord& seed, openvdb::DoubleGrid::ValueOnIter& cIter)
{
	bool all = true;

	for (openvdb::DoubleGrid::ValueOnIter iter = cIter; iter.test(); ++iter){
		if (!_proteinBlocked[iter.getCoord()] && _proteinInterested[iter.getCoord()]){
			if (_ccLabel[iter.getCoord()] == -1){
				seed = iter.getCoord();
				all = false;
				break;
			}
		}
	}

	return all;
}

int CavityDetection::SpanSeed(openvdb::Coord seed, int label, openvdb::Vec3d& pos, int& voxelNum, int& csVoxelNum, bool test)
{
	bool csSeedFound = false;
	std::vector<openvdb::Coord> csSeedVec;

	std::list<openvdb::Coord> grayNodes;
	grayNodes.push_back(seed);
	// Every push, mark label
	_ccLabel[seed] = label;

	pos.setZero();
	voxelNum = 0;

	while (!grayNodes.empty()){
		openvdb::Coord coord = grayNodes.back();
		grayNodes.pop_back();

		if (_membraneInterested[coord]) {
			csSeedFound = true;
			csSeedVec.push_back(coord);
		}

		pos += coord.asVec3d();
		++voxelNum;

		// push adjacent interesting unblocked ones
		// check all direction
		for (int i = 0; i < 27; ++i){
			openvdb::Coord adjCoord = coord;
			int ii = i;

			int a, b, c;
			a = ii / 9; ii = ii - a * 9;
			b = ii / 3; ii = ii - b * 3;
			c = ii;

			if (a == 1)
				adjCoord.x()++;
			else if (a == 2)
				adjCoord.x()--;
			if (b == 1)
				adjCoord.y()++;
			else if (b == 2)
				adjCoord.y()--;
			if (c == 1)
				adjCoord.z()++;
			else if (c == 2)
				adjCoord.z()--;
			
			// test
			if (_proteinInterested[adjCoord] && !_proteinBlocked[adjCoord])
				if (_ccLabel[adjCoord] == -1){
					grayNodes.push_back(adjCoord);
					_ccLabel[adjCoord] = label;
				}
		}
	}

	if (csSeedFound) {
		// there are many seeds. find the piece with largest spanned component.
		std::set<openvdb::Coord> visited;
		std::vector<int> pieces;
		
		for (int i = 0; i < csSeedVec.size(); ++i) {
			std::set<openvdb::Coord>::iterator iter = visited.find(csSeedVec[i]);
			if (iter == visited.end()) {
				int num = 0;
				std::list<openvdb::Coord> csGrayNodes;
				csGrayNodes.push_back(csSeedVec[i]);
				visited.insert(csSeedVec[i]);

				while (!csGrayNodes.empty()) {
					openvdb::Coord coord = csGrayNodes.back();
					csGrayNodes.pop_back();

					++num;

					for (int i = 0; i < 27; ++i){
						openvdb::Coord adjCoord = coord;
						int ii = i;

						int a, b, c;
						a = ii / 9; ii = ii - a * 9;
						b = ii / 3; ii = ii - b * 3;
						c = ii;

						if (a == 1)
							adjCoord.x()++;
						else if (a == 2)
							adjCoord.x()--;
						if (b == 1)
							adjCoord.y()++;
						else if (b == 2)
							adjCoord.y()--;
						if (c == 1)
							adjCoord.z()++;
						else if (c == 2)
							adjCoord.z()--;

						// test
						if (_membraneInterested[adjCoord] && !_blocked[adjCoord]) {
							std::set<openvdb::Coord>::iterator it = visited.find(adjCoord);
							if (it == visited.end()){
								csGrayNodes.push_back(adjCoord);
								visited.insert(adjCoord);
							}
						}
					}
				}

				pieces.push_back(num);

			}
		
		}
		csVoxelNum = *std::max_element(pieces.begin(), pieces.end());
	}
	else{
		std::cout << "Cross section seed not found. Cannot calculate cross section area." << std::endl;
	}

	return 1;
}

int CavityDetection::MoveInfo()
{
	_preBlocked.clear();
	_preProteinBlocked.clear();
	_preMembraneInterested.clear();
	_preCCLabel.clear();
	_preLabelToNode.clear();
	_preLabelToCCSize.clear();

	_preBlocked = _blocked;
	_preProteinBlocked = _proteinBlocked;
	_preMembraneInterested = _membraneInterested;
	_preCCLabel = _ccLabel;
	_preNumCC = _numCC;
	_preLabelToNode = _labelToNode;
	_preLabelToCCSize = _labelToCCSize;

	return 1;
}

int CavityDetection::AppendPersistenceGraph(int step)
{
	if (step == 0){
		for (int i = 0; i < _numCC; ++i){
			PNodeIter pni = persistenceGraph->InerstNode(i, step, true, _labelToCCSize[i], _labelToCCSizeMembrane[i], CCPos[i], 0);
			_labelToNode[i] = pni;
		}
	}
	else{
		for (int i = 0; i < _numCC; ++i) {
			PNodeIter pni = persistenceGraph->InerstNode(i, step, false, _labelToCCSize[i], _labelToCCSizeMembrane[i], CCPos[i], _preNumCC);
			_labelToNode[i] = pni;
		}
		for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
			if (!_proteinBlocked[iter.getCoord()] && _proteinInterested[iter.getCoord()]){
				PNodeIter node;
				node = _labelToNode[_ccLabel[iter.getCoord()]];
				node->votes[_preCCLabel[iter.getCoord()]] += 1;
			}
		}
		
		for (int k=0; k<_numCC; ++k){
			PNodeIter node;
			node = _labelToNode[k];

			// find best vote
			int tarPreCCLabel;
			int numVotes;
			for (int i = 0; i < node->votes.size(); ++i){
				if (i == 0){
					tarPreCCLabel = 0;
					numVotes = node->votes[0];
				}
				else if (node->votes[i] > numVotes){
					numVotes = node->votes[i];
					tarPreCCLabel = i;
				}
			}

			PNodeIter preNode = _preLabelToNode[tarPreCCLabel];
			node->parent = preNode;
			preNode->isLeaf = false;

		}
	}

	return 1;
}

int CavityDetection::AssembleChildren()
{
	for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni){
		if(!pni->isRoot && pni->parent!= persistenceGraph->nodes.end())
			pni->parent->children.push_back(pni);
	}

	return 1;
}

int CavityDetection::CheckPersistenceGraphCompleteness()
{
	for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni){
		if (!pni->isRoot && pni->parent == persistenceGraph->nodes.end()){
			std::cout << "Warning: Track of CC relation lost!" << std::endl;
			std::cout << "Component " << pni->label << " at level " << pni->step << " is lost!" << std::endl;
		}
	}

	return 1;
}

int CavityDetection::EliminateSmallTree()
{
	bool cont = true;

	while (cont){
		PNodeIter tar;
		cont = false;
		for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni){
			if (!pni->isRoot && pni->parent == persistenceGraph->nodes.end()){
				tar = pni;
				cont = true;
				break;
			}
		}

		if (cont){
			RemoveTree(tar);
		}
	}

	return 1;
}

int CavityDetection::RemoveTree(PNodeIter root){
	std::stack<PNodeIter> s;
	s.push(root);

	while (!s.empty()){
		PNodeIter pni = s.top();
		s.pop();

		// push children
		std::cout << pni->children.size() << std::endl;
		for (int i = 0; i < pni->children.size(); ++i)
			s.push(pni->children[i]);

		persistenceGraph->nodes.erase(pni);

	}

	return 1;
}


int CavityDetection::EliminateNoise()
{
	// for voxels that are stuck.
	for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni){
		if (pni->isLeaf)
			pni->leafTreated = false;
	}

	bool cont = true;
	while (cont){
		PNodeIter tar;
		cont = false;
		for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni){
			if (pni->isLeaf && !pni->leafTreated){
				tar = pni;
				pni->leafTreated = true;
				cont = true;
				break;
			}
		}

		if (cont){
			PNodeIter parentTar;

			while (!tar->isRoot){
				parentTar = tar->parent;
				if (parentTar->children.size() != 1)
					break;

				if (parentTar->approxArea == tar->approxArea)
				{
					RemovePNode(tar);

					tar = parentTar;
				}
				else
					break;
			}
		}
	}

	return 1;
}

int CavityDetection::RemovePNode(PNodeIter pni)
{
	PNodeIter parent = pni->parent;

	parent->isLeaf = true;
	parent->leafTreated = true;

	for (int i = 0; i < parent->children.size(); ++i)
		if (parent->children[i] == pni)
			parent->children.erase(parent->children.begin() + i);

	persistenceGraph->nodes.erase(pni);

	return 1;
}

int CavityDetection::AssemblePGraphData()
{
	int k = 0;
	for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni, ++k){
		pni->idx = k;
		gLocations.push_back(pni->pos);
		gSteps.push_back(1.0*pni->step / persistenceGraph->depth);
	}

	for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni, ++k){
		for (int i = 0; i < pni->children.size(); ++i){
			gIndices.push_back(pni->idx);
			gIndices.push_back(pni->children[i]->idx);
		}
	}

	GetNodeHeight();

	return 1;
}

int CavityDetection::AssembleHomologyComponent(int numDeform)
{
	// assemble priority queue
	std::priority_queue<PNodeIter, std::vector<PNodeIter>, PNICompare> elements;

	for (PNodeIter pni = persistenceGraph->nodes.begin(); pni != persistenceGraph->nodes.end(); ++pni){
		pni->passed = false;
		elements.push(pni);
	}

	while (!elements.empty()){
		PNodeIter seed = elements.top();
		elements.pop();

		//trace back
		if (!seed->passed){
			HomologyComponent hc;
			TraceBack(seed, hc);

			_components.push_back(hc);
		}
	}

	EliminateShortPersistence(m_persistenceThreshold);

	return 1;
}

int CavityDetection::TraceBack(PNodeIter seed, HomologyComponent& hc)
{
	PNodeIter pni = seed;
	hc.killTime = seed->step;
	
	while(!pni->passed && !pni->isRoot){
		pni->passed = true;
		hc.path.push_front(pni);

		pni = pni->parent;
	}

	if (pni->isRoot){
		pni->passed = true;
		hc.path.push_front(pni);
	}
	else{
		hc.path.push_front(pni);
	}

	hc.bornTime = pni->step;

	return 1;
}

int CavityDetection::EliminateShortPersistence(double r)
{
	std::sort(_components.begin(), _components.end(), HCCompareDescend());

	int idx=0;
	for (int i = 0; i < _components.size(); ++i, ++idx){
		if (_components[i].lifeTime()*m_moveSpeed*_voxelSize < r)
			break;
	}

	//std::cout << idx<< std::endl;
	_components.erase(_components.begin() + idx, _components.end());

	return 1;
}

int CavityDetection::GetNodeHeight()
{
	RecursiveHeight(persistenceGraph->Root());
	return 1;
}

double CavityDetection::RecursiveHeight(PNodeIter pni)
{
	if (pni->children.empty()) {
		pni->height = 0;

		return 0;
	}
	else if (pni->csArea*pow(_voxelSize, 3) / (2 * _voxelSize) < m_areaMinimum) {
		pni->height = 0;

		return 0;
	}
	else {
		double maxHeight = 0;
		for (int i = 0; i < pni->children.size(); ++i) {
			double rHeight = RecursiveHeight(pni->children[i]) + 1;
			maxHeight = std::max(rHeight, maxHeight);
		}
		pni->height = maxHeight;

		return maxHeight;
	}
}

int CavityDetection::ReconstructMajorComponentTree()
{
	repPGraph = new PGraph();

	std::map<PNodeIter, PNodeIter, PNICompare> oldToNewPni;
	std::map<PNodeIter, PNodeIter, PNICompare> newToOldPni;

	for (int i = 0; i < _components.size(); ++i){
		if (i == 0){
			PNodeIter prePni = repPGraph->nodes.end();
			for (std::list<PNodeIter>::iterator iter = _components[i].path.begin(); iter != _components[i].path.end(); ++iter){
				PNodeIter pni = repPGraph->InsertNode(*iter);
				oldToNewPni[*iter] = pni;
				newToOldPni[pni] = *iter;
				pni->parent = prePni;
				prePni = pni;
			}
		}
		else{
			PNodeIter prePni = oldToNewPni[*(_components[i].path.begin())];
			
			for (std::list<PNodeIter>::iterator iter = ++_components[i].path.begin(); iter != _components[i].path.end(); ++iter){
				PNodeIter pni = repPGraph->InsertNode(*iter);
				oldToNewPni[*iter] = pni;
				newToOldPni[pni] = *iter;
				pni->parent = prePni;
				prePni = pni;
			}
		}
	}

	int kk = 0;
	for (PNodeIter iter = repPGraph->nodes.begin(); iter != repPGraph->nodes.end(); ++iter, ++kk){
		iter->idx = kk;
	}

	std::cout << "Construct children..." << std::endl;
	int initCCNum = 0;
	for (PNodeIter iter = repPGraph->nodes.begin(); iter != repPGraph->nodes.end(); ++iter){
		PNodeIter parent = iter->parent;

		if (parent != repPGraph->nodes.end()){
			parent->children.push_back(iter);
		}
	}

	std::stack<PNodeIter> branchStart;
	branchStart.push(repPGraph->Root());

	while (!branchStart.empty()){
		PNodeIter pni = branchStart.top();
		branchStart.pop();

		TraceBranchUntilBifurcation(pni);

		for (int i = 0; i < pni->children.size(); ++i){
			branchStart.push(pni->children[i]);

		}
	}


	for (int i = 0; i < reps.size(); ++i){
		std::cout <<"Representative:"<<  reps[i]->step << " " << reps[i]->label << " " 
			<< reps[i]->approxArea << " "<< reps[i]->csArea<<" "<< reps[i]->height <<" " << newToOldPni[reps[i]]->idx << std::endl;
	}

	return 1;
}

int CavityDetection::ExtractRepresentatives()
{
	std::stack<PNodeIter> branchStart;
	branchStart.push(repPGraph->Root());

	while (!branchStart.empty()){
		PNodeIter pni=branchStart.top();
		branchStart.pop();

		TraceBranchUntilBifurcation(pni);

		for (int i = 0; i < pni->children.size(); ++i)
			branchStart.push(pni->children[i]);
	}

	return 1;
}

int CavityDetection::TraceBranchUntilBifurcation(PNodeIter& pni)
{
	PNodeIter rPni = pni;
	int count = 0;
	double d = 0;
	while (pni->children.size() == 1){
		if (d > persistenceGraph->depth*0.5) {
			--count;
			rPni = rPni->children[0];
		}
			
		pni = pni->children[0];

		if (pni->approxArea > m_areaMinimum)
			++count;
		
		++d;
	}

	if (pni->children.empty() && rPni->height > m_depthMinimum && rPni->approxArea*1.0 / _proteinNumOfBoundaryVoxels < 0.1){
		rPni->depth = rPni->height;
		reps.push_back(rPni);
	}


	return 1;
}

int CavityDetection::InitFinalCCLabel()
{
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		if (_proteinInterested[iter.getCoord()]){
			_finalCCLabel[iter.getCoord()] = -1;
		}
	}

	return 1;
}

int CavityDetection::SelectStepLabelPair(std::vector<int>& wantedLabels, std::map<int, double>& labelToDepth, int step, std::vector<PNodeIter> reps)
{
	for (int i = 0; i < reps.size(); ++i)
		if (reps[i]->step == step){
			wantedLabels.push_back(reps[i]->label);
			labelToDepth[reps[i]->label] = reps[i]->depth;
		}
	std::sort(wantedLabels.begin(), wantedLabels.end());

	return 1;
}

int CavityDetection::LabelCCAndExtractReps(std::vector<int> wantedLabels, std::map<int, double> labelToDepth, int step, int& finalLabel)
{
	InitLabel();

	openvdb::Coord seed;
	int label = 0;
	int current = 0;
	openvdb::DoubleGrid::ValueOnIter cIter = _proteinGrid->beginValueOn();

	// if wanted is not 0, then we need to extract volume.
	// 1. max to find volume cc
	// 2. pair and count #voxels as volume.
	// 3. record volume correpsonds to final label

	while (!AllLabeled(seed, cIter))
	{		
		std::stack<openvdb::Coord> pocketS; // for calculating volume
		std::map<openvdb::Coord, bool> marked; // for pocketCC to refuse push;
		
		//initialize marked
		for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
			marked[iter.getCoord()] = false;
		}
		// careful!!!

		// if marked by pair, record
		if (current < wantedLabels.size() && wantedLabels[current] == label)
		{
			std::cout << "Representative " << step << " " << label << std::endl;
			++current;

			SpanSeedWithWantedLabel(seed, label, true, finalLabel, pocketS, marked);
			++finalLabel;
			++_totalCCNum;

			// extract volume
			double vol = 0;
			ExtractPocketVolume(pocketS, marked, vol);
			_volume.push_back(vol);
			_depth.push_back(labelToDepth[label]);

			if (vol*pow(_voxelSize,3) < m_volumeMinimum) {
				RemovePocket(pocketS);
				--finalLabel;
				--_totalCCNum;
			}
			else {
				;
				//ExtractPocket(pocketS);
			}
		}
		else
		{
			SpanSeedWithWantedLabel(seed, label, false, 0, pocketS, marked);
		}

		label++;
	}

	_numCC = label;

	return 1;
}

int CavityDetection::SpanSeedWithWantedLabel(openvdb::Coord seed, int label, bool special, int finalLabel, 
												std::stack<openvdb::Coord>& pocketS, std::map<openvdb::Coord, bool>& marked)
{
	std::vector<openvdb::Coord> voxels;

	std::list<openvdb::Coord> grayNodes;
	grayNodes.push_back(seed);
	// Every push, mark label
	_ccLabel[seed] = label;

	int interestedCount = 0;

	if (special){
		_finalCCLabel[seed] = finalLabel;
		voxels.push_back(seed);
		pocketS.push(seed);
		marked[seed] = true;
		++interestedCount;
	}
	while (!grayNodes.empty()){
		openvdb::Coord coord = grayNodes.back();
		grayNodes.pop_back();

		// push adjacent interesting unblocked ones
		// check all direction
		for (int i = 0; i < 27; ++i){
			openvdb::Coord adjCoord = coord;
			int ii = i;

			int a, b, c;
			a = ii / 9; ii = ii - a * 9;
			b = ii / 3; ii = ii - b * 3;
			c = ii;

			if (a == 1)
				adjCoord.x()++;
			else if (a == 2)
				adjCoord.x()--;
			if (b == 1)
				adjCoord.y()++;
			else if (b == 2)
				adjCoord.y()--;
			if (c == 1)
				adjCoord.z()++;
			else if (c == 2)
				adjCoord.z()--;

			// test
			if (_proteinInterested[adjCoord] && !_proteinBlocked[adjCoord])
				if (_ccLabel[adjCoord] == -1){
					grayNodes.push_back(adjCoord);
					_ccLabel[adjCoord] = label;

					if (special){
						_finalCCLabel[adjCoord] = finalLabel;
						voxels.push_back(adjCoord);
						pocketS.push(adjCoord);
						marked[adjCoord] = true;
						++interestedCount;
					}
				}
		}
	}

	if (special){
		_area.push_back(interestedCount);
		_finalLabelToCC[finalLabel] = voxels;
	}

	return 1;
}

int CavityDetection::ExtractPocketVolume(std::stack<openvdb::Coord> pocketS, std::map<openvdb::Coord, bool>& marked, double& volume)
{
	volume = 0;

	while (!pocketS.empty()){
		openvdb::Coord coord = pocketS.top();
		pocketS.pop();
		volume += 1;

		for (int i = 0; i < 27; ++i){
			openvdb::Coord adjCoord = coord;
			int ii = i;

			int a, b, c;
			a = ii / 9; ii = ii - a * 9;
			b = ii / 3; ii = ii - b * 3;
			c = ii;

			if (a == 1)
				adjCoord.x()++;
			else if (a == 2)
				adjCoord.x()--;
			if (b == 1)
				adjCoord.y()++;
			else if (b == 2)
				adjCoord.y()--;
			if (c == 1)
				adjCoord.z()++;
			else if (c == 2)
				adjCoord.z()--;
		
			openvdb::DoubleGrid::Accessor proteinAcc = _proteinGrid->getAccessor();
			openvdb::DoubleGrid::Accessor membraneAcc = _membraneGrid->getAccessor();

			if (proteinAcc.getValue(adjCoord) > membraneAcc.getValue(adjCoord)
				&& proteinAcc.getValue(adjCoord) > 0 && membraneAcc.getValue(adjCoord) < 0
				&& !marked[adjCoord]){
					pocketS.push(adjCoord);
					marked[adjCoord] = true;
			}
		}
	}

	return 1;
}

int CavityDetection::RemovePocket(std::stack<openvdb::Coord> pocketS)
{
	//unmark what is in stack
	while (!pocketS.empty()){
		openvdb::Coord coord = pocketS.top();
		pocketS.pop();
		_finalCCLabel[coord] = -1;
	}

	_volume.pop_back();
	_area.pop_back();
	_depth.pop_back();
}

int CavityDetection::OutputInfo(double time)
{
	// output volume and adjacent atoms
	std::ofstream out("Result.txt", std::ios::out);

	out << "Execution Time: " << time << std::endl;
	out << "Total Active Voxels: " << _totalNumOfActiveVoxels << std::endl;
	out<<std::endl;

	int label = 0;
	for (int i = 0; i < _volume.size(); ++i){
		if (_volume[i] > m_volumeMinimum){
			out << "Pocket " << label << ":" << std::endl;
			out << "Area: " << _area[i]* pow(_voxelSize, 3) / (2 * _voxelSize) << std::endl;
			out << "Volume: " << _volume[i] * pow(_voxelSize, 3) << std::endl;
			out << "Max Depth: " << _depth[i] << std::endl;
			out << "Adjacent Atoms: "<< std::endl;
			int kkk = 0;
			for (std::set<size_t>::iterator iter = _CCLabelToAtoms[i].begin(); iter != _CCLabelToAtoms[i].end(); ++iter){
				out << *iter << " ";
				kkk++;
			}

			out << std::endl;
			label++;
		}
	}

	return 1;
}

int CavityDetection::FillColorVector()
{
	openvdb::Vec3d red(0.8, 0.2,0.2);
	openvdb::Vec3d green(0.2, 0.8, 0.2);
	openvdb::Vec3d blue(0.2, 0.2, 0.8);
	openvdb::Vec3d yellow(0.8, 0.8, 0.2);
	openvdb::Vec3d cyan(0.8, 0.2, 0.8);
	openvdb::Vec3d indigo(0.2, 0.8, 0.8);
	openvdb::Vec3d color1(0.8, 0.5, 0.2);
	openvdb::Vec3d color2(0.5, 0.8, 0.2);
	openvdb::Vec3d color3(0.2, 0.8, 0.5);
	openvdb::Vec3d color4(0.2, 0.5, 0.8);
	openvdb::Vec3d color5(0.8, 0.2, 0.5);
	openvdb::Vec3d color6(0.5, 0.2, 0.8);

	
	colorVector.push_back(red);
	colorVector.push_back(green);
	colorVector.push_back(blue);
	colorVector.push_back(yellow);
	colorVector.push_back(cyan);
	colorVector.push_back(indigo);
	colorVector.push_back(color1);
	colorVector.push_back(color2);
	colorVector.push_back(color3);
	colorVector.push_back(color4);
	colorVector.push_back(color5);
	colorVector.push_back(color6);

	return 1;
}

int CavityDetection::GetCenterAndDist()
{

	center = openvdb::Vec3d(0, 0, 0);
	int num = 0;
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		center += _xform->indexToWorld(iter.getCoord());
		++num;
	}

	center /= num;

	dist = -1;
	for (openvdb::DoubleGrid::ValueOnIter iter = _proteinGrid->beginValueOn(); iter.test(); ++iter){
		double l = (_xform->indexToWorld(iter.getCoord()) - center).length();
		if (l > dist)
			dist = l;
	}
	dist = dist * 2;

	return 1;
}

int CavityDetection::ExtractPocket(std::stack<openvdb::Coord> pocketS)
{
	std::cout<<"Enter ..."<< std::endl;

	double bg = _proteinGrid->background();
	
	openvdb::DoubleGrid::Ptr pocketGrid;
	pocketGrid = _proteinGrid->deepCopy();
	pocketGrid->clear();

	///////
	openvdb::DoubleGrid::Accessor pocketAcc = pocketGrid->getAccessor();
	openvdb::DoubleGrid::Accessor proteinAcc = _proteinGrid->getAccessor();
	openvdb::DoubleGrid::Accessor membraneAcc = _membraneGrid->getAccessor();
	std::map<openvdb::Coord, bool> visited;
	while (!pocketS.empty())
	{
		openvdb::Coord coord = pocketS.top();
		pocketS.pop();

		visited[coord] = true;
		if (proteinAcc.getValue(coord) > membraneAcc.getValue(coord)) {
			pocketAcc.setValue(coord, std::max(membraneAcc.getValue(coord), -proteinAcc.getValue(coord)));
		}

		for (int i = 0; i < 27; ++i)
		{
			openvdb::Coord adjCoord = coord;
			int ii = i;

			int a, b, c;
			a = ii / 9; ii = ii - a * 9;
			b = ii / 3; ii = ii - b * 3;
			c = ii;

			if (a == 1)
				adjCoord.x()++;
			else if (a == 2)
				adjCoord.x()--;
			if (b == 1)
				adjCoord.y()++;
			else if (b == 2)
				adjCoord.y()--;
			if (c == 1)
				adjCoord.z()++;
			else if (c == 2)
				adjCoord.z()--;

			
			if (proteinAcc.getValue(adjCoord) > membraneAcc.getValue(adjCoord)
				&& proteinAcc.getValue(adjCoord) > 0 && membraneAcc.getValue(adjCoord) < 0
				&& !visited[adjCoord])
				pocketS.push(adjCoord);
		}
	}

	OutPutTriangleMesh(pocketGrid);
	std::cout<<"Pocket Mesh Complete..."<<std::endl;
	int asd;
	std::cin >> asd;

	return 1;
}


int CavityDetection::ExtractPGraph()
{
	std::vector<int> edges;

	std::stack<PNodeIter> s;
	s.push(persistenceGraph->Root());

	while (!s.empty()){
		PNodeIter pni = s.top();
		s.pop();

		for (int i = 0; i < pni->children.size(); ++i) {
			edges.push_back(pni->idx);
			edges.push_back(pni->children[i]->idx);

			s.push(pni->children[i]);
		}
	}

	std::ofstream graph("graph.txt");
	for (int i = 0; i < edges.size(); i=i+2) {
		graph << edges[i]+1 << " " << edges[i + 1]+1 << std::endl;
	}
	graph << "*****************************" << std::endl;
	graph << std::endl;
	graph << std::endl;

	//output node color
	for (int i = 0; i < _components.size(); ++i) {
		for (std::list<PNodeIter>::iterator iter = _components[i].path.begin(); iter != --_components[i].path.end(); ++iter) {
			std::list<PNodeIter>::iterator iitt = iter;
			++iitt;
			graph << (*iter)->idx + 1 << " " << (*iitt)->idx + 1 << " " << i + 1 << std::endl;
		}	
	}
		

	return 1;
}

int CavityDetection::Clear()
{
	_laplacian.clear();
	_displacement.clear();

	_blocked.clear();
	_proteinBlocked.clear();
	_membraneInterested.clear();
	_ccLabel.clear(); // connected component label
	_labelToNode.clear();
	_labelToCCSize.clear();

	_preBlocked.clear();
	_preProteinBlocked.clear();
	_preMembraneInterested.clear();
	_preCCLabel.clear();

	_preLabelToNode.clear();
	_preLabelToCCSize.clear();

	_membraneVertices.clear();
	_membraneFaces.clear();

	return 1;
}

PNodeIter PGraph::InerstNode(int l, int s, bool ir, double a, double csa, openvdb::Vec3d p, int pls)
{
	PNode node;
	nodes.push_back(node);
	PNodeIter pni = --nodes.end();

	pni->label = l;
	pni->step = s;
	pni->parent = nodes.end();
	pni->isRoot = ir;
	pni->isLeaf = true;
	pni->pos = p;
	pni->approxArea = a;
	pni->csArea = csa;

	pni->votes.resize(pls);
	for (int kk = 0; kk < pni->votes.size(); ++kk)
		pni->votes[kk] = 0;

	if (s > depth)
		depth = s;

	return pni;
}

PNodeIter PGraph::InsertNode(PNodeIter& pnIter)
{
	PNode node;
	nodes.push_back(node);
	PNodeIter pni = --nodes.end();

	pni->label = pnIter->label;
	pni->step = pnIter->step;
	pni->parent = nodes.end();
	pni->isRoot = pnIter->isRoot;
	pni->isLeaf = true;
	pni->pos = pnIter->pos;
	pni->approxArea = pnIter->approxArea;
	pni->csArea = pnIter->csArea;
	pni->height = pnIter->height;

	if (pni->step > depth)
		depth = pni->step;

	return pni;
}

PNodeIter PGraph::Root()
{
	PNodeIter pni;
	for (PNodeIter iter = nodes.begin(); iter != nodes.end(); ++iter)
		if (iter->parent == nodes.end()){
			pni = iter;
			break;
		}

	return pni;
}