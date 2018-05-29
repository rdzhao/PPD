#pragma once

#ifndef _CAVITY_DETECTION_H_ 
#define _CAVITY_DETECTION_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <ctime>

#include <omp.h>

#include <openvdb/openvdb.h>
#include <openvdb/Grid.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/VolumeToMesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "nanoflann.hpp"
#include "FlannUtils.h"

//#include "Viewer.h"

using namespace nanoflann;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;
typedef Polyhedron_3::HalfedgeDS HalfedgeDS;

typedef Polyhedron_3::Vertex_iterator Vertex_iterator;
typedef Polyhedron_3::Facet_iterator Facet_iterator;

class PGraph;
class PNode;
class PNICompare;
class HomologyComponent;
class HCCompareDescend;
class CavityDetection;

typedef std::list<PNode>::iterator PNodeIter;

class CavityDetection
{
public:
	int InitVDB(int dss, double vs);
	int InitMoleculerSurface(std::string msFile);
	int InitMembrane();
	
	int setParameters(double d = 5, double a = 100, double v = 200);
	int countActiveVoxels();

	int FirstPass();
	int SecondPass();
	
	int Deform();

	int DeformToMarkRepresentativePockets(); 

	int PocketAtomMap(std::string xyzrFile);

	//int View(openvdb_viewer::Viewer& viewer, int counter);
	//int ViewPocket(openvdb_viewer::Viewer& pViewer, openvdb::DoubleGrid::Ptr pocket);

	int GetBoundingBox();

	int OutputInfo(double time);

private:
	int ReadMolecularSurfaceTriangleMesh(std::string msFile);

	int MeanCurvatureDeformation(double& tStep, double& forceRatio, int& i);
	int ReinitializeLevelSet();

	// obj file parser
	int ReadMeshData(std::istream& in);
	int ReadPosition(std::stringstream& ss);
	int ReadFace(std::stringstream& ss);
	int ParseFaceIndex(const std::string& token);
	int MoveToOrigin();
	// obj file parser end

	int GetConvexHull();
	int BuildAndRelaxConvexHull();

	int SelectProteinInterested(int& tav);
	int SelectMembraneInterested();
	int InitProteinBlocked();
	bool AllProteinBlocked();
	bool AllBlocked();

	int OutPutTriangleMesh(openvdb::DoubleGrid::Ptr grid);

	// for persistence graph
	int LabelConnectedComponents();
	int InitLabel();
	bool AllLabeled(openvdb::Coord& seed, openvdb::DoubleGrid::ValueOnIter& cIter);
	int SpanSeed(openvdb::Coord seed, int label, openvdb::Vec3d& pos, int& voxelNum, int& csVoxelNum, bool test = false);
	int MoveInfo();
	int AppendPersistenceGraph(int step);
	int AssembleChildren();
	int CheckPersistenceGraphCompleteness();
	int EliminateNoise(); // some part of protein interested voxels will freeze, we should eliminate it.
	int RemovePNode(PNodeIter pni);
	int EliminateSmallTree();
	int RemoveTree(PNodeIter root);
	int AssemblePGraphData();

	// for component extraction
	int AssembleHomologyComponent(int numDeform); 
	int TraceBack(PNodeIter seed, HomologyComponent& hc);
	int EliminateShortPersistence(double r); // the same as eliminate small depth persistence
	int GetNodeHeight();
	double RecursiveHeight(PNodeIter pni);

	int ReconstructMajorComponentTree();
	int ExtractRepresentatives(); //just mark protein surface components (from coord to label)
	int TraceBranchUntilBifurcation(PNodeIter& pni);

	int InitFinalCCLabel();
	int SelectStepLabelPair(std::vector<int>& wantedLabels, std::map<int, double>& labelToDepth, int step, std::vector<PNodeIter> reps);
	int LabelCCAndExtractReps(std::vector<int> wantedLabels, std::map<int, double> labelToDepth,int step, int& finalLabel);
	int SpanSeedWithWantedLabel(openvdb::Coord seed, int label, bool special, int finalLabel, std::stack<openvdb::Coord>& pocketS, std::map<openvdb::Coord, bool>& marked);
	int ExtractPocketVolume(std::stack<openvdb::Coord> pocketS, std::map<openvdb::Coord, bool>& marked, double& volume);
	int RemovePocket(std::stack<openvdb::Coord> pocketS);

	int FillColorVector();
	int GetCenterAndDist();

	// pocket bound volume
	int ExtractPocket(std::stack<openvdb::Coord> pocketS);

	// graph
	int ExtractPGraph();

	int Clear();

	double _voxelSize;
	double _halfWidth;
	int deformStepSize;
	openvdb::DoubleGrid::Ptr _proteinGrid;
	openvdb::DoubleGrid::Ptr _OriMembraneGrid;
	openvdb::DoubleGrid::Ptr _membraneGrid;
	int _totalNumOfActiveVoxels;
	int _proteinNumOfBoundaryVoxels;

	std::map<openvdb::Coord, double> _laplacian;
	std::map<openvdb::Coord, double> _displacement;

	std::map<openvdb::Coord, bool> _proteinInterested;

	std::map<openvdb::Coord, bool> _blocked;	
	std::map<openvdb::Coord, bool> _proteinBlocked;
	std::map<openvdb::Coord, bool> _membraneInterested;
	std::map<openvdb::Coord, int> _ccLabel; // connected component label
	std::map<openvdb::Coord, int> _ccLabelMembrane; // connected component label for membrane
	int _numCC;
	std::map<int, PNodeIter> _labelToNode;
	std::map<int, int> _labelToCCSize;
	std::map<int, int> _labelToCCSizeMembrane; // for membrane
	
	std::map<openvdb::Coord, bool> _preBlocked;
	std::map<openvdb::Coord, bool> _preProteinBlocked;
	std::map<openvdb::Coord, bool> _preMembraneInterested;
	std::map<openvdb::Coord, int> _preCCLabel;
	int _preNumCC;
	std::map<int, PNodeIter> _preLabelToNode;
	std::map<int, int> _preLabelToCCSize;

	std::map<int, std::vector<int> > _voteForParent;

	std::vector<openvdb::Vec3d> colorVector;

	std::vector<openvdb::Vec3s> _meshVertices;
	std::vector<openvdb::Vec3I> _meshFaces;

	std::vector<openvdb::Vec3s> _membraneVertices;
	std::vector<openvdb::Vec3I> _membraneFaces;

	openvdb::math::Transform::Ptr _xform;
	
	PGraph* persistenceGraph;
	std::map<int, openvdb::Vec3d> CCPos;
	std::vector<openvdb::Vec3d> gLocations;
	std::vector<int> gIndices;
	std::vector<double> gSteps; // re-range to [0,1]

	std::vector<HomologyComponent> _components;
	PGraph* repPGraph;
	std::vector<PNodeIter> reps;
	std::map<openvdb::Coord, int> _finalCCLabel;
	std::map<int, std::vector<openvdb::Coord> > _finalLabelToCC;
	int _totalCCNum;
	std::vector<double> _volume;
	std::vector<double> _area;
	std::vector<double> _depth;

	std::vector<openvdb::Vec3d> _atomLocations;
	std::map<int, std::set<size_t> > _CCLabelToAtoms;

	openvdb::Vec3d center;
	double dist;

	int _numBVM; // # boundary voxels of membrane

	// all related parameters
	double m_persistenceThreshold;
	double m_areaMinimum;
	double m_depthMinimum;
	double m_volumeMinimum;
	double m_moveSpeed;
};

class PGraph
{
public:
	PGraph() { depth = 0; };
	PNodeIter InerstNode(int l, int s, bool ir, double a, double csa, openvdb::Vec3d p, int pls);
	PNodeIter InsertNode(PNodeIter& pnIter);
	
	PNodeIter Root();

	std::list<PNode> nodes;
	int depth;
};

class PNode
{
public:
	PNode() {};

	inline bool operator< (PNode& rhs) { return (step < rhs.step); };

	int idx;
	bool isRoot;
	bool isLeaf;
	bool leafTreated;
	int label;
	int step;
	bool passed;
	double approxArea;
	double csArea; // cross section area
	double depth;
	double height;

	PNodeIter parent;
	std::vector<PNodeIter> children;
	std::vector<int> votes;

	openvdb::Vec3d pos;
};

class PNICompare
{
public:
	bool operator()(const PNodeIter& lhs, const PNodeIter& rhs) const
	{
		bool less;
		if (lhs->step < rhs->step)
			less= true;
		else if (lhs->step > rhs->step)
			less=false;
		else{
			less = (lhs->label < rhs->label);
		}

		return less;
		//return lhs->step < rhs->step;
	}

};

class HomologyComponent
{
public:
	HomologyComponent() {};
	HomologyComponent(const HomologyComponent& hc) : path(hc.path), bornTime(hc.bornTime), killTime(hc.killTime) {};

	inline int lifeTime() const { return killTime - bornTime; };

	std::list<PNodeIter> path;
	int bornTime;
	int killTime;
};

class HCCompareDescend
{
public:
	bool operator()(const HomologyComponent& lhs, const HomologyComponent& rhs) const{
		return lhs.lifeTime() > rhs.lifeTime();
	}
};

#endif // !_CAVITY_DETECTION_H_ 
