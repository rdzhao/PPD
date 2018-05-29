#include "CavityDetection.h"
//#include <chrono>

int main(int argc, char* argv[])
{

	CavityDetection cd;

	std::stringstream ssDepth(argv[3]);
	std::stringstream ssArea(argv[4]);
	std::stringstream ssVolume(argv[5]);
	double depth, area, volume;
	ssDepth >> depth;
	ssArea >> area;
	ssVolume >> volume;


	//std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	cd.InitVDB(200, 1.0);
	cd.InitMoleculerSurface(argv[1]);
	cd.InitMembrane();
	cd.setParameters(depth, area, volume);
	cd.countActiveVoxels();

	cd.FirstPass();
	cd.SecondPass();

	std::cout << "Start collecting atom info..."<< std::endl;
	cd.PocketAtomMap(argv[2]);

	//std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	//std::chrono::duration<double> elapsed = t2 - t1;
	//std::cout << "Execution Time: " << elapsed.count() << std::endl;
	//cd.OutputInfo(elapsed.count());
	cd.OutputInfo(100);

	


	//cd.View();

	return 1;
}