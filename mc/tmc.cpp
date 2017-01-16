//
//

// Libs
#include <iostream>
#include <string>

// Project files
#include "MarchingCubes.h"

int main(int argc, char* argv[])
{
	// Marcing cubes
    tmc::MarchingCubes mc;

	// parameters
	std::string i_file = "";
	std::string objF = "";
	std::string offF = "";
    const double i0 = 801.3;
    mc(i0,i_file,true,true,objF,true,offF);
    return 0;
}
