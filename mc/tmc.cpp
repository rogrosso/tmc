//
//

// Libs
#include <iostream>
#include <string>

// Project files
#include "MarchingCubes.h"

int main(int argc, char* argv[])
{
    if (argc != 4) {
        std::cout << "usage: " << std::endl;
        std::cout << "       tmc input_filename  output_obj_filename output_off_filename \n";
        exit(1);
    }
    tmc::MarchingCubes mc;
    std::string i_file(argv[1]);
    std::string objF(argv[2]);
    std::string offF(argv[3]);
    mc(i_file,objF,offF);
    return 0;
}


