
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <array>

// project
#include "Mesh.h"
#include "DualMarchingCubes.h"

int main()
{
    using ushort = unsigned short;
    p_dmc::DualMarchingCubes dmc;
    p_dmc::Mesh mesh;
    
    // Read volume data from file
    //std::string i_file = "./head_ushort_512_512_641.bin";
    //std::string o_file = "./head_ushort_512_512_641.obj";
    //float i0 = 900; // iso-value
    //dmc.readDataFromFile<ushort>(i_file);
    // dmc.dualMC(i0, mesh);
    
    // generate synthetic volume
    std::string o_file = "./synthetic.obj";
    std::array<int, 3> dim{ 64, 64, 64 };
    float i0 = 0;
    dmc.generateData(dim);

    // compute iso-surface
    dmc.dualMC(i0, mesh);
    // write mesh in obj file format
    std::cout << " ... writing obj\n";
    // write triangle mesh
    //mesh.writeObjTriangles(o_file);
    // write quad mesh
    mesh.writeObjQuads(o_file);

    return 0;
}


