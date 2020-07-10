
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <array>
#include <map>

// project
#include "Mesh.h"
#include "DualMarchingCubes.h"

int main()
{
    using SurfaceCase = p_mc::UniformGrid::SurfaceCase;
    // Implicit Surfaces
    //  GenusTwo = 4,
    //  iWP = 5, (triple symmetry)
    //  pwHybrid = 6,
    //  neovius = 7,
    //  Goursat = 8,
    //  SteinerRoman = 9
    std::map<std::string,int> config;
    config["implicit"] = 5;
    config["p3X3YColor"] = 1; // color based mesh simplification
    config["p3X3YOld"] = 1;  // simplify isolated elements, i.e. no neighbor with same valence pattern
    config["p3333"] = 1; // simplify vertex valence pattern 3333
    std::array<int, 3> dim{128,128,128};
    p_mc::DualMarchingCubes dmc;
    p_mc::Mesh mesh;
    
    // Init: generate a synthetic data set
    std::cout << " ... init" << std::endl;
    dmc.init(dim,static_cast<SurfaceCase>(7));
    const float i0 = 0.0f;

    // compute iso-surface
    std::cout << " ... compute iso-surface" << std::endl;
    dmc.dualMC(i0, mesh, config);
    // write mesh in obj file format
    std::cout << " ... writing obj\n";
    std::string o_file{"dmc_mesh.obj"};
    // write triangle mesh
    //mesh.writeObjTriangles(o_file);
    // write quad mesh
    mesh.writeObjQuads(o_file);

    return 0;
}


