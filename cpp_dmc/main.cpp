#include <iostream>

// Project files
#include "UniformGrid.h"
#include "Volumes.h"
#include "DualMarchingCubes.h"

int main()
{
    using Surface = dmc::Volumes::Surface;
    using Vertex = dmc::Vector;
    using Normal = dmc::Vector;
    dmc::UniformGrid ugrid;
    dmc::Volumes vol;
    dmc::DualMarchingCubes dmc;

    // proble size
    const int nx = 64;
    const int ny = 64;
    const int nz = 64;
    // generate volume data
    vol.scalar<Surface::GenusTwo>(ugrid, nx, ny, nz);
    const double i0 = 0.0;
    // compute iso-surfaces
    std::vector<Vertex> v;
    std::vector<Normal> n;
    std::vector<int> t; // triangles
    std::vector<int> q; // quadrilaterals
    const bool sFlag{true}; // mesh simplification
    std::cout << " ... compute iso-surface" << std::endl;
    dmc.dualMC(i0, ugrid, v, n, t, q, sFlag);

    // write mesh in obj file format
    std::cout << " ... writing obj\n";
    std::ofstream objF;
    // write triangle mesh
    objF.open("dmc_mesh_tris.obj");
    if (!objF.is_open()) {
        std::cout << "ERROR: can't open output file " << std::endl;
    }
    else {
        objF << "#Dual Marching Cubes, triangle mesh\n";
        const int nr_v = v.size();
        for (int i = 0; i < nr_v; i++)
        {
            objF << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
        }
        for (int i = 0; i < nr_v; i++)
        {
            objF << "vn " << n[i][0] << " " << n[i][1] << " " << n[i][2] << std::endl;
        }
        const int nr_t = t.size() / 3;
        for (int i = 0; i < nr_t; i++)
        {

            objF << "f " << (t[3 * i] + 1) << "//" << (t[3 * i] + 1)
            << " " << (t[3 * i + 1] + 1) << "//" << (t[3 * i + 1] + 1)
            << " " << (t[3 * i + 2] + 1) << "//" << (t[3 * i + 2] + 1) << std::endl;
        }
        objF.close();
    }
    // write quad mesh
    objF.open("dmc_mesh_quad.obj");
    if (!objF.is_open()) {
        std::cout << "ERROR: can't open output file " << std::endl;
    }
    else {
        objF << "#Dual Marching Cubes, quad only mesh\n";
        const int nr_v = v.size();
        for (int i = 0; i < nr_v; i++)
        {
            objF << "v " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
        }
        for (int i = 0; i < nr_v; i++)
        {
            objF << "vn " << n[i][0] << " " << n[i][1] << " " << n[i][2] << std::endl;
        }
        const int nr_q = q.size() / 4;
        for (int i = 0; i < nr_q; i++)
        {
            objF << "f " << (q[4 * i] + 1) << "//" << (q[4 * i] + 1)
            << " " << (q[4 * i + 1] + 1) << "//" << (q[4 * i + 1] + 1)
            << " " << (q[4 * i + 2] + 1) << "//" << (q[4 * i + 2] + 1)
            << " " << (q[4 * i + 3] + 1) << "//" << (q[4 * i + 3] + 1) << std::endl;
        }
        objF.close();
    }

    return 0;
}
