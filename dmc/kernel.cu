
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <vector>
#include <array>
#include <map>

// project
#include "UniformGrid.h"
#include "DualMarchingCubes.h"

int main()
{
    //  List of possible surfaces generated implicitly
    //  GenusTwo
    //  iWP (triple symmetry)
    //  pwHybrid
    //  neovius
    //  Goursat
    //  SteinerRoman
    using SurfaceCase = p_mc::UniformGrid::SurfaceCase;

    std::map<std::string,int> config;
    config["valence"] = 0; // computes vertex valences
    config["element-quality"] = 0; // computes element quality
    config["p3X3YColor"] = 1; // color based mesh simplification
    config["p3X3YOld"] = 1;  // simplify isolated elements, i.e. no neighbor with same valence pattern
    config["p3333"] = 1; // simplify vertex valence pattern 3333
    config["halfedge-datastructure"] = 0; // computes halfedge data structure for quad only mesh
    config["non-manifold"] = 0; // compute number of non-manifold edges
    std::array<int, 3> dim{64,64,64};
    p_mc::DualMarchingCubes dmc;

    // Init: generate a synthetic data set
    std::cout << " ... init" << std::endl;
    dmc.init(dim,SurfaceCase::GenusTwo);
    const float i0 = 0.0f;

    // mesh
    using Vertex = p_mc::DualMarchingCubes::Vertex;
    using Normal = p_mc::DualMarchingCubes::Normal;
    using Halfedge = p_mc::DualMarchingCubes::Halfedge;
    using HalfedgeFace = p_mc::DualMarchingCubes::HalfedgeFace;
    using HalfedgeVertex = p_mc::DualMarchingCubes::HalfedgeVertex;
    using Triangle = p_mc::DualMarchingCubes::Triangle;
    using Quadrilateral = p_mc::DualMarchingCubes::Quadrilateral;

    std::vector<Vertex> v;
    std::vector<Normal> n;
    std::vector<Triangle> t;
    std::vector<Quadrilateral> q;
    std::vector<Halfedge> h;
    std::vector<HalfedgeFace> hf;
    std::vector<HalfedgeVertex> hv;

    // compute iso-surface
    std::cout << " ... compute iso-surface" << std::endl;
    dmc.dualMC(i0, v, n, t, q, h, hf, hv, config);

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
        const int nr_t = t.size();
        for (int i = 0; i < nr_t; i++)
        {
            objF << "f " << (t[i][0] + 1) << "//" << (t[i][0] + 1)
            << " " << (t[i][1] + 1) << "//" << (t[i][1] + 1)
            << " " << (t[i][2] + 1) << "//" << (t[i][2] + 1) << std::endl;
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
        const int nr_q = q.size();
        for (int i = 0; i < nr_q; i++)
        {
            objF << "f " << (q[i][0] + 1) << "//" << (q[i][0] + 1)
            << " " << (q[i][1] + 1) << "//" << (q[i][1] + 1)
            << " " << (q[i][2] + 1) << "//" << (q[i][2] + 1)
            << " " << (q[i][3] + 1) << "//" << (q[i][3] + 1) << std::endl;
        }
        objF.close();
    }
    return 0;
}
