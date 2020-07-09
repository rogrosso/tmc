#include "DualMarchingCubes.h"

// CUDA
#include "helper_cuda.h"
#include "CTimer.h"
#include "QualityMeasure.h"
#include "QuadrilateralHashTable.h"
#include "HalfedgeHashTable.h"
#include "EdgeHashTable.h"
#include "VertexHashTable.h"
#include "Vertices.h"
#include "Triangles.h"
#include "Quadrilaterals.h"
#include "Edges.h"
#include "Halfedges.h"
#include "HalfedgeVertices.h"
#include "HalfedgeFaces.h"
#include "MarchingCubesLookupTables.h"
#include "CellIntersection.h"
#include "VertexMap.h"
#include "QuadrilateralMap.h"
#include "HalfedgeMesh.h"
#include "MeshSimplification.h"
#include "FaceColoring.h"

// type aliases
// Introduce convenient aliases here
using namespace p_mc;
using uint = unsigned int;
using uchar = unsigned char;
using ushort = unsigned short;
using ullong = unsigned long long;
using UGrid = p_mc::DualMarchingCubes::UGrid;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      CUDA GLOBAL FUNCTIONS
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      Compute Element Quality and Generate best triangle mesh out of a quadrilateral mesh
//      Use the MaxMin angle criterium
//      Compute triangle angles based on cosine rule
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void quadrilateral_to_triangle(Quadrilaterals q_, Vertices v_, Triangles t_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q)
        return;
    // get vertices
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    const float3 p0 = v_.vertices[v0];
    const float3 p1 = v_.vertices[v1];
    const float3 p2 = v_.vertices[v2];
    const float3 p3 = v_.vertices[v3];

    float a1_ = t_.minAngle(p0, p1, p2);
    float a2_ = t_.minAngle(p0, p2, p3);
    float b1_ = fminf(a1_, a2_);
    float b2_ = fmaxf(a1_, a2_);
    a1_ = t_.minAngle(p1, p3, p0);
    a2_ = t_.minAngle(p1, p2, p3);
    float c1_ = fminf(a1_, a2_);
    float c2_ = fmaxf(a1_, a2_);

    if (b1_ < c1_ || (b1_ == c1_ && b2_ <= c2_))
    {
        t_.addTriangle(2 * tid, v1, v3, v0);
        t_.addTriangle(2 * tid + 1, v1, v2, v3);
    }
    else 
    {
        t_.addTriangle(2 * tid, v0, v1, v2);
        t_.addTriangle(2 * tid + 1, v0, v2, v3);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute quality measure for a triangle mesh
__global__ void mean_ratio_measure(Triangles t_, Vertices v_, QualityMeasure q_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= t_.nr_t)
        return;
    // collect triangle vertices
    const int v0 = t_.triangles[tid].x;
    const int v1 = t_.triangles[tid].y;
    const int v2 = t_.triangles[tid].z;
    const float3 p0 = v_.vertices[v0];
    const float3 p1 = v_.vertices[v1];
    const float3 p2 = v_.vertices[v2];
    // compute quality measure
    q_.mean_ratio(tid, p0, p1, p2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for quadrilaterals
__global__ void init_quadrilateral_hashtable(QuadrilateralHashTable ht_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= ht_.t_size) return;
    ht_.init(tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for vertices, this is a redundant data structure, just to differentiate
__global__ void init_vertex_hashtable(VertexHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size) return;
    ht_.init(tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for halfedges
//__global__ void init_halfedge_hashtable(HalfedgeHashTable ht_)
//{
//    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
//    if (tid >= ht_.t_size) return;
//    ht_.init(tid);
//}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      DUAL MARCHING CUBES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hybrid version
//  - ambiguous cases are processed without lookup table
//  - unambiguous case are processed with the standard mc
__global__ void dual_mc(const float i0, const uint t_size, UGrid ugrid, CellIntersection c_, QuadrilateralHashTable ht_, Vertices v_)
{
    // use a 1d grid
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (t_size <= gl_index)
        return;

    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    if (i_index >= (ugrid.idim - 1) || j_index >= (ugrid.jdim - 1) || k_index >= (ugrid.kdim - 1))
    {
        return;
    }
    
    // scalar values at vertices
    float u[8];
    u[0] = ugrid(i_index, j_index, k_index);
    u[1] = ugrid(i_index + 1, j_index, k_index);
    u[2] = ugrid(i_index, j_index + 1, k_index);
    u[3] = ugrid(i_index + 1, j_index + 1, k_index);
    u[4] = ugrid(i_index, j_index, k_index + 1);
    u[5] = ugrid(i_index + 1, j_index, k_index + 1);
    u[6] = ugrid(i_index, j_index + 1, k_index + 1);
    u[7] = ugrid(i_index + 1, j_index + 1, k_index + 1);

    // 
    uchar i_case{ 0 };
    i_case = i_case + ((uint)(u[0] >= i0));
    i_case = i_case + ((uint)(u[1] >= i0)) * 2;
    i_case = i_case + ((uint)(u[2] >= i0)) * 4;
    i_case = i_case + ((uint)(u[3] >= i0)) * 8;
    i_case = i_case + ((uint)(u[4] >= i0)) * 16;
    i_case = i_case + ((uint)(u[5] >= i0)) * 32;
    i_case = i_case + ((uint)(u[6] >= i0)) * 64;
    i_case = i_case + ((uint)(u[7] >= i0)) * 128;

    if (i_case == 0 || i_case == 255)
        return;
    // intersect cell
    //c_.slice(i0, i_case, i_index, j_index, k_index, u, ugrid, ht_, v_);
    //c_.sliceQ(i0, i_case, i_index, j_index, k_index, u, ugrid, ht_, v_);
    c_.sliceP(i0, i_case, i_index, j_index, k_index, u, ugrid, ht_, v_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      Standard MARCHING CUBES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// standard marching cubes
__global__ void standard_mc(const float i0, const uint t_size, UGrid ugrid, MarchingCubesLookupTables l_tables, VertexHashTable ht_, Vertices v_, Triangles t_)
{
    // use a 1d grid
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (gl_index >= t_size)
        return;

    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    if (i_index >= (ugrid.idim - 1) || j_index >= (ugrid.jdim - 1) || k_index >= (ugrid.kdim - 1))
    {
        return;
    }

    // scalar values at vertices
    float u[8];
    u[0] = ugrid(i_index, j_index, k_index);
    u[1] = ugrid(i_index + 1, j_index, k_index);
    u[2] = ugrid(i_index, j_index + 1, k_index);
    u[3] = ugrid(i_index + 1, j_index + 1, k_index);
    u[4] = ugrid(i_index, j_index, k_index + 1);
    u[5] = ugrid(i_index + 1, j_index, k_index + 1);
    u[6] = ugrid(i_index, j_index + 1, k_index + 1);
    u[7] = ugrid(i_index + 1, j_index + 1, k_index + 1);

    // compute case
    uchar i_case{ 0 };
    i_case = i_case + ((uint)(u[0] >= i0));
    i_case = i_case + ((uint)(u[1] >= i0)) * 2;
    i_case = i_case + ((uint)(u[2] >= i0)) * 4;
    i_case = i_case + ((uint)(u[3] >= i0)) * 8;
    i_case = i_case + ((uint)(u[4] >= i0)) * 16;
    i_case = i_case + ((uint)(u[5] >= i0)) * 32;
    i_case = i_case + ((uint)(u[6] >= i0)) * 64;
    i_case = i_case + ((uint)(u[7] >= i0)) * 128;

    if (i_case == 0 || i_case == 255)
        return;
    
    // compute cell intersection
    // compute unique edge global index
    auto e_glIndex = [](const int e, const int i_idx, const int j_idx, const int k_idx, UGrid& ugrid)
    {
        const unsigned long long gei_pattern_ = 670526590282893600ull;
        const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
        const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
        const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
        const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
        return (3 * ugrid.gl_index(i, j, k) + offs);
    };
    // table listing end vertices of an edge
    const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
    const ushort e_ = l_tables.e_pattern[i_case];
    float3 n[8];
    ugrid.gradient(n, u, i_index, j_index, k_index);
    ushort flag{ 1 };
    int v_addr[12];
    for (int e = 0; e < 12; e++)
    {
        if (flag & e_)
        {
            const int e_id = e_glIndex(e, i_index, j_index, k_index, ugrid);
            // check if vertex was already generated
            const bool v_flag = ht_.addVertex(e_id, v_addr, e);
            if (!v_flag) {
                // create vertex 
                const int v0 = (l_edges_[e] & 0xF);
                const int v1 = (l_edges_[e] >> 4) & 0xF;
                const float l = (i0 - u[v0]) / (u[v1] - u[v0]);
                const float3 vt0{ ugrid.x0 + (i_index + (v0 & 0x1)) * ugrid.dx, ugrid.y0 + (j_index + ((v0 & 0x2) >> 1)) * ugrid.dy, ugrid.z0 + (k_index + ((v0 & 0x4) >> 2)) * ugrid.dz };
                const float3 vt1{ ugrid.x0 + (i_index + (v1 & 0x1)) * ugrid.dx, ugrid.y0 + (j_index + ((v1 & 0x2) >> 1)) * ugrid.dy, ugrid.z0 + (k_index + ((v1 & 0x4) >> 2)) * ugrid.dz };
                float3 ei = make_float3(0, 0, 0);
                float3 ni = make_float3(0, 0, 0);
                ei.x = vt0.x + l * (vt1.x - vt0.x);
                ei.y = vt0.y + l * (vt1.y - vt0.y);
                ei.z = vt0.z + l * (vt1.z - vt0.z);
                ni.x = n[v0].x + l * (n[v1].x - n[v0].x);
                ni.y = n[v0].y + l * (n[v1].y - n[v0].y);
                ni.z = n[v0].z + l * (n[v1].z - n[v0].z);
                const float factor = sqrtf(ni.x * ni.x + ni.y * ni.y + ni.z * ni.z);
                ni.x /= factor;
                ni.y /= factor;
                ni.z /= factor;
                // create vertex
                const int addr_ = v_.addVertex(ei, ni);
                // map index
                ht_.set(v_addr[e], addr_);
            }
        }
        flag <<= 1;
    }
    // construct triangles
    char* t_pattern = &l_tables.t_pattern[ 16 * i_case ];
    for (int t = 0; t < 16; t += 3)
    {
        if (t_pattern[t] == -1)
            break;
        // add triangle to list
        const int v0 = v_addr[static_cast<int>(t_pattern[t])];
        const int v1 = v_addr[static_cast<int>(t_pattern[t+1])];
        const int v2 = v_addr[static_cast<int>(t_pattern[t+2])];
        t_.addTriangle(v0, v1, v2);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      SHARED VERTEX or INDEXED FACE DATA STRUCTURES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id to position in vertex array
// Construct shared vertex list of triangles
__global__ void map_triangles(VertexHashTable ht_, Triangles t_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= t_.nr_t) return;
    const int v0 = t_.triangles[tid].x;
    const int v1 = t_.triangles[tid].y;
    const int v2 = t_.triangles[tid].z;

    t_.triangles[tid].x = ht_.v(v0);
    t_.triangles[tid].y = ht_.v(v1);
    t_.triangles[tid].z = ht_.v(v2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id to position in vertex array
// Construct shared vertex list for quadrilaterals
__global__ void map_quadrilaterals(QuadrilateralHashTable ht_, Quadrilaterals q_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= ht_.size()) return;
    // empty bucket
    if (ht_.empty(tid)) return;
    // add quadrilateral
    //q_.addColoredQuadrilateral(ht_.v0(tid), ht_.v1(tid), ht_.v2(tid), ht_.v3(tid), ht_.color(tid));
    q_.addColoredQuadrilateral(ht_.quadrilateral(tid), ht_.color(tid));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// color elements, inherit coloring from uniform grid, a coloring with 5 colors is always possible
// This method works, because the quad mesh has already a consistent coloring inherited from the uniform grid,
// that is, neighbors have a color different from c. 
__global__ void color_quadrilaterals(Quadrilaterals q_, Halfedges he_, HalfedgeFaces he_f, int c)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (q_.getColor(tid) != c) return; // we are processing color c
    // collect halfedge from face
    const int e0 = he_f.he_e[tid];
    const int e1 = he_.getNext(e0); // he_.he_e[e0].z;
    const int e2 = he_.getNext(e1); // he_.he_e[e1].z;
    const int e3 = he_.getNext(e2); // he_e[e2].z;
    // collect twins
    const int t0 = he_.getTwin(e0); // he_.he_e[e0].w;
    const int t1 = he_.getTwin(e1); // he_.he_e[e1].w;
    const int t2 = he_.getTwin(e2); // he_.he_e[e2].w;
    const int t3 = he_.getTwin(e3); // he_.he_e[e3].w;
    // collect neighboor faces
    int f0{ -1 };
    int f1{ -1 };
    int f2{ -1 };
    int f3{ -1 };
    if (t0 > -1) f0 = he_.getFace(t0); // he_.he_e[t0].y;
    if (t1 > -1) f1 = he_.getFace(t1); // he_.he_e[t1].y;
    if (t2 > -1) f2 = he_.getFace(t2); // he_.he_e[t2].y;
    if (t3 > -1) f3 = he_.getFace(t3); // he_.he_e[t3].y;
    // collect colors
    int c0{ -1 }; if (f0 > -1) c0 = q_.getColor(f0);
    int c1{ -1 }; if (f1 > -1) c1 = q_.getColor(f1);
    int c2{ -1 }; if (f2 > -1) c2 = q_.getColor(f2);
    int c3{ -1 }; if (f3 > -1) c3 = q_.getColor(f3);
    // check which is the first free color from the first 5 colors
    if (c0 != 0 && c1 != 0 && c2 != 0 && c3 != 0)
    {
        q_.setColor(tid, 0);
        return;
    }
    if (c0 != 1 && c1 != 1 && c2 != 1 && c3 != 1)
    {
        q_.setColor(tid, 1);
        return;
    }
    if (c0 != 2 && c1 != 2 && c2 != 2 && c3 != 2)
    {
        q_.setColor(tid, 2);
        return;
    }
    if (c0 != 3 && c1 != 3 && c2 != 3 && c3 != 3)
    {
        q_.setColor(tid, 3);
        return;
    }
    if (c0 != 4 && c1 != 4 && c2 != 4 && c2 != 4)
    {
        q_.setColor(tid, 4);
        return;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// optimize element colors by removing the fifth color, if possible
// remember halfedge data structure
//    he.x = origin vertex
//    he.y = face
//    he.z = next
//    he.w = twin
__global__ void optimize_coloring(Quadrilaterals q_, Halfedges he_, HalfedgeFaces he_f)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (q_.getColor(tid) != 4) return; // we are processing color c
    // collect halfedge from face
    const int e0 = he_f.he_e[tid];
    const int e1 = he_.he_e[e0].z;
    const int e2 = he_.he_e[e1].z;
    const int e3 = he_.he_e[e2].z;
    // collect twins
    const int t0 = he_.he_e[e0].w;
    const int t1 = he_.he_e[e1].w;
    const int t2 = he_.he_e[e2].w;
    const int t3 = he_.he_e[e3].w;
    // collect neighboor faces
    int f0{ -1 };
    int f1{ -1 };
    int f2{ -1 };
    int f3{ -1 };
    if (t0 > -1) f0 = he_.he_e[t0].y;
    if (t1 > -1) f1 = he_.he_e[t1].y;
    if (t2 > -1) f2 = he_.he_e[t2].y;
    if (t3 > -1) f3 = he_.he_e[t3].y;
    // collect colors
    int c0{ -1 }; if (f0 > -1) c0 = q_.getColor(f0);
    int c1{ -1 }; if (f1 > -1) c1 = q_.getColor(f1);
    int c2{ -1 }; if (f2 > -1) c2 = q_.getColor(f2);
    int c3{ -1 }; if (f3 > -1) c3 = q_.getColor(f3);
    // check which is the first free color from the first 5 colors
    if (c0 != 0 && c1 != 0 && c2 != 0 && c3 != 0)
    {
        q_.setColor(tid, 0);
        return;
    }
    if (c0 != 1 && c1 != 1 && c2 != 1 && c3 != 1)
    {
        q_.setColor(tid, 1);
        return;
    }
    if (c0 != 2 && c1 != 2 && c2 != 2 && c3 != 2)
    {
        q_.setColor(tid, 2);
        return;
    }
    if (c0 != 3 && c1 != 3 && c2 != 3 && c3 != 3)
    {
        q_.setColor(tid, 3);
        return;
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// transfer attributes from quadrilaterals to halfedge faces
__global__ void transfer_face_attributes(Quadrilaterals q_, HalfedgeFaces he_f)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    he_f.attributes[tid] = q_.attributes[tid];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      HOST CODE
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      Dual Marching Cubes -- Host code
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// compute valence
void checkVertexValence(p_mc::Vertices v, p_mc::Quadrilaterals q, std::string const& step)
{
    /// data
    p_mc::HalfedgeMesh hm;
    p_mc::Edges e;
    p_mc::EdgeHashTable eht;
    /// computes edges
    hm.edgeHashTable(q, eht);
    std::vector<int2> le;
    std::vector<int> nrF;
    eht.getEdges(le);
    eht.getNrFaces(nrF);
    
    std::vector<int> valence(v.size(), 0);
    for (int i = 0; i < le.size(); i++)
    {   
       if (nrF[i] > 0)
        {
            valence[le[i].x] += 1;
            valence[le[i].y] += 1;
        }

    }
    // assume largest valence is 11
    std::vector<int> valenceDist(12, 0);
    for (auto v : valence)
    {
        if (v <= 11)
        {
            if (v == 0 || v == 1)
            {
                std::cout << "ERROR: valence : " << v << std::endl;
            }
            valenceDist[v] += 1;
        }
        else {
            std::cout << "ERROR: valence larger than 11, valence: " << v << std::endl;
        }
    }
    // write valence distribution to text file, accumulate
    std::string o_file = "./data/measurements/valencedistirbution.txt";
    std::ofstream of;
    of.open(o_file, std::ios_base::app);
    of << "new measure: " << step << std::endl;
    for (auto v : valenceDist)
    {
        of << v << std::endl;
    }
    of.close();
}


/// CPU function to control mesh structure
void checkMeshConsistency(p_mc::Vertices v, p_mc::Quadrilaterals q, std::string const &step)
{
    p_mc::HalfedgeMesh hm;
    p_mc::HalfedgeFaces hef;
    p_mc::HalfedgeVertices hev;
    p_mc::Halfedges he;
    p_mc::Edges e;
    p_mc::EdgeHashTable eht;

    /// computes edges
    hm.edgeHashTable(q, eht);
    std::vector<int2> le;
    std::vector<int> nrF;
    eht.getEdges(le);
    eht.getNrFaces(nrF);
    // count nr. of boundary and non-manifold edges
    int nrBndEdges{ 0 };
    int nrInnerEdges{ 0 };
    int nrNonManifold{ 0 };
    int nrWrongEdges{ 0 };
    for (int i = 0; i < le.size(); i++)
    {
        switch (nrF[i])
        {
        case 0:
            break;
        case 1:
            nrBndEdges++;
            break;
        case 2:
            nrInnerEdges++;
            break;
        case 4:
            nrNonManifold++;
            break;
        default:
            nrWrongEdges++;
            break;
        }
    }
    
    // print
    std::cout << " ... Mesh consistency" << std::endl;
    std::cout << " ... nr. boundary edges: " << nrBndEdges << std::endl;
    std::cout << " ... nr. of manifold inner edges: " << nrInnerEdges << std::endl;
    std::cout << " ... nr. of non-manifold edges: " << nrNonManifold << std::endl;
    std::cout << " ... nr. of WRONG edges: " << nrWrongEdges << std::endl;
    ///// test helfedge mesh
    //CTimer timer;
    //hm.halfedges(v.size(), q, he, hef, hev, timer);
    //// copy data back
    //std::vector<int> lf;
    //std::vector<unsigned char> la;
    //hef.getHalfedgeFaces(lf, la);
    //int nrNonManifoldHalfedgeFaces{ 0 };
    //for (int i = 0; i < la.size(); i++)
    //{
    //    bool flag = la[i] & 0x20;
    //    if (flag) nrNonManifoldHalfedgeFaces++;
    //}
    //std::cout << " ... nr. of non-manifold halfedge faces: " << nrNonManifoldHalfedgeFaces << std::endl;
    //// get quads
    //std::vector<int4> lq;
    //std::vector<unsigned char> lqa;
    //q.getQuadrilaterals(lq, lqa);
    //int nrNonManifoldQuadrilaterals{ 0 };
    //for (int i = 0; i < lqa.size(); i++)
    //{
    //    bool flag = lqa[i] & 0x20;
    //    if (flag) nrNonManifoldQuadrilaterals++;
    //}
    //std::cout << " ... nr. of non-manifold quadrilaterals: " << nrNonManifoldQuadrilaterals << std::endl;
    //// compute vertex valences
    //std::vector<int4> le;
    //he.getHalfedges(le);
    //for (auto e : le)
    //{

    //}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main function to process a volume data set
void 
p_mc::DualMarchingCubes::dualMC(const float i0, Mesh& mesh, std::map<std::string,int>& config)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// config
    bool valenceFlag{ static_cast<bool>(config["valence"]) };
    bool p3X3YColor{ static_cast<bool>(config["p3X3YColor"]) };
    bool p3X3YOld{ static_cast<bool>(config["p3X3YOld"]) };
    bool p3333{ static_cast<bool>(config["p3333"]) };
    bool elementQuality{ static_cast<bool>(config["elementQuality"]) };
    bool objFlag{ static_cast<bool>(config["objFile"]) };

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Problem size
    const int t_size = ugrid.t_size(); // dims[0] * dims[1] * dims[2];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// measure processing time
	CTimer ctimer;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// processing time
    std::cout << " ... compute isosurface\n";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute quadrilaterals
    // 0. alloc lookup tables
	// 1. alloc memory for hash table
	// 2. alloc memory for vertices
	// 3. alloc memory for quads
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 0. Cell intersection object
    CellIntersection c_(r_pattern, t_ambig);
    cudaCheckError();

    // 1. alloc and init hash table
	const int ht_size = static_cast<int>(20e6);
    QuadrilateralHashTable ht_(ht_size);
	//cudaCheckError();
	uint b_size = MC_BLOCKSIZE;
	uint g_size{ (static_cast<uint>(ht_.size()) + b_size - 1) / b_size };
	init_quadrilateral_hashtable << < g_size, b_size >> > (ht_);
    cudaDeviceSynchronize();
	cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 2. alloc and init vertices
    int nr_v{ 0 };
    Vertices v_(15e6);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 3. alloc and init quadrilaterals
    int nr_q{ 0 };
    Quadrilaterals q_(15e6);
    cudaDeviceSynchronize();
	cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 4. compute iso-surface 
    ctimer.start();
	b_size = MC_BLOCKSIZE;
	g_size = (t_size + b_size - 1) / b_size;
    dual_mc << < g_size, b_size >> > (i0, t_size, ugrid, c_, ht_, v_);
	cudaDeviceSynchronize();
    cudaCheckError();


	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 5. compute shared vertex list for quadrilateral mesh
	// indices of quadrilateral vertices have to be mapped to global vertex index in vertex array
	// get number of vertices
	nr_v = v_.size();
    if (nr_v == 0)
    {
        std::cout << " ERROR: no vertices\n";
        return;
    }
	
	// map quadrilateral indices
	b_size = MC_BLOCKSIZE;
	g_size = (ht_.size() + b_size - 1) / b_size;
	map_quadrilaterals <<< g_size, b_size >>> (ht_, q_);
    cudaDeviceSynchronize();
    cudaCheckError();
	// get number of quadrilaterals
	nr_q = q_.size();
    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute processing time
	ctimer.stop();
    ctimer.print(std::string("dual marching cubes and indexed face set"));
    timesMC.push_back(ctimer.getTime(std::string("Dual MC ")));
        

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate halfedge data structure
    HalfedgeHashTable et_;
    Halfedges he_;
    HalfedgeFaces he_f;
    HalfedgeVertices he_v;
    HalfedgeMesh he_m;
    
    // compute halfedge mesh data structure
    int nr_e = he_m.halfedges(nr_v, q_, he_, he_f, he_v, ctimer);
    ctimer.print(std::string("halfedge data structure"));
    timesHE.push_back(ctimer.getTime(std::string("HE MC ")));
    
    std::cout << " ... elements generated by dual mc" << std::endl;
    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of quadrilaterals " << nr_q << std::endl;
    // check vertex valence distribution
    if (valenceFlag)
    {
        checkVertexValence(v_, q_, std::string(" ... dual MC"));
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Check face coloring method
    std::cout << " ... compute element coloring" << std::endl;
    FaceColoring fc;
    fc.colorFaces(q_, he_, he_f, ctimer);
    ctimer.stop();
    ctimer.print(std::string(" compute face coloring"));
    timesColoring.push_back(ctimer.getTime(std::string("Face Colring ")));
    g_size = q_.size();
    transfer_face_attributes << < g_size, b_size >> > (q_, he_f);
    cudaDeviceSynchronize();
    cudaCheckError();
    
    /// check mesh
    //checkMeshConsistency(v_, q_, std::string( " ... coloring" ));
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Remove elements with valence pattern 3-x-3-y and 3-3-3-3
    MeshSimplification ms;
    if (p3X3YColor)
    {
        std::cout << " ... start mesh simplification 3X3Y" << std::endl;
        ms.pattern3X3Y(v_, q_, he_, he_f, he_v, ctimer);
        ctimer.print(std::string(" simplification 3X3Y color based"));
        times3X3Y.push_back(ctimer.getTime(std::string("P3X3Y color based ")));
        /// check mesh
        //checkMeshConsistency(v_, q_, std::string(" ... color based P3X3Y"));
        nr_v = v_.size();
        nr_q = q_.size();
        std::cout << " ... nr. of vertices: " << nr_v << std::endl;
        std::cout << " ... nr. of faces: " << nr_q << std::endl;
        ///  check if vertex valence distribution has changed
        if (valenceFlag)
        {
            checkVertexValence(v_, q_, std::string(" ... P3X3Y color based"));
        }
    }
    if (p3X3YOld)
    {
        ms.pattern3X3YOld(v_, q_, he_, he_f, he_v, ctimer);
        ctimer.print(std::string(" simplification 3X3Y old"));
        times3X3Y.push_back(ctimer.getTime(std::string("P3X3Y old method ")));
        /// check mesh
        //checkMeshConsistency(v_, q_, std::string(" ... old P3X3Y"));
        nr_v = v_.size();
        nr_q = q_.size();
        std::cout << " ... nr. of vertices: " << nr_v << std::endl;
        std::cout << " ... nr. of faces: " << nr_q << std::endl;
        ///  check if vertex valence distribution has changed
        if (valenceFlag)
        {
            checkVertexValence(v_, q_, std::string(" ... P3X3Y old method"));
        }
    }
    if (p3333)
    {
        std::cout << " ... start mesh simplification 3333" << std::endl;
        ms.pattern3333(v_, q_, he_, he_f, he_v, ctimer);
        ctimer.print(std::string(" simplification 3333"));
        times3333.push_back(ctimer.getTime(std::string("P3333 ")));
        //checkMeshConsistency(v_, q_, std::string(" ... P3333"));
        nr_v = v_.size();
        nr_q = q_.size();
        std::cout << " ... nr. of vertices: " << nr_v << std::endl;
        std::cout << " ... nr. of faces: " << nr_q << std::endl;
        ///  check if vertex valence distribution has changed
        if (valenceFlag)
        {
            checkVertexValence(v_, q_, std::string(" ... P3333"));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate a triangle mesh by optimal subdivision of quadrilaterals into trinalges
    nr_v = v_.size();
    nr_q = q_.size();
    const int nr_t = 2 * nr_q;
    Triangles t_(nr_t);
    g_size = (static_cast<uint>(t_.a_size) + b_size - 1) / b_size;
    quadrilateral_to_triangle<<< g_size, b_size>>>(q_, v_, t_);
    cudaDeviceSynchronize();
    cudaCheckError();
    t_.nr_t = nr_t;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate edges from quadrilaterals, this is for rendering purposes
    Edges qe_;
    int nr_qe = he_m.edges(q_, qe_);

    // Resume
    std::cout << " ... Total number of elements generated " << std::endl;
    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of quadrilaterals " << nr_q << std::endl;
    std::cout << " ... total nr. of triangles " << nr_t << std::endl;
    std::cout << " ... total nr. of edges " << nr_qe << std::endl;

	// Copy data from device
	// copy first elements
    mesh.resizeQuadrilaterals(nr_q);
    mesh.resizeTriangles(nr_t);
    mesh.resizeLines(nr_qe);
    cudaMemcpy(mesh.quadrilateralsData(), q_.quadrilaterals, nr_q * sizeof(int4), cudaMemcpyDeviceToHost);
    cudaMemcpy(mesh.trianglesData(), t_.triangles, nr_t * sizeof(int3), cudaMemcpyDeviceToHost); 
    cudaMemcpy(mesh.linesData(), qe_.edges, nr_qe * sizeof(int2), cudaMemcpyDeviceToHost); 

  
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// copy vertices to mesh
    // Note: vertices has to be copied, they are double in mesh and float in CUDA
    mesh.resizeVertices(nr_v);
    float3* v_array = new float3[nr_v];
    float3* n_array = new float3[nr_v];
    cudaMemcpy(v_array, v_.vertices, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    cudaMemcpy(n_array, v_.normals, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
	for (int id = 0; id < nr_v; id++) {
		// copy vertices
		Mesh::Vertex v{ v_array[id].x,v_array[id].y, v_array[id].z };
		Mesh::Normal n{ -n_array[id].x, -n_array[id].y, -n_array[id].z };
		mesh.addVertex(id, v, n);
        const double rC = 202.0 / 255.0;
        const double gC = 200.0 / 255.0;
        const double bC = 201.0 / 255.0;
        const double aC = 1.0;
        Mesh::Attribute a{ rC,gC,bC,aC };
        mesh.attribute(id, a);
	}
    delete[] v_array;
    delete[] n_array;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // copy half edge data structure from device to host
    nr_v = he_v.size(); // number of halfedge vertices 
    nr_q = he_f.size(); // number of halfedge faces
    nr_e = he_.size(); // number of halfedge edges
    std::vector<int4> he_e_array;
    std::vector<int> he_v_array;
    std::vector<int> he_f_array;
    std::vector<uchar> he_f_attributes;
    he_.getHalfedges(he_e_array);
    he_v.getHalfedgeVertices(he_v_array);
    he_f.getHalfedgeFaces(he_f_array, he_f_attributes);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // collect infos
    //info = new float3[nr_v];
    //cudaMemcpy(info, v_.error, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    //info_.reset(info);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // copy face colors from quadrilaterals
    //checkFaceColoring(he_e_array, he_f_array, he_f_attributes);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Check halfedge data structure
    //checkHalfedge(he_e_array, he_f_array, mesh);
    
    // compute elment quality
    if (elementQuality)
    {
        float* qm_data{ nullptr };
        QualityMeasure m_(nr_t);
        //allocQualityMeasure(m_, nr_t);
        g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
        mean_ratio_measure << < g_size, b_size >> > (t_, v_, m_);
        // save quality measure data
        qm_data = new float[nr_t];
        cudaMemcpy(qm_data, m_.q_measure, nr_t * sizeof(float), cudaMemcpyDeviceToHost);
        // write to file
        std::ofstream o_file("./data/models/qualityDMC.bin", std::ofstream::binary);
        o_file.write((char*)&nr_t, sizeof(int));
        o_file.write((char*)qm_data, nr_t * sizeof(float));
        o_file.close();
        delete[] qm_data;
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// done!
	std::cout << " ... done\n";
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// write off file
    //objFlag = false;
    if (objFlag) {
        std::ofstream objF;
        objF.open("./data/models/dualMC_result.obj");
        if (!objF.is_open()) {
            std::cout << "ERROR: can't open output file " << std::endl;
        }
        else {
            objF << "#Dual Marching Cubes\n";
            for (int i = 0; i < nr_v; i++)
            {
                Mesh::Vertex v = mesh.vertex(i);
                Mesh::Attribute a = mesh.attribute(i);
                objF << "v " << v[0] << " " << v[1] << " " << v[2] << " " << a[0] << " " << a[1] << " " << a[2] << std::endl;
            }
            for (auto n : mesh.getNormals())
            {
                objF << "vn " << n[0] << " " << n[1] << " " << n[2] << std::endl;
            }
            //for (auto t : mesh.getTriangles())
            for (auto q : mesh.getQuadrilaterals())
            {
                //objF << "f " << (t[0] + 1) << "//" << (t[0] + 1) << " " << (t[1] + 1) << "//" << (t[1] + 1) << " " << (t[2] + 1) << "//" << (t[2] + 1) << std::endl;
                objF << "f " << (q[0] + 1) << "//" << (q[0] + 1) << " " << (q[1] + 1) << "//" << (q[1] + 1) << " " << (q[2] + 1) << "//" << (q[2] + 1) << " " << (q[3] + 1) << "//" << (q[3] + 1) << std::endl;
            }
            objF.close();
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// standard marching cubes
void
p_mc::DualMarchingCubes::standardMC(const float i0, Mesh& mesh)
{
    std::cout << " ... compute MC surface\n";
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Problem size
    const int t_size = ugrid.t_size(); // dims[0] * dims[1] * dims[2];
    // measure processing time
    CTimer ctimer1;
    CTimer ctimer2;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // alloc and init data structures
    VertexHashTable ht_(25000000);
    Vertices v_(10000000);
    Triangles t_(20000000);
    MarchingCubesLookupTables l_tables(e_pattern, t_pattern);

    uint b_size = MC_BLOCKSIZE;
    uint g_size = (static_cast<uint>(ht_.t_size) + b_size - 1) / b_size;
    init_vertex_hashtable << < g_size, b_size >> > (ht_);
    cudaDeviceSynchronize();
    cudaCheckError();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // processing time
    std::cout << " ... compute iso-surface\n";
    ctimer1.start();

    // compute MC
    g_size = (t_size + b_size - 1) / b_size;
    standard_mc << < g_size, b_size >> > (i0, t_size, ugrid, l_tables, ht_, v_, t_);
    cudaDeviceSynchronize();
    cudaCheckError();

    // collect tirangles
    int nr_v = v_.size(); // size<Vertices>(v_);
    int nr_t = t_.size(); // size<Triangles>(t_);
    t_.nr_t = nr_t;
    g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
    map_triangles << < g_size, b_size >> > (ht_, t_);
    cudaDeviceSynchronize();
    cudaCheckError();

    ctimer1.stop();
    ctimer1.print(std::string("Standard Marching Cubes"));
    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of triangles " << nr_t << std::endl;

    // compute elment quality
    float* qm_data{ nullptr };
    
    QualityMeasure m_(nr_t);
    //allocQualityMeasure(m_, nr_t);
    g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
    mean_ratio_measure << < g_size, b_size >> > (t_, v_, m_);
    // save quality measure data
    qm_data = new float[nr_t];
    cudaMemcpy(qm_data, m_.q_measure, nr_t * sizeof(float), cudaMemcpyDeviceToHost);
    // write to file
    std::ofstream o_file("./data/models/qualityMC.bin", std::ofstream::binary);
    o_file.write((char*)& nr_t, sizeof(int));
    o_file.write((char*)qm_data, nr_t * sizeof(float));
    o_file.close();
    delete[] qm_data;
    
    // copy data into mesh
    // Copy data from device
    float3* v_array = new float3[nr_v];
    float3* n_array = new float3[nr_v];
    int3* t_array = new int3[nr_t];

    cudaMemcpy(v_array, v_.vertices, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    cudaMemcpy(n_array, v_.normals, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    cudaMemcpy(t_array, t_.triangles, nr_t * sizeof(int3), cudaMemcpyDeviceToHost);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate mesh structure
    mesh.resizeVertices(nr_v);
    for (int id = 0; id < nr_v; id++) {
        // copy vertices
        Mesh::Vertex v{ v_array[id].x,v_array[id].y, v_array[id].z };
        Mesh::Normal n{ -n_array[id].x, -n_array[id].y, -n_array[id].z };
        mesh.addVertex(id, v, n);
    }
    // generate two triangles for each quad, directX can render only triangles
    mesh.resizeTriangles(nr_t);
    for (int id = 0; id < nr_t; id++) {
        const int v0 = t_array[id].x;
        const int v1 = t_array[id].y;
        const int v2 = t_array[id].z;
        Mesh::Triangle t1{ v0, v1, v2 };
        mesh.addTriangle(id, t1);
    }
}

void p_mc::DualMarchingCubes::checkFaceColoring(std::vector<int4>& he_e_array, std::vector<int>& he_f_array, std::vector<uchar>& fc_array)
{
    int cnt{ 0 };
    int a0{ 0 };
    int a1{ 0 };
    int a2{ 0 };
    int a3{ 0 };
    int a4{ 0 };
    int cA{ 0 };
    int cSimpl{ 0 };
    auto getColor = [](uchar c) { return static_cast<int>(c&0x1F); };
    // loop over all quadrilaterals
    const int nr_q = static_cast<int>(he_f_array.size());
    for (int f = 0; f < nr_q; f++)
    {
        // collect halfedge
        const int e0 = he_f_array[f];
        const int e1 = he_e_array[e0].z;
        const int e2 = he_e_array[e1].z;
        const int e3 = he_e_array[e2].z;
        // collect all four twin edges
        const int t0 = he_e_array[e0].w;
        const int t1 = he_e_array[e1].w;
        const int t2 = he_e_array[e2].w;
        const int t3 = he_e_array[e3].w;
        // collect faces, if any
        int f0{ -1 };
        int f1{ -1 };
        int f2{ -1 };
        int f3{ -1 };
        if (t0 > -1) f0 = he_e_array[t0].y;
        if (t1 > -1) f1 = he_e_array[t1].y;
        if (t2 > -1) f2 = he_e_array[t2].y;
        if (t3 > -1) f3 = he_e_array[t3].y;
        // collect colors
        int c = getColor(fc_array[f]);
        int c0{ -1 };
        int c1{ -1 };
        int c2{ -1 };
        int c3{ -1 };
        if (f0 > -1) c0 = getColor(fc_array[f0]);
        if (f1 > -1) c1 = getColor(fc_array[f1]);
        if (f2 > -1) c2 = getColor(fc_array[f2]);
        if (f3 > -1) c3 = getColor(fc_array[f3]);
        // compare with own color
        bool flag{ false };
        if (c == c0) flag = true;
        if (c == c1) flag = true;
        if (c == c2) flag = true;
        if (c == c3) flag = true;
        if (c == 0) a0++;
        if (c == 1) a1++;
        if (c == 2) a2++;
        if (c == 3) a3++;
        if (c == 4) a4++;
        if (c > 4) cA++;
        if (c == 4) 
        {
            int k{ 0 };
            if (c0 == 0 || c1 == 0 || c2 == 0 || c3 == 0) k++;
            if (c0 == 1 || c1 == 1 || c2 == 1 || c3 == 1) k++;
            if (c0 == 2 || c1 == 2 || c2 == 2 || c3 == 2) k++;
            if (c0 == 3 || c1 == 3 || c2 == 3 || c3 == 3) k++;
            if (k == 4) cSimpl++;
        }
        if (flag)
        {
            cnt++;
        }
    }
    std::cout << " ... ERROR: " << cnt << " faces do not comply coloring condition" << std::endl;
    std::cout << " ... Color 0: " << a0 << " faces have color 0" << std::endl;
    std::cout << " ... Color 1: " << a1 << " faces have color 1" << std::endl;
    std::cout << " ... Color 2: " << a2 << " faces have color 2" << std::endl;
    std::cout << " ... Color 3: " << a3 << " faces have color 3" << std::endl;
    std::cout << " ... Color 4: " << a4 << " faces have color 4" << std::endl;
    std::cout << " ... Color A: " << cA << " faces have color A" << std::endl;
    std::cout << " ... Color complete: " << cSimpl << " faces have all four colors as neighbors" << std::endl;
}

void p_mc::DualMarchingCubes::checkHalfedge(std::vector<int4>& he_e_array, std::vector<int>& he_f_array, Mesh& mesh)
{
    /** halfedge int4:
    //    he.x = origin vertex
    //    he.y = face
    //    he.z = next
    //    he.w = twin
    */
    std::map<int, int> m_;
    int c{ 0 };
    const int nr_e = static_cast<int>(he_e_array.size());
    for (int e = 0; e < nr_e; e++)
    {
        int4 he0 = he_e_array[e];
        if (he0.w == -1)
        {
            // found boundary
            const int v = he_e_array[e].x;
            Mesh::Attribute a{ 1,0,0,1 };
            mesh.attribute(v, a);
            c++;
            m_.insert(std::make_pair(v, e));
        }
        
        // construct quadrilateral, compare with quads from mesh.
        const int ihe0 = e;
        const int ihe1 = he0.z;
        int4 he1 = he_e_array[ihe1];
        const int ihe2 = he1.z;
        int4 he2 = he_e_array[ihe2];
        const int ihe3 = he2.z;
        int4 he3 = he_e_array[ihe3];
        int4 heTest = he_e_array[ihe3];
        const int f = he0.y;
        Mesh::Quadrilateral q = mesh.quadrilateral(f);
        const int v0 = he0.x;
        const int v1 = he1.x;
        const int v2 = he2.x;
        const int v3 = he3.x;
        // compare vertices
        if (q[0] != v0 && q[0] != v1 && q[0] != v2 && q[0] != v3)
        {
            std::cout << "ERROR: can't find v0\n";
        }
        if (q[1] != v0 && q[1] != v1 && q[1] != v2 && q[1] != v3)
        {
            std::cout << "ERROR: can't find v1\n";
        }
        if (q[2] != v0 && q[2] != v1 && q[2] != v2 && q[2] != v3)
        {
            std::cout << "ERROR: can't find v2\n";
        }
        if (q[3] != v0 && q[3] != v1 && q[3] != v2 && q[3] != v3)
        {
            std::cout << "ERROR: can't find v3\n";
        }

        // check twins
        const int t0 = he0.w;
        const int t1 = he1.w;
        const int t2 = he2.w;
        const int t3 = he3.w;
        int tv0{ -1 };
        int tv1{ -1 };
        int tv2{ -1 };
        int tv3{ -1 };
        if (t0 != -1)
        {
            tv0 = he_e_array[t0].w;
            int ni = he_e_array[t0].z;
            int4 ne = he_e_array[ni];
            if (tv0 != ihe0)
            {
                std::cout << "ERROR: wrong twin for he 0\n";
            }
            if (v0 != ne.x)
            {
                std::cout << "ERROR: wrong vertex config by twin edges\n";
            }

        }
        if (t1 != -1)
        {
            tv1 = he_e_array[t1].w;
            if (tv1 != ihe1)
            {
                std::cout << "ERROR: wrong twin for he 1\n";
            }

        }
        if (t2 != -1)
        {
            tv2 = he_e_array[t2].w;
            if (tv2 != ihe2)
            {
                std::cout << "ERROR: wrong twin for he 0\n";
            }

        }
        if (t3 != -1)
        {
            tv3 = he_e_array[t3].w;
            if (tv3 != ihe3)
            {
                std::cout << "ERROR: wrong twin for he 0\n";
            }

        }


    }
    std::cout << " ... nr. of bnd edges: " << c << std::endl;
    return;
    // find edges
    std::vector<std::array<int, 2>> edges;
    for (auto e : m_)
    {
        const int n = he_e_array[e.second].z;
        const int vn = he_e_array[n].x;
        // check if vertex is in list
        auto search = m_.find(vn);
        if (search != m_.end()) {
            //std::cout << "Found " << search->first << " " << search->second << '\n';
            edges.push_back({ e.first,vn });
        }
        else {
            std::cout << "Not found\n";
        }
    }
    std::cout << "found: " << edges.size() << ", edges" << std::endl;
    // clear double cases
    auto  setKey = [](const int v0, const int v1)
    {
        if (v0 < v1)
            return (static_cast<unsigned long long>(v0) << 32) | (v1 & 0xffffffffL);
        else
            return (static_cast<unsigned long long>(v1) << 32) | (v0 & 0xffffffffL);
    };
    std::map<unsigned long long, std::array<int, 2>> em_;
    for (auto e : edges)
    {
        em_.insert({setKey(e[0] ,e[1]),e});
    }
    std::cout << "after clear: " << em_.size() << ", edges" << std::endl;
    // analyze cases
    /*std::ofstream f;
    f.open("./data/models/nonManifold.txt");
    auto boundary = [](const int i, const int j, const int k, UGrid& u)
    {
        bool flag{ false };
        if (i == 0 || i == u.i_size()) flag = true;
        if (j == 0 || j == u.j_size()) flag = true;
        if (k == 0 || k == u.k_size()) flag = true;
        return flag;
    };
    auto scalar = [](const int i_, const int j_, const int k_, float* v, UGrid& u)
    {
        std::array<double, 8> s;
        for (int k = 0; k < 2; k++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int i = 0; i < 2; i++)
                {
                    int index = (i & 1) | (j & 1) << 1 | (k & 1) << 2;
                    int gl = u.gl_index(i_ + i, j_ + j, k_ + k);
                    s[index] = v[gl];
                }
            }
        }
        return s;
    };*/
    //for (auto e : em_)
    //{
    //    const int v0 = e.second[0];
    //    const int v1 = e.second[1];
    //    const int c0 = info[v0].x;
    //    const int c1 = info[v1].x;
    //    const int g0 = info[v0].z;
    //    const int g1 = info[v1].z;
    //    const int i0 = ugrid.i_index(g0);
    //    const int j0 = ugrid.j_index(g0);
    //    const int k0 = ugrid.k_index(g0);
    //    const int i1 = ugrid.i_index(g1);
    //    const int j1 = ugrid.j_index(g1);
    //    const int k1 = ugrid.k_index(g1);
    //    if (!boundary(i0, j0, k0, ugrid) && !boundary(i1, j1, k1, ugrid))
    //    {
    //        auto s0 = scalar(i0, j0, k0, h_volume, ugrid);
    //        auto s1 = scalar(i1, j1, k1, h_volume, ugrid);
    //        f << s0[0] << "," << s0[1] << "," << s0[2] << "," << s0[3] << "," << s0[4] << "," << s0[5] << "," << s0[6] << "," << s0[7] << "," << i0 << "," << j0 << "," << k0 << std::endl;
    //        f << s1[0] << "," << s1[1] << "," << s1[2] << "," << s1[3] << "," << s1[4] << "," << s1[5] << "," << s1[6] << "," << s1[7] << "," << i1 << "," << j1 << "," << k1 << std::endl;
    //    }


    //    //std::cout << "(" << c0 << ", " << i0 << ", " << j0 << ", " << k0 << "), (" << c1 << ", " << i1 << ", " << j1 << ", " << k1 << "), " << std::endl;
    //}
    //f.close();
}

void p_mc::DualMarchingCubes::writeTimes( ) 
{
    std::vector<std::string> f(5);
    f[0] = "./data/measurements/cuda__MC.txt";
    f[1] = "./data/measurements/cuda__HE.txt";
    f[2] = "./data/measurements/cuda_Coloring.txt";
    f[3] = "./data/measurements/cuda__3X3Y.txt";
    f[4] = "./data/measurements/cuda__3333.txt";
    std::map<std::string, std::vector<std::string>> m;
    m.insert({ f[0], timesMC });
    m.insert({ f[1], timesHE });
    m.insert({ f[2], timesColoring });
    m.insert({ f[3], times3X3Y });
    m.insert({ f[4], times3333 });

    for (auto n : m)
    {
        std::ofstream o_file(n.first.c_str(), std::ios::out | std::ios::app);
        for (auto s : n.second) {
            o_file << s;
        }
        o_file.close();
    }

}



