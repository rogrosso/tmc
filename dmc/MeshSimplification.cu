#include "MeshSimplification.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      GLOBAL (Kernels)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init vertex map to smooth mesh
__global__ void init_VertexMap(p_mc::VertexMap v_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= v_.size()) return;
    v_.init(tid);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init quadrilateral map to smooth mesh
__global__ void init_QuadrilateralMap(p_mc::QuadrilateralMap m_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= m_.size()) return;
    m_.init(tid);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute vertex valence using edges
//__global__ void vertex_valence(p_mc::Edges e_, p_mc::VertexMap v_)
//{
//    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
//    if (tid >= e_.nr_e)
//        return;
//    const int v0 = e_.edges[tid].x;
//    const int v1 = e_.edges[tid].y;
//    atomicAdd(&v_.valence[v0], 1);
//    atomicAdd(&v_.valence[v1], 1);
//}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute vertex valence using halfedges
//    he.x = origin vertex
//    he.y = face
//    he.z = next
//    he.w = twin
__global__ void vertex_valence(p_mc::Halfedges e_, p_mc::VertexMap vm_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e_.size())
        return;
    const int v0 = e_.he_e[tid].x;
    atomicAdd(&vm_.valence[v0], 1);
    if (e_.he_e[tid].w == p_mc::INVALID_INDEX)
    {
        // count valence for vertex pointed by halfedge
        const int nv = e_.he_e[e_.he_e[tid].z].x;
        atomicAdd(&vm_.valence[nv], 1);
        //printf("couting valence at boundary\n");
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Mesh simplification
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mark elements with valence 3X3Y, X, Y >= 5
__global__ void mark_elements_P3X3Y(p_mc::Quadrilaterals q_, p_mc::VertexMap vm_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (q_.isNonManifold(tid)) return;
    q_.clearPatternBits(tid);

    const int valence0 = vm_.valence[q_.v0(tid)];
    const int valence1 = vm_.valence[q_.v1(tid)];
    const int valence2 = vm_.valence[q_.v2(tid)];
    const int valence3 = vm_.valence[q_.v3(tid)];

    bool flag1 = (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5);
    bool flag2 = (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3);
    if (flag1 || flag2)
    {
        q_.setPatternBit(tid);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute vertex count, i.e. neighbor elements are of type 3-X-3-Y sharing the same valence 3 vertex
__global__ void count_vertices_P3X3Y(p_mc::Quadrilaterals q_, p_mc::VertexMap v_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    q_.clearPatternBits(tid);
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;

    const int valence0 = v_.valence[v0];
    const int valence1 = v_.valence[v1];
    const int valence2 = v_.valence[v2];
    const int valence3 = v_.valence[v3];


    if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5) {
        atomicAdd(&v_.count[v0], 1);
        atomicAdd(&v_.count[v2], 1);
        q_.setPatternBit(tid);
    }
    else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3) {
        atomicAdd(&v_.count[v1], 1);
        atomicAdd(&v_.count[v3], 1);
        q_.setPatternBit(tid);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// merge vertices if possible
// - if two elements with the pattern 3-X-3-Y are neighbors, they cannot be removed, this is indicated by count, count > 1 means more than
//   one element are sharering a vertex of valence 3
// - if the element can be removed, move the vertices with valence 3 to the midpoint of the line connecting the other two vertices
//   e.g. one can allways move the element v0 (valence(v0) = 3) or v1 (valence v1 = 3)
__global__ void merge_vertices_P3X3Y(p_mc::Quadrilaterals q_, p_mc::VertexMap m_, p_mc::Vertices v_, p_mc::Halfedges e_, p_mc::HalfedgeFaces f_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (q_.isNonManifold(tid)) return;
    if (!q_.isP3X3Y(tid)) return;
    int4 n = e_.findNeighbors(f_.he_e[tid]);
    if (q_.isNonManifold(n.x)) return;
    if (q_.isNonManifold(n.y)) return;
    if (q_.isNonManifold(n.z)) return;
    if (q_.isNonManifold(n.w)) return;

    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;

    const int valence0 = m_.valence[v0];
    const int valence1 = m_.valence[v1];
    const int valence2 = m_.valence[v2];
    const int valence3 = m_.valence[v3];

    if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5 && m_.count[v0] == 1 && m_.count[v2] == 1)
    {
        // compute new position and normals of v0
        v_.vertices[v0].x = v_.vertices[v1].x + 0.5 * (v_.vertices[v3].x - v_.vertices[v1].x);
        v_.vertices[v0].y = v_.vertices[v1].y + 0.5 * (v_.vertices[v3].y - v_.vertices[v1].y);
        v_.vertices[v0].z = v_.vertices[v1].z + 0.5 * (v_.vertices[v3].z - v_.vertices[v1].z);
        v_.normals[v0].x = v_.normals[v1].x + 0.5 * (v_.normals[v3].x - v_.normals[v1].x);
        v_.normals[v0].y = v_.normals[v1].y + 0.5 * (v_.normals[v3].y - v_.normals[v1].y);
        v_.normals[v0].z = v_.normals[v1].z + 0.5 * (v_.normals[v3].z - v_.normals[v1].z);

        float sz_ = v_.normals[v0].x * v_.normals[v0].x + v_.normals[v0].y * v_.normals[v0].y + v_.normals[v0].z * v_.normals[v0].z;
        sz_ = sqrtf(sz_);

        v_.normals[v0].x = v_.normals[v0].x / sz_;
        v_.normals[v0].y = v_.normals[v0].y / sz_;
        v_.normals[v0].z = v_.normals[v0].z / sz_;

        // mark v2 to be removed
        m_.type[v2] = p_mc::P_REMOVE;

        // set twins of v2 to be v0, to be able to remove element later
        m_.twin[v2] = v0;

        // element has to be removed
        q_.setRemoveBit(tid);
    }
    else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3 && m_.count[v1] == 1 && m_.count[v3] == 1)
    {
        // compute new position and normal of v1
        v_.vertices[v1].x = v_.vertices[v0].x + 0.5 * (v_.vertices[v2].x - v_.vertices[v0].x);
        v_.vertices[v1].y = v_.vertices[v0].y + 0.5 * (v_.vertices[v2].y - v_.vertices[v0].y);
        v_.vertices[v1].z = v_.vertices[v0].z + 0.5 * (v_.vertices[v2].z - v_.vertices[v0].z);
        v_.normals[v1].x = v_.normals[v0].x + 0.5 * (v_.normals[v2].x - v_.normals[v0].x);
        v_.normals[v1].y = v_.normals[v0].y + 0.5 * (v_.normals[v2].y - v_.normals[v0].y);
        v_.normals[v1].z = v_.normals[v0].z + 0.5 * (v_.normals[v2].z - v_.normals[v0].z);

        float sz_ = v_.normals[v1].x * v_.normals[v1].x + v_.normals[v1].y * v_.normals[v1].y + v_.normals[v1].z * v_.normals[v1].z;
        sz_ = sqrtf(sz_);

        v_.normals[v1].x = v_.normals[v1].x / sz_;
        v_.normals[v1].y = v_.normals[v1].y / sz_;
        v_.normals[v1].z = v_.normals[v1].z / sz_;

        // mark v3 to be removed
        m_.type[v3] = p_mc::P_REMOVE;
        // set twins, remove v3, use addres of v1
        m_.twin[v3] = v1;

        // element has to be removed
        q_.setRemoveBit(tid);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// merge vertices baesd on element color
// - Color based simplification: if one neighbor have same vertex valence patter but a higher color, the element can be removed
// - If one neighbor have the same vertex valence pattern and a smaller color, the element cannot be removed
__global__ void merge_vertices_P3X3Y_color(p_mc::Quadrilaterals q_, p_mc::VertexMap vm_, p_mc::Vertices v_, p_mc::Halfedges e_, p_mc::HalfedgeFaces f_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (!q_.isP3X3Y(tid)) return;

    const int c = q_.getColor(tid);
    int4 n = e_.findNeighbors(f_.he_e[tid]);
    if (q_.isNonManifold(n.x)) return;
    if (q_.isNonManifold(n.y)) return;
    if (q_.isNonManifold(n.z)) return;
    if (q_.isNonManifold(n.w)) return;
    // get neighbor colors
    int n_color[4];
    n_color[0] = q_.getColor(n.x);
    n_color[1] = q_.getColor(n.y);
    n_color[2] = q_.getColor(n.z);
    n_color[3] = q_.getColor(n.w);
    // get neighbor vertex valence pattern
    bool n_pattern[4];
    n_pattern[0] = q_.isP3X3Y(n.x);
    n_pattern[1] = q_.isP3X3Y(n.y);
    n_pattern[2] = q_.isP3X3Y(n.z);
    n_pattern[3] = q_.isP3X3Y(n.w);
    // case analysis
    bool flag{ true };
    for (int i = 0; i < 4; i++)
    {
        if (n_pattern[i] && (n_color[i] <= c)) flag = false;
    }
    if (!flag) return;
    //if (c != 0) return;

    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;

    const int valence0 = vm_.valence[v0];
    const int valence1 = vm_.valence[v1];
    const int valence2 = vm_.valence[v2];
    const int valence3 = vm_.valence[v3];

    if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5)
    {
        // compute new position and normals of v0
        v_.vertices[v0].x = v_.vertices[v1].x + 0.5 * (v_.vertices[v3].x - v_.vertices[v1].x);
        v_.vertices[v0].y = v_.vertices[v1].y + 0.5 * (v_.vertices[v3].y - v_.vertices[v1].y);
        v_.vertices[v0].z = v_.vertices[v1].z + 0.5 * (v_.vertices[v3].z - v_.vertices[v1].z);
        v_.normals[v0].x = v_.normals[v1].x + 0.5 * (v_.normals[v3].x - v_.normals[v1].x);
        v_.normals[v0].y = v_.normals[v1].y + 0.5 * (v_.normals[v3].y - v_.normals[v1].y);
        v_.normals[v0].z = v_.normals[v1].z + 0.5 * (v_.normals[v3].z - v_.normals[v1].z);

        float sz_ = v_.normals[v0].x * v_.normals[v0].x + v_.normals[v0].y * v_.normals[v0].y + v_.normals[v0].z * v_.normals[v0].z;
        sz_ = sqrtf(sz_);

        v_.normals[v0].x = v_.normals[v0].x / sz_;
        v_.normals[v0].y = v_.normals[v0].y / sz_;
        v_.normals[v0].z = v_.normals[v0].z / sz_;

        // mark v2 to be removed
        vm_.type[v2] = p_mc::P_REMOVE;

        // set twins of v2 to be v0, to be able to remove element later
        vm_.twin[v2] = v0;

        // element has to be removed
        q_.setRemoveBit(tid);
    }
    else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3)
    {
        // compute new position and normal of v1
        v_.vertices[v1].x = v_.vertices[v0].x + 0.5 * (v_.vertices[v2].x - v_.vertices[v0].x);
        v_.vertices[v1].y = v_.vertices[v0].y + 0.5 * (v_.vertices[v2].y - v_.vertices[v0].y);
        v_.vertices[v1].z = v_.vertices[v0].z + 0.5 * (v_.vertices[v2].z - v_.vertices[v0].z);
        v_.normals[v1].x = v_.normals[v0].x + 0.5 * (v_.normals[v2].x - v_.normals[v0].x);
        v_.normals[v1].y = v_.normals[v0].y + 0.5 * (v_.normals[v2].y - v_.normals[v0].y);
        v_.normals[v1].z = v_.normals[v0].z + 0.5 * (v_.normals[v2].z - v_.normals[v0].z);

        float sz_ = v_.normals[v1].x * v_.normals[v1].x + v_.normals[v1].y * v_.normals[v1].y + v_.normals[v1].z * v_.normals[v1].z;
        sz_ = sqrtf(sz_);

        v_.normals[v1].x = v_.normals[v1].x / sz_;
        v_.normals[v1].y = v_.normals[v1].y / sz_;
        v_.normals[v1].z = v_.normals[v1].z / sz_;

        // mark v3 to be removed
        vm_.type[v3] = p_mc::P_REMOVE;
        // set twins, remove v3, use addres of v1
        vm_.twin[v3] = v1;

        // element has to be removed
        q_.setRemoveBit(tid);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// new vertices
__global__ void remove_vertices_P3X3Y(const int nr_v, p_mc::Vertices v_, p_mc::VertexMap m_, p_mc::Vertices n_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_v)
        return;
    //
    if (m_.type[tid] != p_mc::P_REMOVE) {
        // copy vertex to new list
        const int addr = atomicAdd(n_.t_size, 1);
        n_.vertices[addr] = v_.vertices[tid];
        n_.normals[addr] = v_.normals[tid];
        // keep address for mapping
        m_.map_addr[tid] = addr;
    }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// new elements
__global__ void remove_quadrilaterals_P3X3Y(p_mc::Quadrilaterals q_, p_mc::VertexMap m_, p_mc::Quadrilaterals n_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    // check if quadrilateral has to be removed
    if (q_.isRemove(tid)) return; // element is removed from the list

    // compute new quadrilateral
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    int4 nq_;
    if (m_.type[v0] == p_mc::P_REMOVE) nq_.x = m_.map_addr[m_.twin[v0]];
    else  nq_.x = m_.map_addr[v0];
    if (m_.type[v1] == p_mc::P_REMOVE) nq_.y = m_.map_addr[m_.twin[v1]];
    else  nq_.y = m_.map_addr[v1];
    if (m_.type[v2] == p_mc::P_REMOVE) nq_.z = m_.map_addr[m_.twin[v2]];
    else  nq_.z = m_.map_addr[v2];
    if (m_.type[v3] == p_mc::P_REMOVE) nq_.w = m_.map_addr[m_.twin[v3]];
    else  nq_.w = m_.map_addr[v3];

    // search an address to store element
    const int addr = n_.addQuadrilateral(nq_);
    if (n_.isNonManifold(tid)) n_.setNonManifold(addr);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// use halfedge data structure pattern is 3-3-3-3, mark element as type P_3333, and vertices as type P_3333
// if neighbor element is of the same type, the element can't be removed
// mark neighbor for remove: P_NEIG
__global__ void mark_elements_P3333(p_mc::HalfedgeFaces q_, p_mc::Halfedges he_, p_mc::VertexMap vm_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >=  q_.size()) return;
    if (q_.isNonManifold(tid)) return;
    q_.clearPatternBits(tid);

    // find element vertices
    const int e0 = q_.he_e[tid];
    const int e1 = he_.he_e[e0].z;
    const int e2 = he_.he_e[e1].z;
    const int e3 = he_.he_e[e2].z;
    const int v0 = he_.he_e[e0].x;
    const int v1 = he_.he_e[e1].x;
    const int v2 = he_.he_e[e2].x;
    const int v3 = he_.he_e[e3].x;
    const int valence0 = vm_.valence[v0];
    const int valence1 = vm_.valence[v1];
    const int valence2 = vm_.valence[v2];
    const int valence3 = vm_.valence[v3];
    if (valence0 != 3 || valence1 != 3 || valence2 != 3 || valence3 != 3) return;

    // find element neighbors
    // f0
    int twin = he_.he_e[e0].w;
    if (twin == -1) return; // boundary element
    const int f0 = he_.he_e[twin].y;
    int next = he_.he_e[twin].z;
    next = he_.he_e[next].z;
    const int v4 = he_.he_e[next].x;
    // f1
    twin = he_.he_e[e1].w;
    if (twin == -1) return; // boundary element
    const int f1 = he_.he_e[twin].y;
    next = he_.he_e[twin].z;
    next = he_.he_e[next].z;
    const int v5 = he_.he_e[next].x;
    // f2
    twin = he_.he_e[e2].w;
    if (twin == -1) return; // boundary element
    const int f2 = he_.he_e[twin].y;
    next = he_.he_e[twin].z;
    next = he_.he_e[next].z;
    const int v6 = he_.he_e[next].x;
    // f3
    twin = he_.he_e[e3].w;
    if (twin == -1) return; // boundary element
    const int f3 = he_.he_e[twin].y;
    next = he_.he_e[twin].z;
    next = he_.he_e[next].z;
    const int v7 = he_.he_e[next].x;

    // check neighbors for non-manifold
    if (q_.isNonManifold(f0)) return;
    if (q_.isNonManifold(f1)) return;
    if (q_.isNonManifold(f2)) return;
    if (q_.isNonManifold(f3)) return;

    // compute valence
    const int valence4 = vm_.valence[v4];
    const int valence5 = vm_.valence[v5];
    const int valence6 = vm_.valence[v6];
    const int valence7 = vm_.valence[v7];

    // check if neighor is of the sampe type
    bool flag0 = (valence4 == 3 && valence5 == 3);
    bool flag1 = (valence5 == 3 && valence6 == 3);
    bool flag2 = (valence6 == 3 && valence7 == 3);
    bool flag3 = (valence7 == 3 && valence4 == 3);
    if (flag0 || flag1 || flag2 || flag3) {
        // element can't be removed
        return;
    }
    // valence will be reduced by one by neighbor vertices,
    // if one vertex has valence 3, simplification can't be done,
    // it will result in one vertex with valence 2
    if (valence4 == 3 || valence5 == 3 || valence6 == 3 || valence7 == 3) {
        // element can't be removed
        return;
    }
    // check for special case, where removing element would
    // generate a non-manifold mesh
    if (valence4 == 4 && valence5 == 4 && valence6 == 4 && valence7 == 4) {
        return;
    }

    // mark element as type P_3333
    q_.setPatternBit(tid);
    // mark neighbors, they will be removed if possible
    q_.setRemoveBit(f0);
    q_.setRemoveBit(f1);
    q_.setRemoveBit(f2);
    q_.setRemoveBit(f3);
    // vertex type
    vm_.type[v0] = p_mc::P_3333;
    vm_.type[v1] = p_mc::P_3333;
    vm_.type[v2] = p_mc::P_3333;
    vm_.type[v3] = p_mc::P_3333;
    // vertex mappring
    vm_.twin[v0] = v4;
    vm_.twin[v1] = v5;
    vm_.twin[v2] = v6;
    vm_.twin[v3] = v7;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// remove vertices
__global__ void remove_vertices_P3333(p_mc::Vertices v_, p_mc::VertexMap vm_, p_mc::Vertices n_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= v_.nr_v) return;
    if (vm_.type[tid] == p_mc::P_3333) return;

    const int addr = atomicAdd(n_.t_size, 1);
    n_.vertices[addr] = v_.vertices[tid];
    n_.normals[addr] = v_.normals[tid];
    // keep address for mapping
    vm_.map_addr[tid] = addr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// remove elements
__global__ void remove_quadrilaterals_P3333(p_mc::Quadrilaterals q_, p_mc::HalfedgeFaces f_, p_mc::VertexMap vm_, p_mc::Quadrilaterals n_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (f_.isRemove(tid)) return;

    // compute new quadrilateral
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    int4 nq_;
    if (f_.isP3333(tid))
    {
        nq_.x = vm_.map_addr[vm_.twin[v0]];
        nq_.y = vm_.map_addr[vm_.twin[v1]];
        nq_.z = vm_.map_addr[vm_.twin[v2]];
        nq_.w = vm_.map_addr[vm_.twin[v3]];
    }
    else
    {
        nq_.x = vm_.map_addr[v0];
        nq_.y = vm_.map_addr[v1];
        nq_.z = vm_.map_addr[v2];
        nq_.w = vm_.map_addr[v3];
    }

    // search an address to store element
    const int addr = n_.addQuadrilateral(nq_);
    if (q_.isNonManifold(tid)) n_.setNonManifold(addr);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      HOST CODE
//
//      Mesh simplification
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simplify elements with valence pattern 3X3Y
// Assume the mesh data structures are complete
void p_mc::MeshSimplification::pattern3X3Y(Vertices& v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& hef, HalfedgeVertices& hev, CTimer& timer)
{
    const int nr_v = v.size();
    const int nr_q = q.size();
    const int nr_e = he.size();
    Vertices nv(nr_v);
    Quadrilaterals nq(nr_q);
    VertexMap vm(nr_v);

    // measure time
    timer.start();
    // initialize data structures
    int b_size = MC_BLOCKSIZE;
    int g_size = (nr_v + b_size - 1) / b_size;
    init_VertexMap << < g_size, b_size >> > (vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // compute vertex valence
    g_size = (nr_e + b_size - 1) / b_size;
    vertex_valence << < g_size, b_size >> > (he, vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // mark elemenets and vertices for removal
    g_size = (nr_q + b_size - 1) / b_size;
    mark_elements_P3X3Y << < g_size, b_size >> > (q, vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // merge vertices
    g_size = (nr_q + b_size - 1) / b_size;
    merge_vertices_P3X3Y_color << < g_size, b_size >> > (q, vm, v, he, hef);
    cudaDeviceSynchronize();
    cudaCheckError();
    // remove vertices from list
    g_size = (nr_v + b_size - 1) / b_size;
    remove_vertices_P3X3Y << < g_size, b_size >> > (nr_v, v, vm, nv);
    cudaDeviceSynchronize();
    cudaCheckError();
    // remove quadrilaterals
    g_size = (nr_q + b_size - 1) / b_size;
    remove_quadrilaterals_P3X3Y << < g_size, b_size >> > (q, vm, nq);
    cudaDeviceSynchronize();
    cudaCheckError();
    // measure time
    timer.stop();
    // copy data back
    q.copy(nq);
    v.copy(nv);
    // re-compute halfedge data structure
    //HalfedgeMesh hm;
    //CTimer t;
    //hm.halfedges(v.size(), q, he, hef, hev, t);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simplify elements with valence pattern 3X3Y
// Assume the mesh data structures are complete
void p_mc::MeshSimplification::pattern3X3YOld(Vertices& v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& hef, HalfedgeVertices& hev, CTimer& timer)
{
    const int nr_v = v.size();
    const int nr_q = q.size();
    const int nr_e = he.size();
    Vertices nv(nr_v);
    Quadrilaterals nq(nr_q);
    VertexMap vm(nr_v);

    // measure time
    timer.start();
    // initialize data structures
    int b_size = MC_BLOCKSIZE;
    int g_size = (nr_v + b_size - 1) / b_size;
    init_VertexMap << < g_size, b_size >> > (vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // compute vertex valence
    g_size = (nr_e + b_size - 1) / b_size;
    vertex_valence << < g_size, b_size >> > (he, vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // mark elemenets and vertices for removal
    g_size = (nr_q + b_size - 1) / b_size;
    count_vertices_P3X3Y << < g_size, b_size >> > (q, vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // merge vertices
    g_size = (nr_q + b_size - 1) / b_size;
    merge_vertices_P3X3Y << < g_size, b_size >> > (q, vm, v, he, hef);
    cudaDeviceSynchronize();
    cudaCheckError();
    // remove vertices from list
    g_size = (nr_v + b_size - 1) / b_size;
    remove_vertices_P3X3Y << < g_size, b_size >> > (nr_v, v, vm, nv);
    cudaDeviceSynchronize();
    cudaCheckError();
    // remove quadrilaterals
    g_size = (nr_q + b_size - 1) / b_size;
    remove_quadrilaterals_P3X3Y << < g_size, b_size >> > (q, vm, nq);
    cudaDeviceSynchronize();
    cudaCheckError();
    // measure time
    timer.stop();
    // copy data back
    q.copy(nq);
    v.copy(nv);
    // re-compute halfedge data structure
    //HalfedgeMesh hm;
    //CTimer t;
    //hm.halfedges(v.size(), q, he, hef, hev, t);
}
void p_mc::MeshSimplification::pattern3333(Vertices& v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& hef, HalfedgeVertices& hev, CTimer& timer)
{
    const int nr_v = v.size();
    const int nr_q = q.size();
    const int nr_e = he.size();
    Vertices nv(nr_v);
    Quadrilaterals nq(nr_q);
    VertexMap vm(nr_v);

    // measure time
    timer.start();
    // initialize data structures
    int b_size = MC_BLOCKSIZE;
    int g_size = (nr_v + b_size - 1) / b_size;
    init_VertexMap << < g_size, b_size >> > (vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // compute vertex valence
    g_size = (nr_e + b_size - 1) / b_size;
    vertex_valence << < g_size, b_size >> > (he, vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // mark elements with vertex valence pattern 3333
    g_size = (nr_q + b_size) / b_size;
    mark_elements_P3333 << < g_size, b_size >> > (hef, he, vm);
    cudaDeviceSynchronize();
    cudaCheckError();
    // remove vertices
    g_size = (nr_v + b_size) / b_size;
    remove_vertices_P3333 << < g_size, b_size >> > (v, vm, nv);
    cudaDeviceSynchronize();
    cudaCheckError();
    // remove elements
    g_size = (nr_q + b_size - 1) / b_size;
    remove_quadrilaterals_P3333 << < g_size, b_size >> > (q, hef, vm, nq);
    cudaDeviceSynchronize();
    cudaCheckError();
    // measure time
    timer.stop();

    // copy data back
    q.copy(nq);
    v.copy(nv);
}
