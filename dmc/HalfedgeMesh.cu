#include "HalfedgeMesh.h"


// Defines
#define MC_BLOCKSIZE 128

/// global functions

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for edges
__global__ void init_edge_hashtable(p_mc::EdgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size) return;
    ht_.init(tid);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for edges
__global__ void init_halfedge_hashtable(p_mc::HalfedgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size) return;
    ht_.init(tid);

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Save edges into hash table and keep track which face shares the edge, up to four faces are expected for non-manifold cases
__global__ void collect_edges(p_mc::Quadrilaterals q_, p_mc::EdgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q)
        return;
    // collect vertices
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;

    // add edges to hash table
    ht_.addEdge(v0, v1, tid);
    ht_.addEdge(v1, v2, tid);
    ht_.addEdge(v2, v3, tid);
    ht_.addEdge(v3, v0, tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id for edges, required for rendering purposes
__global__ void map_edges(p_mc::EdgeHashTable ht_, p_mc::Edges e_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.size())
        return;
    if (ht_.empty(tid))
        return;
    // add only if quadrilateral is complete
    e_.addEdge(ht_.edge[tid].x, ht_.edge[tid].y);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// count non-manifold edges
__global__ void non_manifold(p_mc::EdgeHashTable e, int* t_size)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e.size()) return;
    if (e.empty(tid)) return;
    if (e.nrFaces(tid) > 2) // up to four faces can share an edge
    {
        atomicAdd(t_size, 1);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// collect halfedges, faces and vertices
// Loop over quadrilaterals
__global__ void collect_halfedge_elements(p_mc::Quadrilaterals q_, p_mc::Halfedges he_, p_mc::HalfedgeFaces he_f, p_mc::HalfedgeVertices he_v, p_mc::EdgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q)
        return;
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    // get an address and save the half edge
    // there are four halfedges for each quadrilateral
    const int he_addr = 4 * tid;
    //    he.x = origin vertex
    //    he.y = face
    //    he.z = next
    //   he.w = tween
    // 1.
    he_.he_e[he_addr].x = v0;
    he_.he_e[he_addr].y = tid;
    he_.he_e[he_addr].z = he_addr + 1;
    he_.he_e[he_addr].w = -1;
    // 2.
    he_.he_e[he_addr + 1].x = v1;
    he_.he_e[he_addr + 1].y = tid;
    he_.he_e[he_addr + 1].z = he_addr + 2;
    he_.he_e[he_addr + 1].w = -1; // twins will be computed later
    // 3.
    he_.he_e[he_addr + 2].x = v2;
    he_.he_e[he_addr + 2].y = tid;
    he_.he_e[he_addr + 2].z = he_addr + 3;
    he_.he_e[he_addr + 2].w = -1; // twins will be computed later
    // 4.
    he_.he_e[he_addr + 3].x = v3;
    he_.he_e[he_addr + 3].y = tid;
    he_.he_e[he_addr + 3].z = he_addr;
    he_.he_e[he_addr + 3].w = -1; // twins will be computed later

    // collect faces
    he_f.he_e[tid] = he_addr;
    he_f.attributes[tid] = q_.attributes[tid];

    // collect vertices, don't care about race conditions
    // last writing get the index
    he_v.he_e[v0] = he_addr; // v0 is the origin of he_addr
    he_v.he_e[v1] = he_addr + 1; // v1 is the origin of he_addr+1
    he_v.he_e[v2] = he_addr + 2; // v2 is the origin of he_addr+2
    he_v.he_e[v3] = he_addr + 3; // v3 is the origin of he_addr+3

    // set hash tables to compute connectivity
    ht_.addEdge(v0, v1, tid);
    ht_.addEdge(v1, v2, tid);
    ht_.addEdge(v2, v3, tid);
    ht_.addEdge(v3, v0, tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// collect halfedges, faces and vertices
// Loop over quadrilaterals
__global__ void collect_halfedge_elements_he(p_mc::Quadrilaterals q_, p_mc::Halfedges he_, p_mc::HalfedgeFaces he_f, p_mc::HalfedgeVertices he_v, p_mc::HalfedgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q)
        return;
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    // get an address and save the half edge
    // there are four halfedges for each quadrilateral
    const int he_addr = 4 * tid;
    //    he.x = origin vertex
    //    he.y = face
    //    he.z = next
    //   he.w = tween
    // 1.
    he_.he_e[he_addr].x = v0;
    he_.he_e[he_addr].y = tid;
    he_.he_e[he_addr].z = he_addr + 1;
    he_.he_e[he_addr].w = -1;
    // 2.
    he_.he_e[he_addr + 1].x = v1;
    he_.he_e[he_addr + 1].y = tid;
    he_.he_e[he_addr + 1].z = he_addr + 2;
    he_.he_e[he_addr + 1].w = -1; // twins will be computed later
    // 3.
    he_.he_e[he_addr + 2].x = v2;
    he_.he_e[he_addr + 2].y = tid;
    he_.he_e[he_addr + 2].z = he_addr + 3;
    he_.he_e[he_addr + 2].w = -1; // twins will be computed later
    // 4.
    he_.he_e[he_addr + 3].x = v3;
    he_.he_e[he_addr + 3].y = tid;
    he_.he_e[he_addr + 3].z = he_addr;
    he_.he_e[he_addr + 3].w = -1; // twins will be computed later

    // collect faces
    he_f.he_e[tid] = he_addr;
    he_f.attributes[tid] = q_.attributes[tid];

    // collect vertices, don't care about race conditions
    // last writing get the index
    he_v.he_e[v0] = he_addr; // v0 is the origin of he_addr
    he_v.he_e[v1] = he_addr + 1; // v1 is the origin of he_addr+1
    he_v.he_e[v2] = he_addr + 2; // v2 is the origin of he_addr+2
    he_v.he_e[v3] = he_addr + 3; // v3 is the origin of he_addr+3

    // set hash tables to compute connectivity
    ht_.addHalfedgeTwins(v0, v1, he_addr);
    ht_.addHalfedgeTwins(v1, v2, he_addr + 1);
    ht_.addHalfedgeTwins(v2, v3, he_addr + 2);
    ht_.addHalfedgeTwins(v3, v0, he_addr + 3);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// connect half edge twins, mark non-manifold elements, i.e. elements at the boundary and faces with non-manifold edges
__global__ void connect_halfedge_twins(p_mc::EdgeHashTable e_, p_mc::HalfedgeFaces f_, p_mc::Halfedges he_, p_mc::Quadrilaterals q_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e_.size())
        return;
    if (e_.empty(tid))
        return;

    const int v0 = e_.v0(tid);
    const int v1 = e_.v1(tid);
    const int nrFaces = e_.nrFaces(tid);
    switch (nrFaces)
    {
    case 1:
    {
        f_.setNonManifold(e_.f0(tid));
        q_.setNonManifold(e_.f0(tid));
        break;
    }
    case 2:
    {
        const int he0 = f_.he_e[e_.f0(tid)];
        const int he1 = f_.he_e[e_.f1(tid)];
        const int e0 = he_.findEdge(he0, v0, v1);
        const int e1 = he_.findEdge(he1, v0, v1);
        // connect
        he_.he_e[e0].w = e1;
        he_.he_e[e1].w = e0;
        break;
    }
    case 4:
    {
        const int he0 = f_.he_e[e_.f0(tid)];
        const int he1 = f_.he_e[e_.f1(tid)];
        const int he2 = f_.he_e[e_.f2(tid)];
        const int he3 = f_.he_e[e_.f3(tid)];
        const int e0 = he_.findEdge(he0, v0, v1);
        const int e1 = he_.findEdge(he1, v0, v1);
        const int e2 = he_.findEdge(he2, v0, v1);
        const int e3 = he_.findEdge(he3, v0, v1);
        const int v0 = he_.he_e[e0].x;
        const int v1 = he_.he_e[e1].x;
        const int v2 = he_.he_e[e2].x;
        if (v0 != v1)
        {
            // connect
            he_.he_e[e0].w = e1;
            he_.he_e[e1].w = e0;
            he_.he_e[e2].w = e3;
            he_.he_e[e3].w = e2;
        }
        else if (v0 != v2)
        {
            he_.he_e[e0].w = e2;
            he_.he_e[e2].w = e0;
            he_.he_e[e1].w = e3;
            he_.he_e[e3].w = e1;
        }
        else
        {
            he_.he_e[e0].w = e3;
            he_.he_e[e3].w = e0;
            he_.he_e[e1].w = e2;
            he_.he_e[e2].w = e1;
        }
        // mark face
        f_.setNonManifold(e_.f0(tid));
        f_.setNonManifold(e_.f1(tid));
        f_.setNonManifold(e_.f2(tid));
        f_.setNonManifold(e_.f3(tid));
        q_.setNonManifold(e_.f0(tid));
        q_.setNonManifold(e_.f1(tid));
        q_.setNonManifold(e_.f2(tid));
        q_.setNonManifold(e_.f3(tid));
        break;
    }
    } // switch
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// connect half edge twins, mark non-manifold elements, i.e. elements at the boundary and faces with non-manifold edges
__global__ void connect_halfedge_twins_he(p_mc::HalfedgeHashTable e_, p_mc::HalfedgeFaces f_, p_mc::Halfedges he_, p_mc::Quadrilaterals q_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e_.size())
        return;
    if (e_.empty(tid))
        return;

    const int nrFaces = e_.nrFaces(tid);
    switch (nrFaces)
    {
    case 1:
    {
        f_.setNonManifold(he_.getFace(e_.t0(tid)));
        q_.setNonManifold(he_.getFace(e_.t0(tid)));
        break;
    }
    case 2:
    {
        // connect
        const int t0 = e_.t0(tid);
        const int t1 = e_.t1(tid);
        he_.he_e[t0].w = t1;
        he_.he_e[t1].w = t0;
        break;
    }
    case 4:
    {
        const int t0 = e_.t0(tid);
        const int t1 = e_.t1(tid);
        const int t2 = e_.t2(tid);
        const int t3 = e_.t3(tid);
        const int v0 = he_.he_e[t0].x;
        const int v1 = he_.he_e[t1].x;
        const int v2 = he_.he_e[t2].x;
        if (v0 != v1)
        {
            // connect
            he_.he_e[t0].w = t1;
            he_.he_e[t1].w = t0;
            he_.he_e[t2].w = t3;
            he_.he_e[t3].w = t2;
        }
        else if (v0 != v2)
        {
            he_.he_e[t0].w = t2;
            he_.he_e[t2].w = t0;
            he_.he_e[t1].w = t3;
            he_.he_e[t3].w = t1;
        }
        else
        {
            he_.he_e[t0].w = t3;
            he_.he_e[t3].w = t0;
            he_.he_e[t1].w = t2;
            he_.he_e[t2].w = t1;
        }
        // mark face
        f_.setNonManifold(he_.getFace(t0));
        f_.setNonManifold(he_.getFace(t1));
        f_.setNonManifold(he_.getFace(t2));
        f_.setNonManifold(he_.getFace(t3));
        q_.setNonManifold(he_.getFace(t0));
        q_.setNonManifold(he_.getFace(t1));
        q_.setNonManifold(he_.getFace(t2));
        q_.setNonManifold(he_.getFace(t3));
        break;
    } // case 4
    } // switch
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mark elements which are incident to a non-manifold edge
__global__ void mark_nonmanifold_elements(p_mc::EdgeHashTable e_, p_mc::HalfedgeFaces f_, p_mc::Quadrilaterals q_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e_.size())
        return;
    if (e_.empty(tid))
        return;
    switch(e_.nrFaces(tid))
    {
    case 1: // consider boundary edges to be non-manifold, avoid modifications
        //f_.setNonManifold(e_.f0(tid));
        q_.setNonManifold(e_.f0(tid));
        break;
    case 4:
        //f_.setNonManifold(e_.f0(tid));
        q_.setNonManifold(e_.f0(tid));
        //f_.setNonManifold(e_.f1(tid));
        q_.setNonManifold(e_.f1(tid));
        //f_.setNonManifold(e_.f2(tid));
        q_.setNonManifold(e_.f2(tid));
        //f_.setNonManifold(e_.f3(tid));
        q_.setNonManifold(e_.f3(tid));
        break;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// mark elements which are incident to a non-manifold edge
__global__ void mark_nonmanifold_elements_he(p_mc::HalfedgeHashTable e_, p_mc::Halfedges he_, p_mc::HalfedgeFaces f_, p_mc::Quadrilaterals q_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e_.size()) return;
    if (e_.empty(tid)) return;
    const int nrFaces = e_.nrFaces(tid);
    switch (nrFaces)
    {
    case 1:
    {
        f_.setNonManifold(he_.getFace(e_.t0(tid)));
        q_.setNonManifold(he_.getFace(e_.t0(tid)));
        break;
    }
    case 4:
    {
        const int t0 = e_.t0(tid);
        const int t1 = e_.t1(tid);
        const int t2 = e_.t2(tid);
        const int t3 = e_.t3(tid);
        // mark face
        f_.setNonManifold(he_.getFace(t0));
        f_.setNonManifold(he_.getFace(t1));
        f_.setNonManifold(he_.getFace(t2));
        f_.setNonManifold(he_.getFace(t3));
        q_.setNonManifold(he_.getFace(t0));
        q_.setNonManifold(he_.getFace(t1));
        q_.setNonManifold(he_.getFace(t2));
        q_.setNonManifold(he_.getFace(t3));
        break;
    } // case 4
    } // switch
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Class Methods
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void p_mc::HalfedgeMesh::edgeHashTable(Quadrilaterals& q, EdgeHashTable& h)
{
    const int nr_q = q.size();
    h.resize(4 * nr_q);
    int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(h.size()) + b_size - 1) / b_size;
    init_edge_hashtable << < g_size, b_size >> > (h);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
    collect_edges << < g_size, b_size >> > (q, h);
    cudaDeviceSynchronize();
    cudaCheckError();
}

int p_mc::HalfedgeMesh::edges(Quadrilaterals& q, Edges& e)
{
    const int nr_q = q.size();
    EdgeHashTable h;
    h.resize(4 * nr_q);
    e.resize(static_cast<int>(2.5 * nr_q));
    int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(h.size()) + b_size - 1) / b_size;
    init_edge_hashtable << < g_size, b_size >> > (h);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
    collect_edges << < g_size, b_size >> > (q, h);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(h.size()) + b_size - 1) / b_size;
    map_edges << < g_size, b_size >> > (h, e);
    cudaDeviceSynchronize();
    cudaCheckError();
    return e.size();
}

//int p_mc::HalfedgeMesh::halfedges(const int nr_v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& f, HalfedgeVertices& v, CTimer& timer)
//{
//    //CTimer t;
//    //t.start();
//    const int nr_q = q.size();
//    const int nr_e = 4 * nr_q;
//    he.resize(nr_e);
//    f.resize(nr_q);
//    v.resize(nr_v);
//    EdgeHashTable ht(4 * nr_q);
//    const int b_size = MC_BLOCKSIZE;
//    int g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
//    //init_halfedge_hashtable << < g_size, b_size >> > (ht);
//    init_edge_hashtable << < g_size, b_size >> > (ht);
//    cudaDeviceSynchronize();
//    cudaCheckError();
//    //t.stop();
//    //t.print(std::string("Halfedges: init"));
//    // measure timer
//    timer.start();
//    //t.start();
//    // collect all halfedge elements
//    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
//    collect_halfedge_elements << < g_size, b_size >> > (q, he, f, v, ht);
//    cudaDeviceSynchronize();
//    cudaCheckError();
//    //t.stop();
//    //t.print(std::string("Halfedges: collect halfedge elements"));
//    //t.start();
//    // connect halfedges and mark faces as non-manifold
//    g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
//    connect_halfedge_twins << < g_size, b_size >> > (ht, f, he, q);
//    cudaDeviceSynchronize();
//    cudaCheckError();
//    // measure time
//    timer.stop();
//    //t.stop();
//    //t.print(std::string("Halfedges: connect twins, mark non-manifold"));
//
//    // return number of computed halfedges
//    return nr_e;
//}

int p_mc::HalfedgeMesh::halfedges(const int nr_v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& f, HalfedgeVertices& v, CTimer& timer)
{
    //CTimer t;
    //t.start();
    const int nr_q = q.size();
    const int nr_e = 4 * nr_q;
    he.resize(nr_e);
    f.resize(nr_q);
    v.resize(nr_v);
    HalfedgeHashTable ht(static_cast<int>(100./70. * 4 * nr_q)); // use the 70% rule, there are 4 halfedge for each quad, multiply by 100/70
    const int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
    init_halfedge_hashtable << < g_size, b_size >> > (ht);
    cudaDeviceSynchronize();
    cudaCheckError();
    //t.stop();
    //t.print(std::string("Halfedges: init"));
    // measure timer
    timer.start();
    //t.start();
    // collect all halfedge elements
    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
    //collect_halfedge_elements << < g_size, b_size >> > (q, he, f, v, ht);
    collect_halfedge_elements_he << < g_size, b_size >> > (q, he, f, v, ht);
    cudaDeviceSynchronize();
    cudaCheckError();
    //t.stop();
    //t.print(std::string("Halfedges: collect halfedge elements"));
    //t.start();
    // connect halfedges and mark faces as non-manifold
    g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
    connect_halfedge_twins_he << < g_size, b_size >> > (ht, f, he, q);
    cudaDeviceSynchronize();
    cudaCheckError();
    // measure time
    timer.stop();
    //t.stop();
    //t.print(std::string("Halfedges: connect twins, mark non-manifold"));
    // Mark non-manifold elements
    //g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
    //mark_nonmanifold_elements_he << < g_size, b_size >> > (ht, he, f, q);
    //cudaDeviceSynchronize();
    //cudaCheckError();

    // return number of computed halfedges
    return nr_e;
}

int p_mc::HalfedgeMesh::nonManifold(Quadrilaterals& q)
{
    const int nr_q = q.size();
    EdgeHashTable ht(4 * nr_q);
    const int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
    init_edge_hashtable << < g_size, b_size >> > (ht);
    cudaDeviceSynchronize();
    cudaCheckError();
    // compute edges hash table
    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
    collect_edges << < g_size, b_size >> > (q, ht);
    cudaDeviceSynchronize();
    cudaCheckError();
    // count non-manifold edges
    int* t_size{ nullptr };
    cudaMalloc(&t_size, sizeof(int));
    cudaMemset(t_size, 0, sizeof(int));
    cudaCheckError();
    g_size = (static_cast<uint>(ht.size()) + b_size - 1) / b_size;
    non_manifold << < g_size, b_size >> > (ht, t_size);
    cudaDeviceSynchronize();
    cudaCheckError();
    int nr_n{ 0 };
    cudaMemcpy(&nr_n, t_size, sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaCheckError();
    // free memory
    cudaFree(t_size);
    cudaCheckError();
    return nr_n;
}
