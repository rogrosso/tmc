#include "VertexValence.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// initialize hash table for edges
__global__ void init_valence_hashtable(p_mc::ValenceHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size) return;
    ht_.init(tid);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Save edges into hash table and keep track which face shares the edge, up to four faces are expected for non-manifold cases
__global__ void collect_quads_edges(p_mc::Quadrilaterals q_, p_mc::ValenceHashTable ht_)
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
    ht_.addEdge(v0, v1);
    ht_.addEdge(v1, v2);
    ht_.addEdge(v2, v3);
    ht_.addEdge(v3, v0);
}
__global__ void collect_tris_edges(p_mc::Triangles t_, p_mc::ValenceHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= t_.nr_t)
        return;
    // collect vertices
    const int v0 = t_.triangles[tid].x;
    const int v1 = t_.triangles[tid].y;
    const int v2 = t_.triangles[tid].z;

    // add edges to hash table
    ht_.addEdge(v0, v1);
    ht_.addEdge(v1, v2);
    ht_.addEdge(v2, v0);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// count vertex valence
__global__ void vertex_valence(p_mc::ValenceHashTable h_, unsigned int* v_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= h_.size()) return;
    if (h_.empty(tid)) return;
    // collect vertices
    const unsigned int v0 = h_.v0(tid);
    const unsigned int v1 = h_.v1(tid);

    // add edges to hash table
    atomicAdd(&v_[v0], 1);
    atomicAdd(&v_[v1], 1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute valence distribution
__global__ void vertex_valence_distribution(const int nr_v, unsigned int* v_, unsigned int* valence_distribution)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_v) return;
    // collect vertices
    const unsigned int index = v_[tid];

    // add edges to hash table
    atomicAdd(&valence_distribution[index], 1);
}



void p_mc::VertexValence::vertexValence(const int nr_v, Quadrilaterals& q_, std::vector<uint>& valence)
{
    const int nr_q{ q_.size() };
    const int nr_h{ 4 * nr_q };
    ValenceHashTable h(nr_h);
    const int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(nr_h) + b_size - 1) / b_size;
    init_valence_hashtable << < g_size, b_size >> > (h);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
    collect_quads_edges << < g_size, b_size >> > (q_, h);
    cudaDeviceSynchronize();
    cudaCheckError();
    // allocate an integer array to compute valence of vertices
    unsigned int* v{ nullptr };
    cudaMalloc(&v, nr_v * sizeof(unsigned int));
    cudaMemset(v, 0, nr_v * sizeof(unsigned int));
    cudaCheckError();
    // compute vertex valence
    g_size = (static_cast<uint>(nr_h) + b_size - 1) / b_size;
    vertex_valence << < g_size, b_size >> > (h, v);
    // compute valence distribution
    unsigned int* valence_distribution{ nullptr };
    const int max_expected_valence{ 100 }; // max valence 99, starting at 0
    cudaMalloc(&valence_distribution, max_expected_valence * sizeof(unsigned int));
    cudaMemset(valence_distribution, 0, max_expected_valence * sizeof(unsigned int));
    g_size = (static_cast<uint>(nr_v) + b_size - 1) / b_size;
    vertex_valence_distribution << < g_size, b_size >> > (nr_v, v, valence_distribution);
    cudaCheckError();
    // copy data to device
    valence.resize(max_expected_valence);
    cudaMemcpy(valence.data(), valence_distribution, max_expected_valence * sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaCheckError();
    // check largest index with non-vanishing valence
    int index{ 0 };
    for (int i = 1; i < max_expected_valence; i++)
    {
        if (valence[i] > 0) index = i;
    }
    valence.resize(index + 1);
    // free memory
    cudaFree(v);
    cudaFree(valence_distribution);
}


void p_mc::VertexValence::vertexValence(const int nr_v, Triangles& t_, std::vector<uint>& valence)
{
    const int nr_t{ t_.size() };
    const int nr_h{ 4 * nr_t };
    ValenceHashTable h(nr_h);
    const int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(nr_h) + b_size - 1) / b_size;
    init_valence_hashtable << < g_size, b_size >> > (h);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
    collect_tris_edges << < g_size, b_size >> > (t_, h);
    cudaDeviceSynchronize();
    cudaCheckError();
    // allocate an integer array to compute valence of vertices
    unsigned int* v{ nullptr };
    cudaMalloc(&v, nr_v * sizeof(unsigned int));
    cudaMemset(v, 0, nr_v * sizeof(unsigned int));
    cudaCheckError();
    // compute vertex valence
    g_size = (static_cast<uint>(nr_h) + b_size - 1) / b_size;
    vertex_valence << < g_size, b_size >> > (h, v);
    // compute valence distribution
    unsigned int* valence_distribution{ nullptr };
    const int max_expected_valence{ 100 }; // max valence 99, starting at 0
    cudaMalloc(&valence_distribution, max_expected_valence * sizeof(unsigned int));
    cudaMemset(valence_distribution, 0, max_expected_valence * sizeof(unsigned int));
    g_size = (static_cast<uint>(nr_v) + b_size - 1) / b_size;
    vertex_valence_distribution << < g_size, b_size >> > (nr_v, v, valence_distribution);
    cudaCheckError();
    // copy data to device
    valence.resize(max_expected_valence);
    cudaMemcpy(valence.data(), valence_distribution, max_expected_valence * sizeof(unsigned int), cudaMemcpyDeviceToHost);
    cudaCheckError();
    // check largest index with non-vanishing valence
    int index{ 0 };
    for (int i = 1; i < max_expected_valence; i++)
    {
        if (valence[i] > 0) index = i;
    }
    valence.resize(index + 1);
    // free memory
    cudaFree(v);
    cudaFree(valence_distribution);
}
