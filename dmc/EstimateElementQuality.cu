#include "EstimateElementQuality.h"

__global__ void quadQuality(const int nr_q, p_mc::Vertices v_, p_mc::Quadrilaterals q_, p_mc::ElementQuality eQ)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_q) return;
    // collect vertices
    int4 q;
    q.x = q_.v0(tid);
    q.y = q_.v1(tid);
    q.z = q_.v2(tid);
    q.w = q_.v3(tid);
    float3 v[4];
    v[0] = v_[q.x];
    v[1] = v_[q.y];
    v[2] = v_[q.z];
    v[3] = v_[q.w];
    eQ.quads(tid, v);
}

__global__ void triQuality(const int nr_t, p_mc::Vertices v_, p_mc::Triangles t_, p_mc::ElementQuality eQ)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_t) return;
    // collect vertices
    int3 t;
    t.x = t_.v0(tid);
    t.y = t_.v1(tid);
    t.z = t_.v2(tid);
    float3 v[3];
    v[0] = v_[t.x];
    v[1] = v_[t.y];
    v[2] = v_[t.z];
    eQ.tris(tid, v);
}


void p_mc::EstimateElementQuality::q(Vertices v_, Quadrilaterals q_, ElementQuality eQ)
{
    const int nr_q = q_.size();
    const int nr_e = eQ.size();
    if (nr_q != nr_e) return; // wrong size
    const int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(nr_q) + b_size) / (b_size - 1);
    quadQuality << < g_size, b_size >> > (nr_q, v_, q_, eQ);
}

void p_mc::EstimateElementQuality::q(Vertices v_, Triangles t_, ElementQuality eQ)
{
    const int nr_t = t_.size();
    const int nr_e = eQ.size();
    if (nr_t != nr_e) return; // wrong size
    const int b_size = MC_BLOCKSIZE;
    int g_size = (static_cast<uint>(nr_t) + b_size) / (b_size - 1);
    triQuality << < g_size, b_size >> > (nr_t, v_, t_, eQ);
}
