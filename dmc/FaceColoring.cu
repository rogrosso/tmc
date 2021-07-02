#include "FaceColoring.h"



__global__ void color_distribution(p_mc::Quadrilaterals q, p_mc::FaceColoring::ColorCount fc)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q.nr_q) return;
    const int c = q.getColor(tid);
    if (c < p_mc::C) return;
    fc.addColor(c, tid);
}
__global__ void color_faces(const int nr_q, const int c, p_mc::FaceColoring::ColorCount fc, p_mc::Quadrilaterals q_, p_mc::Halfedges he_, p_mc::HalfedgeFaces he_f)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_q) return;
    const int q_id = fc.getFace(c, tid);
    // collect halfedge from face
    const int e0 = he_f.he_e[q_id];
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
        q_.setColor(q_id, 0);
        return;
    }
    if (c0 != 1 && c1 != 1 && c2 != 1 && c3 != 1)
    {
        q_.setColor(q_id, 1);
        return;
    }
    if (c0 != 2 && c1 != 2 && c2 != 2 && c3 != 2)
    {
        q_.setColor(q_id, 2);
        return;
    }
    if (c0 != 3 && c1 != 3 && c2 != 3 && c3 != 3)
    {
        q_.setColor(q_id, 3);
        return;
    }
    if (c0 != 4 && c1 != 4 && c2 != 4 && c2 != 4)
    {
        q_.setColor(q_id, 4);
        return;
    }
}
__global__ void optimize_face_coloring(p_mc::Quadrilaterals q_, p_mc::Halfedges he_, p_mc::HalfedgeFaces he_f)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q) return;
    if (q_.getColor(tid) != 4) return; // we are processing color 4
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

__host__ void p_mc::FaceColoring::colorFaces(Quadrilaterals& q, Halfedges he, HalfedgeFaces& hef, CTimer& timer)
{
    const int nr_q = q.size();
    timer.start();
    ColorCount fc(nr_q);
    int b_size = MC_BLOCKSIZE;
    int g_size = (nr_q + b_size - 1) / b_size;
    color_distribution << < g_size, b_size >> > (q, fc);
    for (int c = 5; c < 24; c++)
    {
        const int nr_f = fc.size(c);
        if (nr_f > 0)
        {
            g_size = (nr_f + b_size - 1) / b_size;
            color_faces << < g_size, b_size >> > (nr_f, c, fc, q, he, hef);
            cudaDeviceSynchronize();
            cudaCheckError();
        }
    }
    g_size = (nr_q + b_size - 1) / b_size;
    optimize_face_coloring << < g_size, b_size >> > (q, he, hef);
    cudaDeviceSynchronize();
    cudaCheckError();
    timer.stop();
}
