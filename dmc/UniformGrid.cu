#include "UniformGrid.h"


struct uGrid {
    int nx;
    int ny;
    int nz;
    float dx;
    float dy;
    float dz;
    float xmin;
    float ymin;
    float zmin;
    __device__ int size() { return nx * ny * nz; }
    __device__ int i_index(const int gl_index) { return (gl_index % nx); }
    __device__ int j_index(const int gl_index) { return ((gl_index / nx) % ny); }
    __device__ int k_index(const int gl_index) { return (gl_index / (nx * ny)); }
    __device__ float3 cellVertex(const int i, const int j, const int k) { return make_float3(xmin + i * dx, ymin + j * dy, zmin + k * dz); }
};

__global__ void volume(uGrid ugrid, float* d_scalar, p_mc::UniformGrid::SurfaceCase sc)
{
    using SurfaceCase = p_mc::UniformGrid::SurfaceCase;
    float pi{ 3.14159265358979323846f };
    // use a 1d grid
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (ugrid.size() <= gl_index)
        return;

    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    if (i_index >= ugrid.nx || j_index >= ugrid.ny || k_index >= ugrid.nz)
    {
        return;
    }
    float val = 0.0f;
    float3 v = ugrid.cellVertex(i_index, j_index, k_index);
    auto sq = [](const float x) { return x * x;  };
    auto qu = [](const float x) { return x * x * x * x; };
    auto torus_h = [sq](float3 pos, float3 center, float2 param)
    {
        const float c = sq(param.x);
        const float a = sq(param.y);
        const float x = sq(pos.x - center.x);
        const float y = sq(pos.y - center.y);
        const float z = sq(pos.z - center.z);
        return sq(x + y + z + c - a) - 4 * c * (x + y);
    };
    auto torus_v = [sq](float3 pos, float3 center, float2 param)
    {
        const float c = sq(param.x);
        const float a = sq(param.y);
        const float x = sq(pos.x - center.x);
        const float y = sq(pos.y - center.y);
        const float z = sq(pos.z - center.z);
        return sq(x + y + z + c - a) - 4 * c * (x + z);
    };
    auto genusTwo = [sq](const float3 pos)
    {
        float alpha = 1.0;
        float x = (pos.x + 1.0f) / 2.0f;
        float y = (pos.y + 1.0f) / 2.0f;
        float z = (pos.z + 1.0f) / 2.0f;
        x = alpha * (4.0f * x - 2.0f);
        y = alpha * (4.0f * y - 2.0f);
        z = alpha * (4.0f * z - 2.0f);
        float t = 2 * y * (y * y - 3 * x * x) * (1 - z * z);
        t += (x * x + y * y) * (x * x + y * y);
        t -= (9 * z * z - 1) * (1 - z * z);
        return t;
    };
    auto iWP = [pi](const float3 p)
    {
        const float alpha = 5.01;
        //const float alpha = 1.01;
        const float x = alpha * (p.x + 1) * pi;
        const float y = alpha * (p.y + 1) * pi;
        const float z = alpha * (p.z + 1) * pi;
        return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z); // iso-value = 0

    };
    auto pwHybrid = [pi](const float3 p)
    {
        //const float alpha = 3.01;
        const float alpha = 1.01;
        const float x = alpha * (p.x + 1) * pi;
        const float y = alpha * (p.y + 1) * pi;
        const float z = alpha * (p.z + 1) * pi;
        return 4.0f * (cosf(x) * cosf(y) + cosf(y) * cosf(z) + cosf(z) * cosf(x)) -  3* cosf(x) * cosf(y) * cosf(z) + 0.8f; // iso-value = 0
    };
    auto neovius = [pi](const float3 p)
    {
        const float alpha = 1;
        const float x = alpha * (p.x + 1) * pi;
        const float y = alpha * (p.y + 1) * pi;
        const float z = alpha * (p.z + 1) * pi;
        return 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z); // iso_value = 0.0
    };
    auto goursat = [sq,qu](const float3 p)
    {
        const float a = -1.0f;
        const float b = 0.0f;
        const float c = 0.5f;
        return qu(p.x) + qu(p.y) + qu(p.z) + a * (sq(p.x) + sq(p.y) + sq(p.z)) + b * (sq(p.x) + sq(p.y) + sq(p.z)) + c;

    };
    auto steinerRoman = [sq, qu](const float3 p)
    {
        const float alpha = 1.5f;
        const float x = alpha * p.x;
        const float y = alpha * p.y;
        const float z = alpha * p.z;
        return sq(x*x + y*y + z*z - 1.0f) - (sq(z - 1) - 2.0f * x*x) * (sq(z + 1) - 2 * y*y);
    };
    switch (sc) {
    case SurfaceCase::Sphere:
        val = v.x * v.x + v.y * v.y + v.z * v.z - 0.16;
        break;
    case SurfaceCase::Torus:
    {
        const float2 param{ 0.3,0.15 };
        const float3 center{ 0,0,0 };
        val = torus_h(v, center, param);
        break;
    }
    case SurfaceCase::TwoHoledTorus:
    {
        const float2 p1{ 0.3,0.15 };
        const float t1 = 0.38;
        const float t2 = 0.2;
        const float delta = 0.38;
        const float vt1 = torus_h(v, make_float3(-t1, 0, 0), p1);
        const float vt2 = torus_h(v, make_float3(t2, delta, 0), p1);
        val = fminf(vt1, vt2);
        break;
    }
    case SurfaceCase::FourHoledTorus:
    {
        const float2 p2{ 0.3,0.15 };
        const float t = 0.38;
        const float v1 = torus_h(v, make_float3(-t, 0, 0), p2);
        const float v2 = torus_h(v, make_float3(t, 0, 0), p2);
        const float v3 = torus_v(v, make_float3(0, 0, -t), p2);
        const float v4 = torus_v(v, make_float3(0, 0, t), p2);
        val = fminf(v1, v2);
        val = fminf(val, v3);
        val = fminf(val, v4);
        break;
    }
    case SurfaceCase::GenusTwo:
        val = genusTwo(v);
        break;
    case SurfaceCase::Goursat:
        val = goursat(v);
        break;
    case SurfaceCase::iWP:
        val = iWP(v);
        break;
    case SurfaceCase::pwHybrid:
        val = pwHybrid(v);
        break;
    case SurfaceCase::neovius:
        val = neovius(v);
        break;
    case SurfaceCase::SteinerRoman:
        val = steinerRoman(v);
        break;
    default:
        val = genusTwo(v);
    }

    d_scalar[gl_index] = val;
}


__host__ void p_mc::UniformGrid::generateVolume(const std::array<int, 3>& dim, SurfaceCase sc)
{
    // volume size
    idim = dim[0];
    jdim = dim[1];
    kdim = dim[2];

    // domain
    float xmin = -1.0f;
    float ymin = -1.0f;
    float zmin = -1.0f;
    float xmax = 1.0f;
    float ymax = 1.0f;
    float zmax = 1.0f;

    // define grid size
    dx = (xmax - xmin) / (idim - 1.0f);
    dy = (ymax - ymin) / (jdim - 1.0f);
    dz = (zmax - zmin) / (kdim - 1.0f);
    x0 = xmin;
    y0 = ymin;
    z0 = zmin;

    uGrid vol;
    vol.nx = idim;
    vol.ny = jdim;
    vol.nz = kdim;
    vol.dx = dx;
    vol.dy = dy;
    vol.dz = dz;
    vol.xmin = x0;
    vol.ymin = y0;
    vol.zmin = z0;

    // allocate data
    cudaMalloc(&d_scalar, t_size() * sizeof(float));
    cudaCheckError();
    d_scalar_.reset(d_scalar, cudaFree);
    // compute volume
    const size_t size_ = t_size();
    uint b_size = MC_BLOCKSIZE;
    uint g_size = (static_cast<uint>(size_) + b_size - 1) / b_size;
    //volume << < g_size, b_size >> > (ugrid, sCase);
    volume << < g_size, b_size >> > (vol, d_scalar, sc);
    cudaDeviceSynchronize();
    cudaCheckError();
}
