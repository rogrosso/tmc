#pragma once

// libs
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <list>
#include <vector>
#include <array>
#include <stack>
#include <map>
#include <utility>
#include <bitset>
#include <algorithm>
#include <random>
#include <chrono>
#include <cmath>
#include <memory>

// cuda stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "helper_cuda.h"

namespace p_mc {
    /// <summary>
    /// A uniform grid.
    /// </summary>
    struct UniformGrid {
        /// defines
        using uint = unsigned int;
        /// data
        int idim{ 0 };
        int jdim{ 0 };
        int kdim{ 0 };
        float x0{ 0 };
        float y0{ 0 };
        float z0{ 0 };
        float dx{ 0 };
        float dy{ 0 };
        float dz{ 0 };
        float* d_scalar{ nullptr };
        std::shared_ptr<float> d_scalar_{ nullptr };
        /// <summary>
        /// Contructors
        /// </summary>
        /// <returns></returns>
        __host__ UniformGrid() {}
        /// <summary>
        /// Copy Constructor
        /// </summary>
        /// <param name="ug">UniformGrid</param>
        /// <returns></returns>
        __host__ UniformGrid(const UniformGrid& ug)
        {
            idim = ug.idim;
            jdim = ug.jdim;
            kdim = ug.kdim;
            x0 = ug.x0;
            y0 = ug.y0;
            z0 = ug.z0;
            dx = ug.dx;
            dy = ug.dy;
            dz = ug.dz;
            d_scalar = ug.d_scalar;
            d_scalar_ = ug.d_scalar_;
        }
        /// <summary>
        /// Construct a uniform grid object
        /// </summary>
        /// <param name="nx">x-size</param>
        /// <param name="ny">y-size</param>
        /// <param name="nz">z-size</param>
        /// <param name="xmin">min x-coord</param>
        /// <param name="ymin">min y-coord</param>
        /// <param name="zmin">min z-coord</param>
        /// <param name="dx">spacing in x-coord</param>
        /// <param name="dy">spacing in y-coord</param>
        /// <param name="dz">spacing in z-coord</param>
        /// <param name="data">scalar data</param>
        /// <returns></returns>
        __host__ UniformGrid(const int nx, const int ny, const int nz,
            const float xmin, const float ymin, const float zmin,
            const float dx, const float dy, const float dz,
            const float* data)
        {
            this->idim = nx;
            this->jdim = ny;
            this->kdim = nz;
            x0 = xmin;
            y0 = ymin;
            z0 = zmin;
            this->dx = dx;
            this->dy = dz;
            this->dz = dz;
            cudaMalloc(&d_scalar, nx * ny * nz * sizeof(float));
            cudaMemcpy((void*)d_scalar, (void*)data, nx * ny * nz * sizeof(float), cudaMemcpyHostToDevice);
            cudaCheckError();
            d_scalar_.reset(d_scalar, cudaFree);
        }
        /// <summary>
        /// Destructor
        /// </summary>
        /// <returns></returns>
        __host__ ~UniformGrid()
        {
            idim = 0;
            jdim = 0;
            kdim = 0;
            x0 = 0;
            y0 = 0;
            z0 = 0;
            dx = 0;
            dy = 0;
            dz = 0;
            d_scalar_.reset();
            d_scalar = nullptr;
        }
        /// <summary>
        /// Assignment operator
        /// </summary>
        /// <param name="ug">uniform grid object to be copied</param>
        /// <returns></returns>
        __host__ UniformGrid& operator=(const UniformGrid& ug)
        {
            idim = ug.idim;
            jdim = ug.jdim;
            kdim = ug.kdim;
            x0 = ug.x0;
            y0 = ug.y0;
            z0 = ug.z0;
            dx = ug.dx;
            dy = ug.dy;
            dz = ug.dz;
            d_scalar = ug.d_scalar;
            d_scalar_ = ug.d_scalar_;
            return *this;
        }
        /// <summary>
        /// reset volume data
        /// </summary>
        void clear()
        {
            idim = 0;
            jdim = 0;
            kdim = 0;
            x0 = 0;
            y0 = 0;
            z0 = 0;
            dx = 0;
            dy = 0;
            dz = 0;
            d_scalar = nullptr;
            d_scalar_.reset();
        }
        /// <summary>
        /// Set scalar data. Copy scalar data to device, other data
        /// of uniform grid was already set.
        /// </summary>
        /// <param name="data"></param>
        /// <returns></returns>
        __host__ void setScalar(const float* data)
        {
            cudaMalloc(&d_scalar, t_size() * sizeof(float));
            cudaMemcpy((void*)d_scalar, (void*)data, t_size() * sizeof(float), cudaMemcpyHostToDevice);
            cudaCheckError();
            d_scalar_.reset(d_scalar, cudaFree);
        }
        /// <summary>
        /// Compute global index
        /// </summary>
        /// <param name="i">i-index</param>
        /// <param name="j">j-index</param>
        /// <param name="k">k-index</param>
        /// <returns>global index</returns>
        __host__ __device__ int gl_index(const int i, const int j, const int k) {
            return (k * jdim * idim + j * idim + i);
        }
        /// <summary>
        /// Compute i-index from global index
        /// </summary>
        /// <param name="gl_index">global index</param>
        /// <returns>i-index</returns>
        __host__ __device__ int i_index(const int gl_index) {
            return (gl_index % idim);
        }
        /// <summary>
        /// Compute j-index from global index
        /// </summary>
        /// <param name="gl_index">global index</param>
        /// <returns>j-index</returns>
        __host__ __device__ int j_index(const int gl_index) {
            return ((gl_index / idim) % jdim);
        }
        /// <summary>
        /// Compute k-index from global index
        /// </summary>
        /// <param name="gl_index">global index</param>
        /// <returns>k-index</returns>
        __host__ __device__ int k_index(const int gl_index) {
            return (gl_index / (idim * jdim));
        }
        /// <summary>
        /// Total number of vertices in uniform grid, i.e. total size of scalar data
        /// </summary>
        /// <returns>total size of scalar field</returns>
        __host__ __device__ int t_size() { return idim * jdim * kdim; }
        /// <summary>
        /// Set grid size
        /// </summary>
        /// <param name="dims">grid size</param>
        /// <returns></returns>
        __host__ __device__ void size(int dims[3])
        {
            dims[0] = idim;
            dims[1] = jdim;
            dims[2] = kdim;
        }
        /// <summary>
        /// Grid size in x-dimension
        /// </summary>
        /// <returns>i-size</returns>
        __host__ __device__ int i_size() { return idim; }
        /// <summary>
        /// Grid size in y-dimension
        /// </summary>
        /// <returns>j-size</returns>
        __host__ __device__ int j_size() { return jdim; }
        /// <summary>
        /// Grid size in z-dimension
        /// </summary>
        /// <returns>k-size</returns>
        __host__ __device__ int k_size() { return kdim; }
        /// <summary>
        /// Set grid size
        /// </summary>
        /// <param name="x_size">size in x-dimension</param>
        /// <param name="y_size">size in y-dimension</param>
        /// <param name="z_size">size in z-dimension</param>
        /// <returns></returns>
        __host__ void size(const int x_size, const int y_size, const int z_size)
        {
            idim = x_size;
            jdim = y_size;
            kdim = z_size;
        }
        /// <summary>
        /// compute unique edge global index.
        /// </summary>
        /// <param name="e">local edge index</param>
        /// <param name="i_idx">i-index of cell</param>
        /// <param name="j_idx">j-index of cell</param>
        /// <param name="k_idx">k-index of cell</param>
        /// <returns></returns>
        __host__ __device__ int e_glIndex(const int e, const int i_idx, const int j_idx, const int k_idx)
        {
            const unsigned long long gei_pattern_ = 670526590282893600ull;
            const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
            const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
            const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
            const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
            return (3 * gl_index(i, j, k) + offs);
        }
        /// <summary>
        /// Set coordinates oriding of uniform grid, required to compute
        /// vertex coordinates in physical space
        /// </summary>
        /// <param name="x">min x-coordinate</param>
        /// <param name="y">min y-coordinate</param>
        /// <param name="z">min z-coordinate</param>
        /// <returns></returns>
        __host__ void origin(const float x, const float y, const float z)
        {
            x0 = x;
            y0 = y;
            z0 = z;
        }
        /// <summary>
        /// Set grid spacing in all three dimensions
        /// </summary>
        /// <param name="x">spacing in x-dimension</param>
        /// <param name="y">spacing in y-dimension</param>
        /// <param name="z">spacing in z-dimension</param>
        /// <returns></returns>
        __host__ void spacing(const float x, const float y, const float z)
        {
            dx = x;
            dy = y;
            dz = z;
        }
        /// evaluate volume data at grid vertex (i,j,k)-index
        __device__ float operator() (const int i, const int j, const int k) { return d_scalar[k * jdim * idim + j * idim + i]; }
        /// evaluate volume data at grid vertex
        __device__ float evaluate(const int i, const int j, const int k) { return d_scalar[k * jdim * idim + j * idim + i]; }
        /// trilinear interpolation of position given local coordinates and cell vertices
        __device__ float3 trilinear(const float3 p[8], const float u, const float v, const float w)
        {
            float3 po;
            po.x = (1 - w) * ((1 - v) * (p[0].x + u * (p[1].x - p[0].x)) + v * (p[2].x + u * (p[3].x - p[2].x))) + w * ((1 - v) * (p[4].x + u * (p[5].x - p[4].x)) + v * (p[6].x + u * (p[7].x - p[6].x)));
            po.y = (1 - w) * ((1 - v) * (p[0].y + u * (p[1].y - p[0].y)) + v * (p[2].y + u * (p[3].y - p[2].y))) + w * ((1 - v) * (p[4].y + u * (p[5].y - p[4].y)) + v * (p[6].y + u * (p[7].y - p[6].y)));
            po.z = (1 - w) * ((1 - v) * (p[0].z + u * (p[1].z - p[0].z)) + v * (p[2].z + u * (p[3].z - p[2].z))) + w * ((1 - v) * (p[4].z + u * (p[5].z - p[4].z)) + v * (p[6].z + u * (p[7].z - p[6].z)));
            return po;
        }

        /// trilinear interpolation of scalar given local coordinates and scalar values at the vertices
        __device__ float trilinear(const float p[8], const float u, const float v, const float w)
        {
            return (1 - w) * ((1 - v) * (p[0] + u * (p[1] - p[0])) + v * (p[2] + u * (p[3] - p[2]))) + w * ((1 - v) * (p[4] + u * (p[5] - p[4])) + v * (p[6] + u * (p[7] - p[6])));
        }
        /// compute the gradient of the scalar field with central difference
        __device__ void gradient(float3 n[8], const float s[8], const int i, const int j, const int k)
        {
            auto index = [](const int dim, const int ii) { return (ii < 0 ? 0 : ii >= dim ? (dim - 1) : ii); };
            // vertex 0
            n[0].x = 0.5f * (s[1] - evaluate(index(idim, i - 1), j, k)) / dx;
            n[0].y = 0.5f * (s[2] - evaluate(i, index(jdim, j - 1), k)) / dy;
            n[0].z = 0.5f * (s[4] - evaluate(i, j, index(kdim, k - 1))) / dz;

            // vertex 1
            n[1].x = 0.5f * (evaluate(index(idim, i + 2), j, k) - s[0]) / dx;
            n[1].y = 0.5f * (s[3] - evaluate(i + 1, index(jdim, j - 1), k)) / dy;
            n[1].z = 0.5f * (s[5] - evaluate(i + 1, j, index(kdim, k - 1))) / dz;

            // vertex 2
            n[2].x = 0.5f * (s[3] - evaluate(index(idim, i - 1), j + 1, k)) / dx;
            n[2].y = 0.5f * (evaluate(i, index(jdim, j + 2), k) - s[0]) / dy;
            n[2].z = 0.5f * (s[6] - evaluate(i, j + 1, index(kdim, k - 1))) / dz;

            // vertex 3
            n[3].x = 0.5f * (evaluate(index(idim, i + 2), j + 1, k) - s[2]) / dx;
            n[3].y = 0.5f * (evaluate(i + 1, index(jdim, j + 2), k) - s[1]) / dy;
            n[3].z = 0.5f * (s[7] - evaluate(i + 1, j + 1, index(kdim, k - 1))) / dz;

            // vertex 4
            n[4].x = 0.5f * (s[5] - evaluate(index(idim, i - 1), j, k + 1)) / dx;
            n[4].y = 0.5f * (s[6] - evaluate(i, index(jdim, j - 1), k + 1)) / dy;
            n[4].z = 0.5f * (evaluate(i, j, index(kdim, k + 2)) - s[0]) / dz;

            // vertex 5
            n[5].x = 0.5f * (evaluate(index(idim, i + 2), j, k + 1) - s[4]) / dx;
            n[5].y = 0.5f * (s[7] - evaluate(i + 1, index(jdim, j - 1), k + 1)) / dy;
            n[5].z = 0.5f * (evaluate(i + 1, j, index(kdim, k + 2)) - s[1]) / dz;

            // vertex 6
            n[6].x = 0.5f * (s[7] - evaluate(index(idim, i - 1), j + 1, k + 1)) / dx;
            n[6].y = 0.5f * (evaluate(i, index(jdim, j + 2), k + 1) - s[4]) / dy;
            n[6].z = 0.5f * (evaluate(i, j + 1, index(kdim, k + 2)) - s[2]) / dz;

            // vertex 7
            n[7].x = 0.5f * (evaluate(index(idim, i + 2), j + 1, k + 1) - s[6]) / dx;
            n[7].y = 0.5f * (evaluate(i + 1, index(jdim, j + 2), k + 1) - s[5]) / dy;
            n[7].z = 0.5f * (evaluate(i + 1, j + 1, index(kdim, k + 2)) - s[3]) / dz;
        }

        /// compute the vertices of a cell in the uniform grid by index
        __device__ void cell_vertices(float3 v[8], const int i, const int j, const int k)
        {
            v[0].x = x0 + i * dx;
            v[0].y = y0 + j * dy;
            v[0].z = z0 + k * dz;

            v[1].x = v[0].x + dx;
            v[1].y = v[0].y;
            v[1].z = v[0].z;

            v[2].x = v[0].x;
            v[2].y = v[0].y + dy;
            v[2].z = v[0].z;

            v[3].x = v[0].x + dx;
            v[3].y = v[0].y + dy;
            v[3].z = v[0].z;

            v[4].x = v[0].x;
            v[4].y = v[0].y;
            v[4].z = v[0].z + dz;

            v[5].x = v[0].x + dx;
            v[5].y = v[0].y;
            v[5].z = v[0].z + dz;

            v[6].x = v[0].x;
            v[6].y = v[0].y + dy;
            v[6].z = v[0].z + dz;

            v[7].x = v[0].x + dx;
            v[7].y = v[0].y + dy;
            v[7].z = v[0].z + dz;
        }
        /// compute vertex coordinate by vertex index
        __device__ float3 cell_vertex(const int i, const int j, const int k)
        {
            return make_float3(x0 + i * dx, y0 + j * dy, z0 + k * dz);
        }
        /// compute cell index from vertex position
        __device__ int3 cell_index(const float3& p)
        {
            int3 index;
            index.x = static_cast<int>((p.x - x0) / dx);
            index.y = static_cast<int>((p.y - y0) / dy);
            index.z = static_cast<int>((p.z - z0) / dz);
            return index;
        }
        /// <summary>
        /// Generate uniform grid for predefined scalar functions
        /// </summary>
        enum class SurfaceCase: int {
            Sphere = 0, Torus = 1, TwoHoledTorus = 2, FourHoledTorus = 3,
            GenusTwo = 4,
            iWP = 5,
            pwHybrid = 6,
            neovius = 7,
            Goursat = 8,
            SteinerRoman = 9
        };
        /// <summary>
        /// Given a predefinde function and size of the grid, generate scalar data
        /// </summary>
        /// <param name="dim"></param>
        /// <param name="sc"></param>
        /// <returns></returns>
        __host__ void generateVolume(const std::array<int, 3>& dim, SurfaceCase sc);
        /// <summary>
        /// Read scalar field on a uniform grid from file. Scalar values might
        /// be given as unsigned short (CT- or MRI-data) or as floats.
        /// </summary>
        template<typename T>
        __host__ void readDataFromFile(const std::string& i_file)
        {
            std::ifstream ifile;
            ifile.open(i_file, std::ios::binary);
            if (!ifile.is_open()) {
                exit(1);
            }
            ifile.read(reinterpret_cast<char*>(&idim), sizeof(int));
            ifile.read(reinterpret_cast<char*>(&jdim), sizeof(int));
            ifile.read(reinterpret_cast<char*>(&kdim), sizeof(int));
            ifile.read(reinterpret_cast<char*>(&dx), sizeof(float));
            ifile.read(reinterpret_cast<char*>(&dy), sizeof(float));
            ifile.read(reinterpret_cast<char*>(&dz), sizeof(float));
            x0 = 0;
            y0 = 0;
            z0 = 0;

            size_t size_ = t_size(); // static_cast<size_t>(m_nx)* static_cast<size_t>(m_ny)* static_cast<size_t>(m_nz);
            std::vector<T> t_buff(size_);
            ifile.read(reinterpret_cast<char*>(t_buff.data()), size_ * sizeof(T));
            ifile.close();
            float* h_volume = new float[size_];
            int pos{ 0 };
            for (auto v : t_buff)
            {
                h_volume[pos] = static_cast<float>(v);
                pos++;
            }
            //
            setScalar(h_volume);
            delete[] h_volume;
        }
    };
}// namespace
