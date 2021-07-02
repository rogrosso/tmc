#pragma once

// C++ libs
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"

namespace p_mc {
    /// <summary>
    /// Triangles in a shared vertex data structure.
    /// </summary>
    struct Triangles {
        /// size of buffer
        int a_size{ 0 }; // size of buffer
        /// triangle index buffer
        int3* triangles{ nullptr };
        std::shared_ptr<int3> triangles_{ nullptr };
        /// atomic counter
        int* t_size{ nullptr };
        std::shared_ptr<int> t_size_{ nullptr };
        /// total number of triangles, contains same value a t_size
        int nr_t{ 0 };
        /// constructors
        __host__ Triangles() {}
        __host__ Triangles(const int sz) : a_size{ sz }, nr_t{ 0 }
        {
            cudaMalloc(&triangles, sz * sizeof(int3));
            cudaCheckError();
            cudaMalloc(&t_size, sizeof(int));
            cudaCheckError();
            cudaMemset(t_size, 0, sizeof(int));
            cudaCheckError();
            triangles_.reset(triangles, cudaFree);
            t_size_.reset(t_size, cudaFree);
        }
        /// destructor
        __host__ ~Triangles()
        {
            a_size = 0;
            triangles = nullptr;
            triangles_.reset();
            t_size = nullptr;
            t_size_.reset();
            nr_t = 0;
        }
        /// get size of triangle index buffer
        __host__ __device__ int capacity() { return a_size; }
        /// total number of triangles
        __host__ int size()
        {
            cudaMemcpy(&nr_t, t_size, sizeof(int), cudaMemcpyDeviceToHost);
            return nr_t;
        }

        /// set default value to atomic counter
        __host__ void initAtomicCounter()
        {
            cudaMemset(t_size, 0, sizeof(int));
            nr_t = 0;
        }
        /// add a triangle to index buffer
        __device__ void addTriangle(const int pos, const int v0, const int v1, const int v2)
        {
            triangles[pos].x = v0;
            triangles[pos].y = v1;
            triangles[pos].z = v2;
        }
        /// add triangle to index buffer
        __device__ int addTriangle(const int v0, const int v1, const int v2)
        {
            const int pos = atomicAdd(t_size, 1);
            this->triangles[pos] = { v0, v1, v2 };
            return pos;
        }
        /// compute minimum angle of a triangle based on cosine rule
        __device__ float minAngle(const float3 v0, const float3 v1, const float3 v2)
        {
            const float a = sqrtf((v1.x - v0.x) * (v1.x - v0.x) + (v1.y - v0.y) * (v1.y - v0.y) + (v1.z - v0.z) * (v1.z - v0.z));
            const float b = sqrtf((v2.x - v1.x) * (v2.x - v1.x) + (v2.y - v1.y) * (v2.y - v1.y) + (v2.z - v1.z) * (v2.z - v1.z));
            const float c = sqrtf((v0.x - v2.x) * (v0.x - v2.x) + (v0.y - v2.y) * (v0.y - v2.y) + (v0.z - v2.z) * (v0.z - v2.z));
            const float A = acosf((b * b + c * c - a * a) / (2 * b * c));
            const float B = acosf((a * a + c * c - b * b) / (2 * a * c));
            const float C = acosf((b * b + a * a - c * c) / (2 * b * a));

            return fminf(fminf(A, B), C);
        }
        /// copy data
        __host__ void copy(Triangles& t)
        {
            nr_t = t.size();
            if (nr_t > a_size) {
                // needs to resize buffers
                cudaMalloc(&triangles, nr_t * sizeof(int3));
                cudaCheckError();
                triangles_.reset(triangles, cudaFree);
                a_size = nr_t;
            }
            cudaMemcpy((void*)triangles, t.triangles, nr_t * sizeof(float3), cudaMemcpyDeviceToDevice);
            cudaMemcpy((void*)t_size, t.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
            cudaCheckError();
        }
        /// access vertex index
        __device__ int v0(const int pos)
        {
            return triangles[pos].x;
        }
        __device__ int v1(const int pos)
        {
            return triangles[pos].y;
        }
        __device__ int v2(const int pos)
        {
            return triangles[pos].z;
        }
    };
} // namespace p_mc
