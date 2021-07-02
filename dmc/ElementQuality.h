#pragma once

// CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// C++
#include <memory>
#include <vector>

// Project files
#include "helper_cuda.h"
#include "Vertices.h"
#include "Quadrilaterals.h"
#include "Triangles.h"

namespace p_mc {
    struct ElementQuality {
        // types
        using uint = unsigned int;
        //
        float factor{ 2 * sqrt(3.f) };
        int t_size{ 0 };
        float* quality{ nullptr };
        std::shared_ptr<float> quality_{ nullptr };
        __host__ ElementQuality() {}
        __host__ ElementQuality(const int sz) : t_size{ sz }
        {
            cudaMalloc(&quality, t_size * sizeof(float));
            cudaCheckError();
            cudaMemset(quality, 0, t_size * sizeof(float));
            quality_.reset(quality, cudaFree);
        }
        __host__ ~ElementQuality()
        {
            t_size = 0;
            quality_.reset();
            quality = nullptr;
        }
        __host__ __device__ int size() { return t_size;  }
        __device__ float3 sum(const float3& a, const float3& b)
        {
            return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
        }
        __device__ float3 diff(float3 a, float3 b)
        {
            return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
        }
        __device__ float3 cross(float3 v1, float3 v2)
        {
            return make_float3(v1.z * v2.y - v1.y * v2.z, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
        }
        __device__ float dot(float3 a, float3 b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }
        __device__ float norm(float3 v)
        {
            return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        }
        __device__ float norm2(float3 v)
        {
            return (v.x * v.x + v.y * v.y + v.z * v.z);
        }
        __device__ void quads(const int index, float3 v[4])
        {
            const float a0 = norm(cross(diff(v[1], v[0]), diff(v[3], v[0])));
            const float a1 = norm(cross(diff(v[2], v[1]), diff(v[0], v[1])));
            const float a2 = norm(cross(diff(v[3], v[2]), diff(v[1], v[2])));
            const float a3 = norm(cross(diff(v[0], v[3]), diff(v[2], v[3])));
            const float q0 = 2 * a0 / (norm2(diff(v[1], v[0])) + norm2(diff(v[3], v[0])));
            const float q1 = 2 * a1 / (norm2(diff(v[2], v[1])) + norm2(diff(v[0], v[1])));
            const float q2 = 2 * a2 / (norm2(diff(v[3], v[2])) + norm2(diff(v[1], v[2])));
            const float q3 = 2 * a3 / (norm2(diff(v[0], v[3])) + norm2(diff(v[2], v[3])));
            float q = q0 < q1 ? q0 : q1;
            q = q < q2 ? q : q2;
            q = q < q3 ? q : q3;
            quality[index] = q;
        }
        __device__ void tris(const int index, const float3 v[3])
        {
            const float a = norm(cross(diff(v[1], v[0]), diff(v[2], v[0])));
            const float l0 = norm2(diff(v[1], v[0]));
            const float l1 = norm2(diff(v[2], v[1]));
            const float l2 = norm2(diff(v[0], v[2]));
            float q{ 0 };
            if (l0 != 0 && l1 != 0 && l2 != 0)
            {
                q = factor * a / (l0 + l1 + l2);
            }
            quality[index] = q;
        }
        __host__ void getQuality(std::vector<float>& q)
        {
            q.resize(t_size);
            cudaMemcpy(q.data(), quality, t_size * sizeof(float), cudaMemcpyDeviceToHost);
        }

    };
} // namespace p_mc
