#pragma once

// C++ libs
#include <memory>

// cuda stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"


namespace p_mc {
    /// <summary>
    /// Computes quality measure for triangles
    /// </summary>
    struct QualityMeasure {
        /// buffer size
        int a_size{ 0 };
        /// element quality
        float* q_measure{ nullptr };
        std::shared_ptr<float> q_measure_;
        /// Constructors
        __host__ QualityMeasure() {  }
        __host__ QualityMeasure(const int sz) : a_size{ sz }
        {
            cudaMalloc(&q_measure, a_size * sizeof(float));
            cudaCheckError();
            q_measure_ = std::shared_ptr<float>(q_measure, cudaFree);
        }
        // Desctructor
        __host__ ~QualityMeasure()
        {
            a_size = 0;
            q_measure = nullptr;
            q_measure_.reset();
        }
        /// buffer size
        __host__ __device__ int buffSize() { return a_size; }
        /// set and get value
        __device__ float& operator[] (const int i) { return q_measure[i]; }
        __device__ const float& operator[] (const int i) const { return q_measure[i]; }
        /// set value by address
        __device__ void set(const int pos, const float val) { q_measure[pos] = val; }
        /// compute square of twice the area of a triangle
        __device__ float area(const float3 v0, const float3 v1, const float3 v2)
        {
            const float x1 = v1.x - v0.x;
            const float x2 = v2.x - v0.x;
            const float y1 = v1.y - v0.y;
            const float y2 = v2.y - v0.y;
            const float z1 = v1.z - v0.z;
            const float z2 = v2.z - v0.z;
            const float a1 = (y1 * z2 - z1 * y2) * (y1 * z2 - z1 * y2);
            const float a2 = (x1 * z2 - z1 * x2) * (x1 * z2 - z1 * x2);
            const float a3 = (x1 * y2 - y1 * x2) * (x1 * y2 - y1 * x2);
            return sqrt(a1 + a2 + a3);
        }
        /// computes the mean ratio quality measure for triangles
        __device__ void mean_ratio(const int pos, const float3 v0, const float3  v1, const float3 v2)
        {
            const float l0 = (v1.x - v0.x) * (v1.x - v0.x) + (v1.y - v0.y) * (v1.y - v0.y) + (v1.z - v0.z) * (v1.z - v0.z);
            const float l1 = (v2.x - v1.x) * (v2.x - v1.x) + (v2.y - v1.y) * (v2.y - v1.y) + (v2.z - v1.z) * (v2.z - v1.z);
            const float l2 = (v0.x - v2.x) * (v0.x - v2.x) + (v0.y - v2.y) * (v0.y - v2.y) + (v0.z - v2.z) * (v0.z - v2.z);
            const float l = l0 + l1 + l2;
            const float a = area(v0, v1, v2);
            //return 2.0f * sqrtf(3.0f) * a / l;
            q_measure[pos] = 2.0f * sqrtf(3.0f) * a / l;
        }
    };
}// namespace p_mc
