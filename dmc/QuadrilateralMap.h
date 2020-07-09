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
    /// Class used to mark quadrilaterals, required for removal and mesh simplification.
    /// This class is no longer required, marks are directly encoded in the quadrilaterals.
    /// </summary>
    struct QuadrilateralMap {
        /// number of elements
        int nr_q{ 0 };
        /// element type
        int* type{ nullptr };
        std::shared_ptr<int> type_{ nullptr };
        /// constructors
        __host__ QuadrilateralMap() {}
        __host__ QuadrilateralMap(const int size_) : nr_q{ size_ }
        {
            cudaMalloc(&type, size_ * sizeof(int));
            cudaCheckError();
            type_.reset(type, cudaFree);
        }
        /// destructor
        __host__ ~QuadrilateralMap()
        {
            nr_q = 0;
            type = nullptr;
            type_.reset();
        }
        /// buffer size
        __host__ int capacity() { return nr_q; }
        /// number of elements
        __host__ __device__ int size() { return nr_q; }
        /// change number of elements, if necessary resize buffers
        __host__ void resize(const int size_)
        {
            if (size_ != nr_q)
            {
                nr_q = size_;
                cudaMalloc(&type, size_ * sizeof(int));
                cudaCheckError();
                type_.reset(type, cudaFree);
            }
        }
        __device__ void init(const int pos)
        {
            type[pos] = P_NEUTRAL; // no special type of element
        }
    };
}// namespace p_mc