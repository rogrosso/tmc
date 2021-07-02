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
    /// Class to keep track of vertex indices and mark vertex types required
    /// by the mesh simplification algorithms
    /// </summary>
    struct VertexMap {
        /// total number of elements
        int nr_e{ 0 };
        /// vertex valence
        int* valence{ nullptr };
        std::shared_ptr<int> valence_{ nullptr };
        /// count elements sharing this vertex
        int* count{ nullptr };
        std::shared_ptr<int> count_{ nullptr };
        /// set vertex type
        int* type{ nullptr };
        std::shared_ptr<int> type_{ nullptr };
        /// vertex to which it will be mapped
        int* twin{ nullptr };
        std::shared_ptr<int> twin_{ nullptr };
        /// mapping addres by vertex removal
        int* map_addr{ nullptr };
        std::shared_ptr<int> map_addr_{ nullptr };
        /// constructors
        __host__ VertexMap() { }
        __host__ VertexMap(const int size_) : nr_e{ size_ }
        {
            cudaMalloc(&valence, size_ * sizeof(int));
            cudaMalloc(&count, size_ * sizeof(int));
            cudaMalloc(&type, size_ * sizeof(int));
            cudaMalloc(&twin, size_ * sizeof(int));
            cudaMalloc(&map_addr, size_ * sizeof(int));
            cudaCheckError();
            valence_.reset(valence, cudaFree);
            count_.reset(count, cudaFree);
            type_.reset(type, cudaFree);
            twin_.reset(twin, cudaFree);
            map_addr_.reset(map_addr, cudaFree);
        }
        /// desctructor
        //__host__ __device__ ~VertexMap()
        __host__ ~VertexMap()
        {
            nr_e = 0;
            map_addr_.reset();
            cudaCheckError();
            valence_.reset();
            cudaCheckError();
            count_.reset();
            cudaCheckError();
            type_.reset();
            cudaCheckError();
            twin_.reset();
            cudaCheckError();
            valence = nullptr;
            count = nullptr;
            type = nullptr;
            twin = nullptr;
            map_addr = nullptr;
        }
        /// size of buffers
        __host__ int capacity() { return nr_e; }
        /// total number of elements
        __host__ __device__ int size() { return nr_e; }
        /// change number of elments, if necessary resize buffers
        __host__ void resize(const int size_)
        {
            if (size_ != nr_e)
            {
                nr_e = size_;
                cudaMalloc(&valence, size_ * sizeof(int));
                valence_.reset(valence, cudaFree);
                cudaCheckError();
                cudaMalloc(&count, size_ * sizeof(int));
                count_.reset(count, cudaFree);
                cudaCheckError();
                cudaMalloc(&type, size_ * sizeof(int));
                type_.reset(type, cudaFree);
                cudaCheckError();
                cudaMalloc(&twin, size_ * sizeof(int));
                twin_.reset(type, cudaFree);
                cudaCheckError();
                cudaMalloc(&map_addr, size_ * sizeof(int));
                map_addr_.reset(map_addr, cudaFree);
                cudaCheckError();
            }
        }
        /// set default values
        __device__ void init(const int pos)
        {
            valence[pos] = 0;
            count[pos] = 0;
            type[pos] = P_NEUTRAL; // no special type
            twin[pos] = INVALID_INDEX; // invalid index
            map_addr[pos] = INVALID_INDEX; // invalid index
        }
    };
} // namespace p_mc
