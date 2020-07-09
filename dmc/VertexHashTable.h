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
    /// Hash table to construct a shared vertex triangle mesh using the 
    /// standard Marching Cubes algorithms. Triangles are also constructed
    /// for rendering purposes.
    /// </summary>
    struct VertexHashTable {
        using uint = unsigned int;
        /// size of buffers
        int t_size{ 0 };
        /// key buffer
        int* key{ nullptr };
        std::shared_ptr<int> key_{ nullptr };
        /// index buffer
        int* index{ nullptr };
        std::shared_ptr<int> index_{ nullptr };
        /// constructors
        __host__ VertexHashTable() {}
        __host__ VertexHashTable(const int sz) : t_size{ sz }
        {
            cudaMalloc(&key, sz * sizeof(int));
            cudaMalloc(&index, sz * sizeof(int));
            cudaCheckError();
            key_.reset(key, cudaFree);
            index_.reset(index, cudaFree);
        }
        /// <summary>
        /// Destructor
        /// </summary>
        /// <returns></returns>
        __host__ ~VertexHashTable()
        {
            t_size = 0;
            key = nullptr;
            key_.reset();
            index = nullptr;
            index_.reset();
        }
        /// get size of buffers
        __host__ __device__ int size() { return t_size; }
        /// set default values
        __device__ void init(const int pos)
        {
            key[pos] = EMPTY_BUCKET_32;
            index[pos] = INVALID_INDEX;
        }
        /// simple hash function
        __device__ int hash(const int k)
        {
            return ((3 * k) % 300000007) % t_size;
        }
        /// a bit mix hash function
        __device__ uint hash_function(uint key)
        {
            key = (key ^ 61) ^ (key >> 16);
            key = key + (key << 3);
            key = key ^ (key >> 4);
            key = key * 0x27d4eb2d;
            key = key ^ (key >> 15);
            return key;
        }
        /// <summary>
        /// Add vertex into the hash table
        /// </summary>
        /// <param name="k">key</param>
        /// <param name="addr">buffer to keep address of bucket, where vertex was added</param>
        /// <param name="pos">position in vertex buffer</param>
        /// <returns></returns>
        __device__ bool addVertex(const int k, int* addr, const int pos)
        {
            const int start_address = hash(k);
            int h = start_address;
            int e = 1;
            for (int loop = 0; loop < 128; loop++) {
                const int old = atomicCAS(&key[h], EMPTY_BUCKET_32, k);
                //if (old == EMPTY_BUCKET_32 || old == k)
                if (old == EMPTY_BUCKET_32)
                {
                    // vertex has to be created or it already exists, return address in field
                    addr[pos] = h;
                    return false;
                }
                else if (old == k)
                {
                    addr[pos] = h;
                    return true;
                }
                else
                {
                    // step with quadratic probing
                    h = (h + e * e) % t_size;
                    e = e + 1;
                    if (h == start_address) {
                        addr[pos] = -1;
                        return false;
                    }
                }
            }
            addr[pos] = -1;
            return false; // something went wrong
        }
        /// access vertex index
        __device__ int v(const int pos) { return index[pos]; }
        /// <summary>
        /// Set vertex index
        /// </summary>
        /// <param name="pos">bucket address</param>
        /// <param name="val">index to be set</param>
        /// <returns></returns>
        __device__ void set(const int pos, const int val) { index[pos] = val; }
    };
} // namespace p_mc