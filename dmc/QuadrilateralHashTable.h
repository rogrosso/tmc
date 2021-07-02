#pragma once

// C++ libs
#include <memory>

// CUDA stuff
// cuda stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"

namespace p_mc {
    /// <summary>
    /// Hash table used to construct a shared vertex quadrilateral mesh.
    /// </summary>
    struct QuadrilateralHashTable {
        /// buffer size
        int t_size{ 0 };
        /// buffer containing the keys
        int* keyBuff{ nullptr };
        std::shared_ptr<int> keyBuff_;
        /**
         * buffer containing vertex indices
         * four consecutive indices build a quadrilateral
         */
        int* indexBuff{ nullptr };
        std::shared_ptr<int> indexBuff_;
        /// buffer face coloring
        int* colorBuff{ nullptr };
        std::shared_ptr<int> colorBuff_;
        /// constructors
        __host__ QuadrilateralHashTable() {}
        __host__ QuadrilateralHashTable(const int sz) : t_size{ sz }
        {
            cudaMalloc(&keyBuff, sz * sizeof(int));
            cudaMalloc(&indexBuff, 4 * sz * sizeof(int));
            cudaMalloc(&colorBuff, sz * sizeof(int));
            cudaCheckError();
            keyBuff_.reset(keyBuff, cudaFree);
            indexBuff_.reset(indexBuff, cudaFree);
            colorBuff_.reset(colorBuff, cudaFree);
        }
        /// desctructor, a host and a device version
        __host__ ~QuadrilateralHashTable()
        {
            t_size = 0;
            keyBuff = nullptr;
            keyBuff_.reset();
            indexBuff = nullptr;
            indexBuff_.reset();
            colorBuff = nullptr;
            colorBuff_.reset();
        }
        /// total size of keys and vertex buffers
        __host__ __device__ int size() { return t_size; }

        /// set default values
        __device__ void init(const int pos)
        {
            keyBuff[pos] = EMPTY_BUCKET_32;
            indexBuff[4 * pos] = INVALID_INDEX;
            indexBuff[4 * pos + 1] = INVALID_INDEX;
            indexBuff[4 * pos + 2] = INVALID_INDEX;
            indexBuff[4 * pos + 3] = INVALID_INDEX;
            colorBuff[pos] = INVALID_COLOR; // this is an invalid color, there are only 24 colors
        }
        /// test if bucket is empty
        __device__ bool empty(const int pos)
        {
            return keyBuff[pos] == EMPTY_BUCKET_32;
        }
        /// add a vertex index in the vertex buffer by key
        __device__ bool addVertex(const int key, const int pos, const int v)
        {
            const int start_address = hash_function(key);
            // open hashing
            int h = start_address;
            int e = 1;
            for (int loop = 0; loop < 128; loop++) {
                const int old = atomicCAS(&keyBuff[h], EMPTY_BUCKET_32, key);
                if (old == EMPTY_BUCKET_32 || old == key)
                {
                    indexBuff[4 * h + pos] = v;
                    return true;
                }
                else {
                    // step with quadratic probing
                    h = (h + e * e) % t_size;
                    e = e + 1;
                    if (h == start_address) {
                        printf("ERROR: can't find free bucket for %d\n", key);
                        return false;
                    }
                }
            }
            return false;
        }
        /// add a vertex index in the vertex buffer by key, also add face color
        __device__ bool addVertex(const int key, const int pos, const int v, const int c)
        {
            const int start_address = hash_function(key);
            // open hashing
            int h = start_address;
            int e = 1;
            for (int loop = 0; loop < 128; loop++) {
                const int old = atomicCAS(&keyBuff[h], EMPTY_BUCKET_32, key);
                if (old == EMPTY_BUCKET_32 || old == key)
                {
                    indexBuff[4 * h + pos] = v;
                    colorBuff[h] = c;
                    return true;
                }
                else {
                    // step with quadratic probing
                    h = (h + e * e) % t_size;
                    e = e + 1;
                    if (h == start_address) {
                        printf("ERROR: can't find free bucket for %d\n", key);
                        return false;
                    }
                }
            }
            return false;
        }
        /// hash function
       /* __device__ int hash_function(const int key)
        {
            return ((3 * key) % 300000007) % t_size;
        }*/
        __device__ int hash_function(const int key)
        {
            /*const int addr = ((3ll * static_cast<unsigned long long>(key)) % 300000007ll) % t_size;
            if (addr < 0) {
                printf("Negative address\n");
            }*/
            return ((3ll * static_cast<unsigned long long>(key)) % 300000007ll) % t_size;
        }
        /// <summary>
        /// Returns indices of all four vertices constituting the quadrilateral
        /// </summary>
        /// <param name="pos">bucket address</param>
        /// <returns>indices of vertices</returns>
        __device__ int4 quadrilateral(const int pos)
        {
            int4 q;
            q.x = indexBuff[4 * pos];
            q.y = indexBuff[4 * pos + 1];
            q.z = indexBuff[4 * pos + 2];
            q.w = indexBuff[4 * pos + 3];
            return q;
        }
        __device__ int v0(const int pos) { return indexBuff[4 * pos]; }
        __device__ int v1(const int pos) { return indexBuff[4 * pos + 1]; }
        __device__ int v2(const int pos) { return indexBuff[4 * pos + 2]; }
        __device__ int v3(const int pos) { return indexBuff[4 * pos + 3]; }
        /// <summary>
        /// Returns face color
        /// </summary>
        /// <param name="pos">bucket address</param>
        /// <returns>face color</returns>
        __device__ int color(const int pos) { return colorBuff[pos]; }
    };

} // namespace p_mc
