#pragma once

// C++ libs
#include <memory>
#include <vector>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"

namespace p_mc {
    struct ValenceHashTable {
        using ulong = unsigned long long;
        int t_size{ 0 };
        /// keys
        ulong* key{ nullptr };
        std::shared_ptr<ulong> key_{ nullptr };
        /// end vertex indices
        int2* edge{ nullptr };
        std::shared_ptr<int2> edge_{ nullptr };

        /// constructors
        __host__ ValenceHashTable() {}
        __host__ ValenceHashTable(const int sz) : t_size{ sz }
        {
            cudaMalloc(&key, t_size * sizeof(ulong));
            cudaMalloc(&edge, t_size * sizeof(int2));
            cudaCheckError();
            key_.reset(key, cudaFree);
            edge_.reset(edge, cudaFree);
        }
        /// destructor
        __host__ ~ValenceHashTable()
        {
            t_size = 0;
            key_.reset();
            key = nullptr;
            edge_.reset();
            edge = nullptr;
        }
        /// get buffer size
        __host__ __device__ int size() { return t_size; }
        /// change buffer size
        __host__ void resize(const int sz)
        {
            if (sz != t_size)
            {
                cudaMalloc(&key, sz * sizeof(ulong));
                cudaCheckError();
                key_.reset(key, cudaFree);
                cudaMalloc(&edge, sz * sizeof(int2));
                cudaCheckError();
                edge_.reset(edge, cudaFree);
                t_size = sz;
            }
        }
        /// init buffers
        __device__ void init(const int pos)
        {
            key[pos] = EMPTY_BUCKET_64;
            edge[pos].x = INVALID_INDEX;
            edge[pos].y = INVALID_INDEX;
        }
        /// <summary>
        ///  check if bucket is empty
        /// </summary>
        /// <param name="pos">address in key array</param>
        /// <returns>true if bucket was not used</returns>
        __device__ bool empty(const int pos)
        {
            return key[pos] == EMPTY_BUCKET_64;
        }
        /// add edge to table
        __device__ bool addEdge(const int v0, const int v1)
        {
            ulong k = setKey(v0, v1);
            int h{ hash1(k) };
            int e{ 1 };
            for (int i = 0; i < 128; i++)
            {
                const ulong old = atomicCAS(&key[h], EMPTY_BUCKET_64, k);
                if (old == EMPTY_BUCKET_64)
                {
                    edge[h].x = v0;
                    edge[h].y = v1;
                    return true;
                }
                else if (old == k)
                {
                    return true;
                }
                else
                {
                    // step quadratic
                    h = (h + e * e) % t_size;
                    e++;
                }
            }
            return false;
        }
        /// <summary>
        /// set 64bit key from two 32bit integers, the smallest integer
        /// is saved in the last 32 bits.
        /// </summary>
        /// <param name="v0">vertex index</param>
        /// <param name="v1">vertex index</param>
        /// <returns>64bit key</returns>
        __device__ ulong setKey(const int v0, const int v1)
        {
            if (v0 < v1)
                return (static_cast<ulong>(v0) << 32) | (v1 & 0xffffffffL);
            else
                return (static_cast<ulong>(v1) << 32) | (v0 & 0xffffffffL);
        }
        /// simple hash function for 64bit key
        __device__ int hash1(const ulong k)
        {
            return static_cast<int>(k % t_size);
        }
        /// murmur like hash function
        __device__ int hash2(const int v0, const int v1)
        {
            ulong h = setKey(v0, v1);
            h ^= h >> 33;
            h *= 0xff51afd7ed558ccdL;
            h ^= h >> 33;
            h *= 0xc4ceb9fe1a85ec53L;
            h ^= h >> 33;
            return static_cast<int>(h % t_size);
        }
        /// 64bit hash function
        __device__ ulong hash64shift(unsigned long long key)
        {
            key = (~key) + (key << 21); // key = (key << 21) - key - 1;
            key = key ^ (key >> 24);
            key = (key + (key << 3)) + (key << 8); // key * 265
            key = key ^ (key >> 14);
            key = (key + (key << 2)) + (key << 4); // key * 21
            key = key ^ (key >> 28);
            key = key + (key << 31);
            return key;
        }
        /// Input 64bit key and return 32bit address
        __device__ int hash3(ulong key)
        {
            key = (~key) + (key << 18); // key = (key << 18) - key - 1;
            key = key ^ (key >> 31);
            key = key * 21; // key = (key + (key << 2)) + (key << 4);
            key = key ^ (key >> 11);
            key = key + (key << 6);
            key = key ^ (key >> 22);
            return static_cast<int>(key % t_size);
        }
        /// access twins
        __device__ int v0(const int pos) { return edge[pos].x; }
        __device__ int v1(const int pos) { return edge[pos].y; }
    };
}// namespace p_mc
