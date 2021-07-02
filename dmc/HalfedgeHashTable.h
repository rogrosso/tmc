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
    /// <summary>
    /// A hash table to store halfedges, required to
    /// construct the halfedge data structure out of
    /// a shared vertex representation of a mesh.
    /// </summary>
    struct HalfedgeHashTable {
        using ulong = unsigned long long;
        /// buffers siz
        int t_size{ 0 };
        /// keys
        ulong* key{ nullptr };
        std::shared_ptr<ulong> key_{ nullptr };
        /// twins indices
        int* twin{ nullptr };
        std::shared_ptr<int> twin_{ nullptr };
        /// indices of faces sharing the same edge
        //int* face{ nullptr };
        //std::shared_ptr<int> face_{ nullptr };
        /// number of faces sharing
        int* nr_face{ nullptr };
        std::shared_ptr<int> nr_face_{ nullptr };
        /// constructors
        __host__ HalfedgeHashTable() {}
        __host__ HalfedgeHashTable(const int sz) : t_size{ sz }
        {
            cudaMalloc(&key, t_size * sizeof(ulong));
            cudaMalloc(&twin, 4 * t_size * sizeof(int));
            //cudaMalloc(&face, 4 * t_size * sizeof(int));
            cudaMalloc(&nr_face, t_size * sizeof(int));
            cudaCheckError();
            key_.reset(key, cudaFree);
            twin_ = std::shared_ptr<int>(twin, cudaFree);
            //face_.reset(face, cudaFree);
            nr_face_.reset(nr_face, cudaFree);
        }
        /// destructor
        __host__ ~HalfedgeHashTable()
        {
            t_size = 0;
            key_.reset();
            key = nullptr;
            twin_.reset();
            twin = nullptr;
            nr_face_.reset();
            nr_face = nullptr;
        }
        /// get buffer size
        __host__ __device__ int size() { return t_size; }
        /// change buffer size, data get lost.
        __host__ void resize(const int sz)
        {
            if (sz > t_size)
            {
                cudaMalloc(&key, sz * sizeof(ulong));
                cudaMalloc(&twin, 4 * sz * sizeof(int));
                cudaMalloc(&nr_face, sz * sizeof(int));
                cudaCheckError();
                key_.reset(key, cudaFree);
                twin_.reset(twin, cudaFree);
                nr_face_.reset(nr_face, cudaFree);
            }
            t_size = sz;
        }
        /// init buffers
        __device__ void init(const int pos)
        {
            key[pos] = EMPTY_BUCKET_64;
            twin[4 * pos] = INVALID_INDEX;
            twin[4 * pos + 1] = INVALID_INDEX;
            twin[4 * pos + 2] = INVALID_INDEX;
            twin[4 * pos + 3] = INVALID_INDEX;
            nr_face[pos] = 0;
        }
        /// test if bucket is empty
        __device__ bool empty(const int pos)
        {
            return key[pos] == EMPTY_BUCKET_64;
        }
        /// <summary>
        /// Add halfedge twins to an edge described by two vertex indices.
        /// Due to non-manifold cases, there might be one twin for a boundary
        /// edge, two for inner edges and four twins for the non-manifold case.
        /// </summary>
        /// <param name="v0">vertex index</param>
        /// <param name="v1">vertex index</param>
        /// <param name="he">index of halfedge</param>
        /// <returns>returns true if halfedge was added to hash table</returns>
        __device__ bool addHalfedgeTwins(const int v0, const int v1, const int he)
        {
            const ulong k = setKey(v0, v1);
            int h{ hash1(k) };
            int e{ 1 };
            for (int i = 0; i < 128; i++)
            {
                const ulong old = atomicCAS(&key[h], EMPTY_BUCKET_64, k);
                if (old == EMPTY_BUCKET_64 || old == k)
                {
                    const int pos = atomicAdd(&nr_face[h], 1);
                    twin[4 * h + pos] = he;
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
        /// set 64bit key from two 32bit integers
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
        __device__ ulong hash64shift(ulong key)
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
        /// <summary>
        /// Access twins sharing this edge, there might be
        /// up to four halfedges (non-manifold case).
        /// Boundary edges and empty buckes have a non-valid
        /// index (-1 in this case)
        /// </summary>
        /// <param name="pos">address in hast table</param>
        /// <returns>halfedge index</returns>
        __device__ int t0(const int pos) { return twin[4 * pos]; }
        __device__ int t1(const int pos) { return twin[4 * pos + 1]; }
        __device__ int t2(const int pos) { return twin[4 * pos + 2]; }
        __device__ int t3(const int pos) { return twin[4 * pos + 3]; }
        /// get number of faces sharing an edge
        __device__ int nrFaces(const int pos) { return nr_face[pos]; }
        /// access hash table entries, copy data from device to host
        __host__ void getTwins(std::vector<int>& t)
        {
            t.resize(4 * t_size);
            cudaMemcpy(t.data(), twin, 4 * t_size * sizeof(int), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
    };
}// namespace p_mc
