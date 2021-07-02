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
    struct EdgeHashTable {
        using ulong = unsigned long long;
        int t_size{ 0 };
        /// keys
        ulong* key{ nullptr };
        std::shared_ptr<ulong> key_{ nullptr };
        /// twins indices
        int2* edge{ nullptr };
        std::shared_ptr<int2> edge_{ nullptr };
        /// non manifold cases up to four faces can share an edge
        int* face{ nullptr };
        std::shared_ptr<int> face_;
        /// count how many faces share this edge
        int* nr_face{ nullptr };
        std::shared_ptr<int> nr_face_{ nullptr };
        /// constructors
        __host__ EdgeHashTable() {}
        __host__ EdgeHashTable(const int sz) : t_size{ sz }
        {
            cudaMalloc(&key, t_size * sizeof(ulong));
            cudaMalloc(&edge, t_size * sizeof(int2));
            cudaMalloc(&face, 4 * t_size * sizeof(int));
            cudaMalloc(&nr_face, t_size * sizeof(int));
            cudaCheckError();
            key_.reset(key, cudaFree);
            edge_.reset(edge, cudaFree);
            face_.reset(face, cudaFree);
            nr_face_.reset(nr_face, cudaFree);
        }
        /// destructor
        __host__ ~EdgeHashTable()
        {
            t_size = 0;
            key_.reset();
            key = nullptr;
            edge_.reset();
            edge = nullptr;
            face_.reset();
            face = nullptr;
            nr_face_.reset();
            nr_face = nullptr;
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
                cudaMalloc(&face, 4 * sz * sizeof(int));
                cudaCheckError();
                face_.reset(face, cudaFree);
                cudaMalloc(&nr_face, sz * sizeof(int));
                cudaCheckError();
                nr_face_.reset(nr_face, cudaFree);
                t_size = sz;
            }
        }
        /// init buffers
        __device__ void init(const int pos)
        {
            key[pos] = EMPTY_BUCKET_64;
            edge[pos].x = INVALID_INDEX;
            edge[pos].y = INVALID_INDEX;
            face[4*pos] = INVALID_INDEX;
            face[4*pos + 1] = INVALID_INDEX;
            face[4*pos + 2] = INVALID_INDEX;
            face[4*pos + 3] = INVALID_INDEX;
            nr_face[pos] = 0;
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
        /// Save edge into hash table and simultaneously save the index of
        /// the face sharing this edge. It might give one face for a boundary
        /// edge, two for inner edges and four faces shearing the edge, if it
        /// is the non-manifold case.
        /// </summary>
        /// <param name="v0">vertex index</param>
        /// <param name="v1">vertex index</param>
        /// <param name="f">Index of face sharing edge</param>
        /// <returns>returns true, if edge was added, otherwise returns false</returns>
        __device__ bool addEdge(const int v0, const int v1, const int f)
        {
            const ulong k = setKey(v0, v1);
            int h{ hash1(k) };
            int e{ 1 };
            for (int i = 0; i < 128; i++)
            {
                const ulong old = atomicCAS(&key[h], EMPTY_BUCKET_64, k);
                if (old == EMPTY_BUCKET_64)
                { // this bucket is empty, first use
                    edge[h].x = v0;
                    edge[h].y = v1;
                    const int pos = atomicAdd(&nr_face[h], 1);
                    face[4*h + pos] = f;
                    return true;
                }
                else if (old == k)
                { // the bucket contains already data, set face sharing this edge
                    const int pos = atomicAdd(&nr_face[h],1);
                    face[4*h + pos] = f;
                    return true;
                }
                else
                {   // the bucket is used by other edge
                    // step quadratic to find an empty bucket
                    h = (h + e * e) % t_size;
                    e++;
                }
            }
            return false;
        }
        /// <summary>
        /// set 64bit key from two 32bit integers, the smalest interger
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
        /// compute number of faces incident to edge v0,v1
        __device__ int nrFaces(const int pos)
        {
            return nr_face[pos];
        }
        /// <summary>
        /// Return index of face sharing this edge
        /// </summary>
        /// <param name="pos">index of bucket in hash table</param>
        /// <returns>face index</returns>
        __device__ int f0(const int pos) { return face[4*pos]; }
        __device__ int f1(const int pos) { return face[4*pos + 1]; }
        __device__ int f2(const int pos) { return face[4*pos + 2]; }
        __device__ int f3(const int pos) { return face[4*pos + 3]; }
        /// <summary>
        /// Host method to read edges from hash table in device.
        /// Each edge consists of the two indices of the end vertices.
        /// </summary>
        /// <param name="e"></param>
        /// <returns></returns>
        __host__ void getEdges(std::vector<int2>& e)
        {
            e.resize(t_size);
            cudaMemcpy(e.data(), edge, t_size * sizeof(int2), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
        /// <summary>
        /// Host function to read face indices from device, there might be up to four
        /// faces for each edge. Empty addresses has an invalid address, -1 in this case
        /// </summary>
        /// <param name="f">vector to store face ids.</param>
        /// <returns></returns>
        __host__ void getFaces(std::vector<int>& f)
        {
            f.resize(4 * t_size);
            cudaMemcpy(f.data(), face, 4 * t_size * sizeof(int), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
        /// <summary>
        /// Returns the number of faces sharing an edge
        /// </summary>
        /// <param name="f">vector to store data from device</param>
        /// <returns></returns>
        __host__ void getNrFaces(std::vector<int>& f)
        {
            f.resize(t_size);
            cudaMemcpy((void*)f.data(), nr_face, t_size * sizeof(int), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
    };
}// namespace p_mc
