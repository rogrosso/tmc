#pragma once

// C++ libs
#include <vector>
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"

namespace p_mc {
    /// <summary>
    /// Halfedge data structure. Halfedges are represented as
    /// follows: halfedge index is a int4:
    ///    he.x = origin vertex
    ///    he.y = face
    ///    he.z = next
    ///    he.w = twin
    /// </summary>
    struct Halfedges {
        using ulong = unsigned long long;
        /// total number of halfedges
        int nr_he{ 0 };
        /// halfedge int4:
        ///    he.x = origin vertex
        ///    he.y = face
        ///    he.z = next
        ///    he.w = twin
        int4* he_e{ nullptr };
        std::shared_ptr<int4> he_e_{ nullptr };

        /// constrcutors
        __host__ Halfedges() { }
        __host__ Halfedges(const int sz) : nr_he{ sz }
        {
            cudaMalloc(&he_e, sz * sizeof(int4));
            cudaCheckError();
            he_e_.reset(he_e, cudaFree);
        }
        /// desctructor
        __host__ ~Halfedges()
        {
            nr_he = 0;
            he_e = nullptr;
            he_e_.reset();
        }
        /// get buffer size
        __host__ int capacity() { return nr_he; }
        /// get total number of halfedges
        __host__ __device__ int size() { return nr_he; }
        /// change buffer size
        __host__ void resize(const int sz)
        {
            if (nr_he != sz)
            {
                nr_he = sz;
                cudaMalloc(&he_e, nr_he * sizeof(int4));
                cudaCheckError();
                he_e_.reset(he_e, cudaFree);
            }
        }
        /// add data to halfedge
        __device__ void addHalfedge(const int pos, const int v, const int f, const int he, const int t)
        {
            he_e[pos].x = v;
            he_e[pos].y = f;
            he_e[pos].z = he;
            he_e[pos].w = t;
        }
        /// add halfedge, set default value for twin
        /// twins are connected in a second kernel
        __device__ void addHalfedge(const int pos, const int v, const int f, const int he)
        {
            he_e[pos].x = v;
            he_e[pos].y = f;
            he_e[pos].z = he;
            he_e[pos].w = -1;
        }
        /// check if halfedge has a neighbor
        __device__ bool hasTwin(const int pos)
        {
            return he_e[pos].w == -1;
        }
        /// set twin
        __device__ void setTwin(const int pos, const int twin)
        {
            if (pos > -1 && twin > -1) he_e[pos].w = twin;
        }
        /// get twin
        __device__ int getTwin(const int pos)
        {
            return he_e[pos].w;
        }
        /// get next
        __device__ int getNext(const int pos)
        {
            return he_e[pos].z;
        }
        /// get face
        __device__ int getFace(const int pos)
        {
            return he_e[pos].y;
        }
        /// get origin vertex
        __device__ int getOrigin(const int pos)
        {
            return he_e[pos].x;
        }
        /// read halfedge data structure out of device memory
        __host__ void getHalfedges(std::vector<int4>& he)
        {
            he.resize(nr_he);
            cudaMemcpy(he.data(), he_e, nr_he * sizeof(int4), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
        /// given a start halfedge and two vertices, find the corresponding edge, vertices might not be oriented
        __device__ int findEdge(const int he0, const int v0, const int v1)
        {
            ulong key = setKey(v0, v1);
            const int hv0 = he_e[he0].x;
            const int he1 = he_e[he0].z;
            const int hv1 = he_e[he1].x;
            ulong k = setKey(hv0, hv1);
            if (k == key) return he0;
            const int he2 = he_e[he1].z;
            const int hv2 = he_e[he2].x;
            k = setKey(hv1, hv2);
            if (k == key) return he1;
            const int he3 = he_e[he2].z;
            const int hv3 = he_e[he3].x;
            k = setKey(hv2, hv3);
            if (k == key) return he2;
            k = setKey(hv3, hv0);
            if (k == key) return he3;
            return -1; // did not find the halfedge
        }
        /// given a start halfedge, find all neighbors of face
        __device__ int4 findNeighbors(const int he0)
        {
            int4 n;
            // neighbor 0
            n.x = he_e[he_e[he0].w].y;
            // neighbor 1
            int next = he_e[he0].z;
            n.y = he_e[he_e[next].w].y;
            // neighbor 2
            next = he_e[next].z;
            n.z = he_e[he_e[next].w].y;
            // neighbor 3
            next = he_e[next].z;
            n.w = he_e[he_e[next].w].y;
            return n;
        }
        /// given a start halfedge, find all face vertices
        __device__ int4 findVertices(const int he0)
        {
            int4 v;
            // vertex 0
            v.x = he_e[he0].x;
            // vertex 1
            int next = he_e[he0].z;
            v.y = he_e[next].x;
            // vertex 2
            next = he_e[next].z;
            v.y = he_e[next].x;
            // vertex 3
            next = he_e[next].z;
            v.y = he_e[next].x;
            return v;
        }
        /// construct a 64bit key by concatenating two ordered integers
        __device__ ulong setKey(const int v0, const int v1)
        {
            if (v0 < v1)
                return (static_cast<ulong>(v0) << 32) | (v1 & 0xffffffffL);
            else
                return (static_cast<ulong>(v1) << 32) | (v0 & 0xffffffffL);
        }
    };
} // namespace p_mc
