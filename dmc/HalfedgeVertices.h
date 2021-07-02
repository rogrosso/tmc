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
    /// A halfedge vertex
    /// </summary>
    struct HalfedgeVertices {
        /// total number of vertices
        int nr_v;
        /// index buffer, contains index of
        /// halfedge starting at vertex
        int* he_e{ nullptr };
        std::shared_ptr<int> he_e_{ nullptr };

        /// constructors
        __host__ HalfedgeVertices() {}
        __host__ HalfedgeVertices(const int sz) : nr_v{ sz }
        {
            cudaMalloc(&he_e, sz * sizeof(int));
            cudaCheckError();
            he_e_.reset(he_e, cudaFree);
        }
        /// destructor
        __host__ ~HalfedgeVertices()
        {
            nr_v = 0;
            he_e = nullptr;
            he_e_.reset();
        }
        /// size of buffer
        __host__ int capacity() { return nr_v; }
        /// total number of vertices
        __host__ __device__ int size() { return nr_v; }
        /// change size of buffer
        __host__ void resize(const int sz)
        {
            if (sz != nr_v)
            {
                nr_v = sz;
                cudaMalloc(&he_e, sz * sizeof(int));
                he_e_.reset(he_e, cudaFree);
            }
        }
        /// add vertex
        __device__ void addVertex(const int pos, const int e)
        {
            he_e[pos] = e;
        }
        /// read halfedge data structure out of device memory
        __host__ void getHalfedgeVertices(std::vector<int>& he_v)
        {
            he_v.resize(nr_v);
            cudaMemcpy(he_v.data(), he_e, nr_v * sizeof(int), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
    };
} // namespace p_mc
