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
    /// Edge data structure. An edge consist of the two indices of the 
    /// end vertices.
    /// </summary>
    struct Edges {
        /// size of buffer
        int a_size{ 0 };
        /// index buffer
        int2* edges{ nullptr };
        std::shared_ptr<int2> edges_{ nullptr };
        /// atomic counter
        int* t_size{ nullptr };
        std::shared_ptr<int> t_size_{ nullptr };
        /// total number of edges
        int nr_e{ 0 };
        /// constructors
        __host__ Edges() {}
        __host__ Edges(const int sz) : a_size{ sz }, nr_e{ 0 }
        {
            cudaMalloc(&edges, sz * sizeof(int2));
            cudaMalloc(&t_size, sizeof(int));
            cudaCheckError();
            cudaMemset(t_size, 0, sizeof(int));
            cudaCheckError();
            edges_.reset(edges, cudaFree);
            t_size_.reset(t_size, cudaFree);
        }
        /// destructor
        __host__ ~Edges()
        {
            a_size = 0;
            edges_.reset();
            t_size_.reset();
            edges = nullptr;
            t_size = nullptr;
        }
        /// get size of buffer
        __host__ int capacity() { return a_size; }
        /// get number of edges
        __host__ int size()
        {
            cudaMemcpy(&nr_e, t_size, sizeof(int), cudaMemcpyDeviceToHost);
            return nr_e;
        }
        /// set default value to atomic counter
        __host__ void initAtomicCounter()
        {
            cudaMemset(t_size, 0, sizeof(int));
            nr_e = 0;
        }
        /// add an edge
        __device__ int addEdge(const int v0, const int v1)
        {
            const int pos = atomicAdd(t_size, 1);
            edges[pos].x = v0;
            edges[pos].y = v1;
            return pos;
        }

        /// copy data
        __host__ void copy(Edges& e)
        {
            nr_e = e.size();
            if (nr_e > a_size)
            {
                cudaMalloc(&edges, nr_e * sizeof(int2));
                cudaCheckError();
                edges_.reset(edges, cudaFree);
            }
            cudaMemcpy((void*)edges, e.edges, nr_e * sizeof(int2), cudaMemcpyDeviceToDevice);
            cudaMemcpy((void*)t_size, e.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
            cudaCheckError();
        }
        /// increase size if necessary, data get lost
        __host__ void resize(const int sz)
        {
            if (sz > a_size)
            {
                cudaMalloc(&edges, sz * sizeof(int2));
                cudaMalloc(&t_size, sizeof(int));
                cudaCheckError();
                cudaMemset(t_size, 0, sizeof(int));
                cudaCheckError();
                edges_.reset(edges, cudaFree);
                t_size_.reset(t_size, cudaFree);
            }
        }
    };
} // namespace p_mc