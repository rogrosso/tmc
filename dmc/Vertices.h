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
    /// Vertices in a shared vertex and halfedge data structure.
    /// This class store coordinates and normals.
    /// </summary>
    struct Vertices {
        /// size of buffers
        int a_size;
        /// vertex buffer
        float3* vertices{ nullptr };
        std::shared_ptr<float3> vertices_{ nullptr };
        /// normal buffer
        float3* normals{ nullptr };
        std::shared_ptr<float3> normals_{ nullptr };
        /// erros for point back projection
        //float3* error{ nullptr };
        //std::shared_ptr<float3> error_{ nullptr };
        /// atomic counter
        int* t_size{ nullptr };
        std::shared_ptr<int> t_size_{ nullptr };
        /// total number of vertices
        int nr_v{ 0 };
        /// constructions
        __host__ Vertices() {}
        __host__ Vertices(const int sz) : a_size{ sz }, nr_v{ 0 }
        {
            cudaMalloc(&vertices, sz * sizeof(float3));
            cudaMalloc(&normals, sz * sizeof(float3));
            //cudaMalloc(&error, sz * sizeof(float3));
            cudaMalloc(&t_size, sizeof(int));
            cudaMemset(t_size, 0, sizeof(int));

            vertices_ = std::shared_ptr<float3>(vertices, cudaFree);
            normals_ = std::shared_ptr<float3>(normals, cudaFree);
            //error_ = std::shared_ptr<float3>(error, cudaFree);
            t_size_ = std::shared_ptr<int>(t_size, cudaFree);
        }
        /// destructor
        __host__ ~Vertices()
        {
            vertices_.reset();
            normals_.reset();
            //error_.reset();
            t_size_.reset();
            vertices = nullptr;
            normals = nullptr;
            //error = nullptr;
            t_size = nullptr;
            a_size = 0;
            nr_v = 0;
        }
        /// size of buffers: capacity in C++ vector
        __host__ int capacity() { return a_size; }
        /// number of vertices
        __host__ int size()
        {
            cudaMemcpy(&nr_v, t_size, sizeof(int), cudaMemcpyDeviceToHost);
            return nr_v;
        }
        /// init atomic counter to reuse buffers
        __host__ void initAtomicCounter()
        {
            cudaMemset(t_size, 0, sizeof(int));
            nr_v = 0;
        }
        /// copy data 
        __host__ void copy(Vertices& v)
        {
            nr_v = v.size();
            if (nr_v > a_size) {
                cudaMalloc(&vertices, nr_v * sizeof(float3));
                vertices_.reset(vertices, cudaFree); // free device memory
                cudaMalloc(&normals, nr_v * sizeof(float3));
                normals_.reset(normals, cudaFree); // free device memory 
                //cudaMalloc(&error, nr_v * sizeof(float3));
                //error_.reset(error, cudaFree); // free device memory
                cudaCheckError();
                a_size = nr_v;
            }
            cudaMemcpy((void*)vertices, v.vertices, nr_v * sizeof(float3), cudaMemcpyDeviceToDevice);
            cudaMemcpy((void*)normals, v.normals, nr_v * sizeof(float3), cudaMemcpyDeviceToDevice);
            //cudaMemcpy((void*)error, v.error, nr_v * sizeof(float3), cudaMemcpyDeviceToDevice);
            cudaMemcpy((void*)t_size, v.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
            cudaCheckError();
        }
        /// add a vertex to buffer
        __device__ int addVertex(const float3 v, const float3 n)
        {
            const int pos = atomicAdd(t_size, 1);
            vertices[pos] = v;
            normals[pos] = n;
            return pos;
        }
        /// add a vertex to buffer, consider the approximation error
        /*__device__ int addVertex(const float3 v, const float3 n, const float3 e)
        {
            const int pos = atomicAdd(t_size, 1);
            vertices[pos] = v;
            normals[pos] = n;
            error[pos] = e;
            return pos;
        }*/
        /// add info given the position in array
        //__device__ void addInfo(const int pos, const float3 info)
        //{
            //error[pos] = info;
        //}
    };
} // namespace p_mc
