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
    /// A quadrilateral in a shared vertex data structure. It consists of the
    /// the four indices of the vertices. It keeps track of face color and
    /// different marks required for mesh simplification.
    /// </summary>
    struct Quadrilaterals {
        using uchar = unsigned char;
        /// size of buffer
        int a_size{ 0 };
        /// quadrilateral index buffer
        int4* quadrilaterals{ nullptr };
        std::shared_ptr<int4> quadrilaterals_{ nullptr };
        /// <summary>
        /// Store face color, non-manifold case and tags for mesh simplification
        /// There are at most 24 colors, thus information is stored as follows
        ///   - first five bits for color
        ///   - sixth bit for non-manifold case
        ///   - last two bits for vertex valence pattern of face
        /// </summary>
        uchar* attributes{ nullptr };
        std::shared_ptr<uchar> attributes_;
        /// atomic counter
        int* t_size{ nullptr };
        std::shared_ptr<int> t_size_{ nullptr };
        /// total number of quadrilaterals
        int nr_q{ 0 };
        /// constructors
        __host__ Quadrilaterals() {}
        __host__ Quadrilaterals(const int sz) : a_size{ sz }, nr_q{ 0 }
        {
            cudaMalloc(&quadrilaterals, sz * sizeof(int4));
            cudaMalloc(&attributes, sz * sizeof(uchar));
            cudaMalloc(&t_size, sizeof(int));
            cudaCheckError();
            cudaMemset(t_size, 0, sizeof(int));
            cudaCheckError();
            quadrilaterals_.reset(quadrilaterals, cudaFree);
            attributes_.reset(attributes, cudaFree);
            t_size_.reset(t_size, cudaFree);
        }
        /// destructor
        __host__ ~Quadrilaterals()
        {
            a_size = 0;
            quadrilaterals_.reset();
            attributes_.reset();
            t_size_.reset();
            quadrilaterals = nullptr;
            attributes = nullptr;
            t_size = nullptr;
            nr_q = 0;
        }
        /// buffer size
        __host__ int capacity() { return a_size; }
        /// number of quadrilaterals
        __host__ int size()
        {
            cudaMemcpy(&nr_q, t_size, sizeof(int), cudaMemcpyDeviceToHost);
            return nr_q;
        }
        /// set default value to atomic counter
        __host__ void initAtomicCounter()
        {
            cudaMemset(t_size, 0, sizeof(int));
            nr_q = 0;
        }
        /// add a quadrilateral to index buffer, position if array is computed via atomicAdd
        __device__ int addQuadrilateral(const int4 q)
        {
            if (q.x >= 0 && q.y >= 0 && q.z >= 0 && q.w >= 0)
            {
                const int pos = atomicAdd(t_size, 1);
                quadrilaterals[pos] = q;
                initAttribute(pos, INVALID_COLOR);
                return pos;
            }
            else
            {
                return -1;
            }
        }
        /// add a quadrilateral to index buffer, position if array is computed via atomicAdd
        __device__ int addQuadrilateral(const int v0, const int v1, const int v2, const int v3)
        {
            if (v0 > -1 && v1 > -1 && v2 > -1 && v3 > -1)
            {
                int pos = atomicAdd(t_size, 1);
                quadrilaterals[pos] = { v0,v1,v2,v3 };
                initAttribute(pos, INVALID_COLOR);
                return pos;
            }
            else
            {
                return -1;
            }
        }
        /// add a quadrilateral to index buffer, position in array is known
        __device__ void addQuadrilateral(const int pos, const int4 q)
        {
            quadrilaterals[pos] = q;
            setColor(pos, INVALID_COLOR);
        }
        /// add a quadrilateral to index buffer, position in array is known
        __device__ void addQuadrilateral(const int pos, const int v0, const int v1, const int v2, const int v3)
        {
            quadrilaterals[pos].x = v0;
            quadrilaterals[pos].y = v1;
            quadrilaterals[pos].z = v2;
            quadrilaterals[pos].w = v3;
            setColor(pos, INVALID_COLOR);
        }
        /// add a quadrilateral to index buffer, position in array is computed via atomicAdd
        __device__ int addColoredQuadrilateral(const int4 q, const int c)
        {
            if (q.x >= 0 && q.y >= 0 && q.z >= 0 && q.w >= 0)
            {
                int pos = atomicAdd(t_size, 1);
                quadrilaterals[pos] = q;
                //setColor(pos, c);
                initAttribute(pos, c); // this value is been set for the very first time!
                return pos;
            }
            else {
                return -1;
            }
        }
        /// add a quadrilateral to index buffer, position in array is computed via atomicAdd
        __device__ int addColoredQuadrilateral(const int v0, const int v1, const int v2, const int v3, const int c)
        {
            if (v0 > -1 && v1 > -1 && v2 > -1 && v3 > -1)
            {
                int pos = atomicAdd(t_size, 1);
                quadrilaterals[pos] = { v0,v1,v2,v3 };
                //setColor(pos, c);
                initAttribute(pos, c); // this value is been set for the very first time!
                return pos;
            }
            else
            {
                return -1;
            }
        }
        /// add a quadrilateral with a face color, position in array is known
        __device__ void addColoredQuadrilateral(const int pos, const int4 q, const int c)
        {
            quadrilaterals[pos] = q;
            setColor(pos, c);
        }
        /// add a quadrilateral with a face color, position in array is known
        __device__ void addColoredQuadrilateral(const int pos, const int v0, const int v1, const int v2, const int v3, const int c)
        {
            quadrilaterals[pos].x = v0;
            quadrilaterals[pos].y = v1;
            quadrilaterals[pos].z = v2;
            quadrilaterals[pos].w = v3;
            setColor(pos,c);
        }
        /// copy data
        __host__ void copy(Quadrilaterals& q)
        {
            nr_q = q.size();
            if (nr_q > a_size) {
                // needs to resize buffers
                cudaMalloc(&quadrilaterals, nr_q * sizeof(int4));
                cudaCheckError();
                quadrilaterals_.reset(quadrilaterals, cudaFree);
                cudaMalloc(&attributes, nr_q * sizeof(uchar));
                cudaCheckError();
                attributes_.reset(attributes, cudaFree);
                a_size = nr_q;
            }
            cudaMemcpy((void*)quadrilaterals, q.quadrilaterals, nr_q * sizeof(int4), cudaMemcpyDeviceToDevice);
            cudaCheckError();
            cudaMemcpy((void*)attributes, q.attributes, nr_q * sizeof(uchar), cudaMemcpyDeviceToDevice);
            cudaCheckError();
            cudaMemcpy((void*)t_size, q.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
            cudaCheckError();
        }
        __device__ void initAttribute(const int pos)
        {
            setColor(pos, INVALID_COLOR);
            unsetNonManifold(pos);
            clearPatternBits(pos);
        }
        __device__ void initAttribute(const int pos, const int color)
        {
            setColor(pos, color);
            unsetNonManifold(pos);
            clearPatternBits(pos);
        }
        /// <summary>
        /// Set color. There are 24 colors stored in the first 5 bits
        /// of the variabla attributes
        /// </summary>
        /// <param name="pos">index of quadrilateral</param>
        /// <param name="c"></param>
        __device__ void setColor(const int pos, const int c)
        {
            //int ct = attributes[pos] & 0xE0;
            //attributes[pos] = ct | c;
            attributes[pos] = (attributes[pos] & 0xE0) | c;
        }
        /// <summary>
        /// Get the color of quadrilateral at given position in array
        /// </summary>
        /// <param name="pos">index in attributes array</param>
        /// <returns></returns>
        __device__ int getColor(const int pos)
        {
            return static_cast<int>(attributes[pos] & 0x1F);
        }
        /// <summary>
        /// Unset color of quadrilateral, set an invalid color, use color 31
        /// for this application there are only 24 colors.
        /// </summary>
        /// <param name="pos">position is array</param>
        __device__ void unsetColor(const int pos)
        {
            attributes[pos] = attributes[pos] | 0x1F; // set an invalid color
        }
        /// <summary>
        /// Set to 1 if quadrilateral has a non-manifold edge
        /// </summary>
        /// <param name="pos"></param>
        __device__ void setNonManifold(const int pos)
        {
            attributes[pos] = (attributes[pos] | 0x20);
        }
        /// <summary>
        /// Check if quadrilateral has at least one edge which is non-manifold
        /// </summary>
        /// <param name="pos">position is array</param>
        /// <returns>1 if true</returns>
        __device__ bool isNonManifold(const int pos)
        {
            return attributes[pos] & 0x20;
        }
        /// <summary>
        /// Unset flag indicating that quadrilateral is non-manifold
        /// </summary>
        /// <param name="pos">position if array</param>
        __device__ void unsetNonManifold(const int pos)
        {
            attributes[pos] = (attributes[pos] & 0xDF);
        }
        /// <summary>
        /// Clear the last two bits reserved to handle mesh simplification
        /// </summary>
        /// <param name="pos">position in attirbutes array to be cleared</param>
        __device__ void clearPatternBits(const int pos)
        {
            attributes[pos] = attributes[pos] & 0x3F;
        }
        /// <summary>
        /// Check if elements has the vertex valence pattern for mesh simplification
        /// Bit 7 is set, if vertex valence pattern is true
        /// </summary>
        /// <param name="pos">position in attributes array</param>
        /// <returns>true, if valence pattern fulfilled</returns>
        __device__ bool isP3X3Y(const int pos)
        {
            return (attributes[pos] & 0x40);
        }
        __device__ bool isP3333(const int pos)
        {
            return (attributes[pos] & 0x40);
        }
        /// <summary>
        /// Set bit indicating that vertex valence pattern for simplification is fulfilled
        /// </summary>
        /// <param name="pos">position in array</param>
        __device__ void setPatternBit(const int pos)
        {
            attributes[pos] = attributes[pos] | 0x40;
        }
        /// <summary>
        /// Unset the vertex valence pattern of the quadrilateral to neutral.
        /// </summary>
        /// <param name="pos">position in attributes array</param>
        __device__ void unsetPatternBit(const int pos)
        {
            attributes[pos] = attributes[pos] & 0xBF;
        }
        /// <summary>
        /// If element can be removed, the last bit in attributes is set
        /// </summary>
        /// <param name="pos">position in attributes array</param>
        __device__ void setRemoveBit(const int pos)
        {
            attributes[pos] = attributes[pos] | 0x80;
        }
        /// <summary>
        /// If the element should not be removed, the last bit in attributes is unset
        /// </summary>
        /// <param name="pos">position in attributes array</param>
        /// <returns>no return value</returns>
        __device__ void unsetRemoveBit(const int pos)
        {
            attributes[pos] = attributes[pos] & 0x7F;
        }
        /// <summary>
        /// Check if last bit in attributes is set, if element can be removed from mesh
        /// </summary>
        /// <param name="pos">position in attributes array</param>
        /// <returns>true if element can be removed from mesh</returns>
        __device__ bool isRemove(const int pos)
        {
            return attributes[pos] & 0x80;
        }
        /// <summary>
        /// Index of vertex v0 in quadrilateral
        /// </summary>
        /// <param name="pos"></param>
        /// <returns>index of v0</returns>
        __device__ int v0(const int pos)
        {
            return quadrilaterals[pos].x;
        }
        /// <summary>
        /// Index of vertex v1 in quadrilateral
        /// </summary>
        /// <param name="pos">position in array</param>
        /// <returns>index of v1</returns>
        __device__ int v1(const int pos)
        {
            return quadrilaterals[pos].y;
        }
        /// <summary>
        /// Index of vertex v2 in quadrilateral
        /// </summary>
        /// <param name="pos">position is array</param>
        /// <returns>index of v2</returns>
        __device__ int v2(const int pos)
        {
            return quadrilaterals[pos].z;
        }
        /// <summary>
        /// Index of vertex v3 in quadrilateral
        /// </summary>
        /// <param name="pos">position in array</param>
        /// <returns>index of v3</returns>
        __device__ int v3(const int pos)
        {
            return quadrilaterals[pos].w;
        }
        /// convenience methods, copy quadrilaterals from device to host.
        __host__ void getQuadrilaterals(std::vector<int4>& q, std::vector<uchar>& a)
        {
            q.resize(nr_q);
            cudaMemcpy(q.data(), quadrilaterals, nr_q * sizeof(int4), cudaMemcpyDeviceToHost);
            cudaCheckError();
            a.resize(nr_q);
            cudaMemcpy(a.data(), attributes, nr_q * sizeof(uchar), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
    };
} // namespace p_mc
