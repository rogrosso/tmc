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
    /// Represent a face in a halfedge data structure
    /// </summary>
    struct HalfedgeFaces {
        using uchar = unsigned char;
        /// total number of faces
        int nr_f{ 0 };
        /// index buffer
        int* he_e{ nullptr };
        std::shared_ptr<int> he_e_{ nullptr };
        /// <summary>
        /// Store face color, non-manifold case and tags for mesh simplification
        /// There are at most 24 colors, thus information is stored as follows
        ///   - first five bits for color
        ///   - sixth bit for non-manifold case
        ///   - last two bits for vertex valence pattern of face
        /// </summary>
        uchar* attributes{ nullptr };
        std::shared_ptr<uchar> attributes_;
        /// constructors
        __host__ HalfedgeFaces() {}
        __host__ HalfedgeFaces(const int sz) : nr_f{ sz }
        {
            cudaMalloc(&he_e, sz * sizeof(int));
            cudaCheckError();
            he_e_ = std::shared_ptr<int>(he_e, cudaFree);
            cudaMalloc(&attributes, sz * sizeof(uchar));
            cudaCheckError();
            attributes_.reset(attributes, cudaFree);
        }
        /// destructor
        __host__ ~HalfedgeFaces()
        {
            nr_f = 0;
            he_e_.reset();
            attributes_.reset();
            he_e = nullptr;
            attributes = nullptr;
        }
        /// size of buffer
        __host__ int capacity() { return nr_f; }
        /// total number of faces
        __host__ __device__ int size() { return nr_f; }
        /// change size of buffer
        __host__ void resize(const int sz)
        {
            if (sz != nr_f)
            {
                nr_f = sz;
                cudaMalloc(&he_e, nr_f * sizeof(int));
                cudaCheckError();
                he_e_.reset(he_e, cudaFree);
                cudaMalloc(&attributes, nr_f * sizeof(uchar));
                cudaCheckError();
                attributes_.reset(attributes, cudaFree);
            }
        }
        /// add vertex
        __device__ void addFace(const int pos, const int f)
        {
            he_e[pos] = f;
            setColor(pos, INVALID_COLOR);
            unsetNonManifold(pos);
            clearPatternBits(pos);
        }
        /// <summary>
        /// Add a face with color, this is probably the common method
        /// </summary>
        /// <param name="pos"></param>
        /// <param name="f"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        __device__ void addFace(const int pos, const int f, const int c)
        {
            he_e[pos] = f;
            setColor(pos, c);
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
        /// read halfedge data structure out of device memory
        __host__ void getHalfedgeFaces(std::vector<int>& f)
        {
            f.resize(nr_f);
            cudaMemcpy(f.data(), he_e, nr_f * sizeof(int), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
        /// read halfedge data structure out of device memory
        __host__ void getHalfedgeFaces(std::vector<int>& f,std::vector<uchar>& a)
        {
            f.resize(nr_f);
            cudaMemcpy(f.data(), he_e, nr_f * sizeof(int), cudaMemcpyDeviceToHost);
            cudaCheckError();
            a.resize(nr_f);
            cudaMemcpy(a.data(), attributes, nr_f * sizeof(uchar), cudaMemcpyDeviceToHost);
            cudaCheckError();
        }
    };
} // namespace p_mc
