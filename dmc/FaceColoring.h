#pragma once

// C++ libs
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"
#include "CTimer.h"
#include "Vertices.h"
#include "Quadrilaterals.h"
#include "Halfedges.h"
#include "HalfedgeFaces.h"
#include "HalfedgeVertices.h"
#include "HalfedgeMesh.h"
#include "Edges.h"
#include "EdgeHashTable.h"
#include "VertexMap.h"
#include "QuadrilateralMap.h"

namespace p_mc {
    /// <summary>
    /// Faces in the DMC mesh will be colored using five colors. 
    /// The initial coloring is inherited from an edge coloring
    /// of the input grid consisting of 24 colors, from 0 t0 23.
    /// This class will eliminate 19 colors, the colors 5 to 23.
    /// C is the index of the first color to be eliminated, cSz
    /// is the total number of colors to be eliminated by this class.
    /// </summary>
    constexpr int C{ 5 };
    constexpr int cSz{ 19 };
    /// <summary>
    /// Computes the face coloring of the DMC mesh consisting of
    /// a total of 5 colors.
    /// </summary>
    struct FaceColoring {
        /// <summary>
        /// Helper class to count how many
        /// faces have a given color.
        /// </summary>
        struct ColorCount {
            int nr_f{ 0 };
            /// <summary>
            /// List of indices of a given color
            /// </summary>
            int* colors[cSz];
            std::vector<std::shared_ptr<int>> colors_;
            /// <summary>
            /// atomic counter, given size of index array
            /// </summary>
            int* t_size[cSz];
            std::vector<std::shared_ptr<int>> t_size_;
            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="nr_q">nr. of quadrilaterals in the DMC mesh</param>
            /// <returns></returns>
            __host__ ColorCount(const int nr_q)
            {
                nr_f = 3 * nr_q / cSz;
                colors_.resize(cSz);
                t_size_.resize(cSz);
                for (int i = 0; i < cSz; i++)
                {
                    cudaMalloc(&colors[i], nr_f * sizeof(int));
                    p_mc::cudaError(__FILE__, __LINE__);
                    colors_[i].reset(colors[i], cudaFree);
                    //colors_[i] = std::make_shared<int>(colors[i], cudaFree);
                    cudaMalloc(&t_size[i], sizeof(int));
                    p_mc::cudaError(__FILE__, __LINE__);
                    cudaMemset(t_size[i], 0, sizeof(int));
                    t_size_[i].reset(t_size[i], cudaFree);
                }
            }
            /// <summary>
            /// Destructor
            /// </summary>
            /// <returns></returns>
            __host__ ~ColorCount()
            {
                for (auto e : colors_)
                {
                    e.reset();
                }
                for (auto e : t_size_)
                {
                    e.reset();
                }
            }
            /// <summary>
            /// Set size of buffers
            /// </summary>
            /// <param name="nr_q">nr. of quadrilaterals in the DMC mesh</param>
            /// <returns></returns>
            __host__ void resize(const int nr_q)
            {
                const int sz = 3 * nr_q / cSz;
                colors_.resize(cSz);
                t_size_.resize(cSz);
                for (int i = 0; i < cSz; i++)
                {
                    cudaMalloc(&colors[i], sz * sizeof(int));
                    p_mc::cudaError(__FILE__, __LINE__);
                    colors_[i].reset(colors[i], cudaFree);
                    //colors_[i] = std::make_shared<int>(colors[i], cudaFree);
                    cudaMalloc(&t_size[i], sizeof(int));
                    p_mc::cudaError(__FILE__, __LINE__);
                    cudaMemset(t_size[i], 0, sizeof(int));
                    t_size_[i].reset(t_size[i], cudaFree);
                }
            }
            /// <summary>
            /// Buffer size
            /// </summary>
            /// <returns></returns>
            __host__ __device__ int bufferSize() { return nr_f; }
            /// <summary>
            /// Add face index in array of all faces of a given color. It is
            /// used to count the number of elements with a certain color.
            /// </summary>
            /// <param name="c">color</param>
            /// <param name="faceId">face index</param>
            /// <returns></returns>
            __device__ void addColor(const int c, const int faceId)
            {
                const int addr = atomicAdd(t_size[c - C], 1);
                colors[c - C][addr] = faceId;
            }
            /// <summary>
            /// Get face id for input color, which is saved at position i in array.
            /// </summary>
            /// <param name="c">Color</param>
            /// <param name="i">position in array for that color</param>
            /// <returns></returns>
            __device__ int getFace(const int c, const int i) { return colors[c - C][i]; }
            __host__ int size(const int c)
            {
                int nr{ 0 };
                cudaMemcpy(&nr, t_size[c - C], sizeof(int), cudaMemcpyDeviceToHost);
                return nr;
            }

        };

        /// <summary>
        /// Classify faces accroding to color, i.e. it counts how many faces
        /// has a certain color. It requires halfedge data structure to find 
        /// neighbor faces and check for equal colors.
        /// </summary>
        /// <param name="q">quadrilateral faces</param>
        /// <param name="he">halfedge edges</param>
        /// <param name="hef">halfedge faces</param>
        /// <returns></returns>
        __host__ void colorFaces(Quadrilaterals& q, Halfedges he, HalfedgeFaces& hef, CTimer& timer);
    };
} // namespace