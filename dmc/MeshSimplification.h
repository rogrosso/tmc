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
    /// Simplify mesh by removing elements with vertex valence pattern
    /// 3X3Y, X,Y >= 5, and 3333. 
    /// </summary>
    struct MeshSimplification {
        /// <summary>
        /// Simplify elements with vertex valence pattern 3X3Y. Use face coloring
        /// to simplify mesh. 
        /// </summary>
        /// <param name="v">vertices</param>
        /// <param name="q">quadrilaterals</param>
        /// <param name="he">halfedges</param>
        /// <param name="hef">halfedge faces</param>
        /// <param name="hev">halfedge vertices</param>
        /// <param name="timer">timer to measure performance</param>
        void pattern3X3Y(Vertices& v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& hef, HalfedgeVertices& hev, CTimer& timer);
        /// <summary>
        /// Simplify elements with vertex valence pattern 3X3Y after the color based method has finished.
        /// The face coloring is at this stage not valid. The method finds isolated elements, i.e. with
        /// no neighbors with the same vertex valence pattern and simplify these elements.
        /// </summary>
        /// <param name="v">vertices</param>
        /// <param name="q">quadrilaterals</param>
        /// <param name="he">halfedges</param>
        /// <param name="hef">halfedge faces</param>
        /// <param name="hev">halfedge vertices</param>
        /// <param name="timer">timer to measure performance</param>
        void pattern3X3YOld(Vertices& v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& hef, HalfedgeVertices& hev, CTimer& timer);
        /// <summary>
        /// Simplify elements with vertex valence pattern 3333
        /// </summary>
        /// <param name="v">vertices</param>
        /// <param name="q">quadrilaterals</param>
        /// <param name="he">halfedges</param>
        /// <param name="hef">halfedge faces</param>
        /// <param name="hev">halfedge vertices</param>
        /// <param name="timer">timer to measure performance</param>
        void pattern3333(Vertices& v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& hef, HalfedgeVertices& hev, CTimer& timer);
    };
} // namespace p_mc