#pragma once

// C++ libs
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
//#include "helper_cuda.h"
#include "helper_cuda.h"
#include "CTimer.h"
#include "Quadrilaterals.h"
#include "Edges.h"
#include "EdgeHashTable.h"
#include "Halfedges.h"
#include "HalfedgeVertices.h"
#include "HalfedgeFaces.h"
#include "HalfedgeHashTable.h"

namespace p_mc {
    /// <summary>
    /// Construct a halfedge data structure from a 
    /// shared vertex mesh
    /// </summary>
    struct HalfedgeMesh {
        using uint = unsigned int;
        /// <summary>
        /// Collect edges from a quadrilateral mesh and store them 
        /// into a hash table for further processing.
        /// </summary>
        /// <param name="q">quadrilaterals</param>
        /// <param name="eht">hash table to store edges</param>
        void edgeHashTable(Quadrilaterals& q, EdgeHashTable& eht);
        /// <summary>
        /// Computed unique edges from a shared vertex quadrilateral mesh
        /// </summary>
        /// <param name="q">quadrilateral mesh</param>
        /// <param name="e">edges</param>
        /// <returns></returns>
        int edges(Quadrilaterals& q, Edges& e);
        /// <summary>
        /// Compute halfedge data structure from a shared vertex
        /// quadrilateral mesh
        /// </summary>
        /// <param name="nr_v">nr. of vertices in the mesh</param>
        /// <param name="q">quadrilaterals</param>
        /// <param name="he">halfedges</param>
        /// <param name="f">halfedge faces</param>
        /// <param name="v">halfedge vertices</param>
        /// <param name="timer">timer to measure performance</param>
        /// <returns>number of halfedges in the mesh</returns>
        int halfedges(const int nr_v, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& f, HalfedgeVertices& v, CTimer& timer);
    };
} // namespace p_mc