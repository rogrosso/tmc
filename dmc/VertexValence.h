#pragma once

// C++
#include <memory>
#include <vector>
// CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project
#include "helper_cuda.h"
#include "ValenceHashTable.h"
#include "Triangles.h"
#include "Quadrilaterals.h"

namespace p_mc {
    struct VertexValence {
        // types
        using uint = unsigned int;
        // methods
        void vertexValence(const int nr_v, Quadrilaterals& q_, std::vector<uint>& valence);
        void vertexValence(const int nr_v, Triangles& t_, std::vector<uint>& valence);
    };
} // namespace p_mc
