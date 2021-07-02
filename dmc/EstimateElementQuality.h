#pragma once

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "Vertices.h"
#include "Triangles.h"
#include "Quadrilaterals.h"
#include "ElementQuality.h"

namespace p_mc {
    struct EstimateElementQuality {
        using uint = unsigned int;
        __host__ void q(Vertices v_, Quadrilaterals q_, ElementQuality eQ);
        __host__ void q(Vertices v_, Triangles t_, ElementQuality eQ);
    };
} // namespace p_mc
