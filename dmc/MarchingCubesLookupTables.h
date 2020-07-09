#pragma once

// C++ libs
#include <memory>
#include <array>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"

namespace p_mc {
    /// <summary>
    /// Collect all DMC kookup tables into a single data structure
    /// for convenience
    /// </summary>
    struct MarchingCubesLookupTables {
        using ushort = unsigned short;
        ushort* e_pattern{ nullptr };
        std::shared_ptr<ushort> e_pattern_{ nullptr };
        char* t_pattern{ nullptr };
        std::shared_ptr<char> t_pattern_{ nullptr };
        /// constructors
        __host__ MarchingCubesLookupTables() {}
        __host__ MarchingCubesLookupTables(const std::array<ushort, 256>& e, const std::array<char, 4096>& t)
        {
            // alloc and init e_pattern
            cudaMalloc(&e_pattern, 256 * sizeof(ushort));
            cudaCheckError();
            cudaMemcpy(e_pattern, &e[0], 256 * sizeof(ushort), cudaMemcpyHostToDevice);
            cudaCheckError();
            // alloc and init t_pattern
            cudaMalloc(&t_pattern, 4096 * sizeof(char));
            cudaCheckError();
            cudaMemcpy(t_pattern, &t[0], 4096 * sizeof(char), cudaMemcpyHostToDevice);
            cudaCheckError();
            e_pattern_ = std::shared_ptr<ushort>(e_pattern, cudaFree);
            t_pattern_ = std::shared_ptr<char>(t_pattern, cudaFree);
        }
        /// desctructor
        __host__ ~MarchingCubesLookupTables()
        {
            e_pattern_.reset();
            t_pattern_.reset();
            e_pattern = nullptr;
            t_pattern = nullptr;
        }
    };
}// namespace p_mc