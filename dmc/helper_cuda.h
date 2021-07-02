#pragma once

// C++ libs
#include <iostream>

// cuda stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace p_mc {
    /// Note: inline constexpr is C++17. The NVIDIA compiler can't handle this.
    /// defines for the cell intersection class and methods
    constexpr int BIT_1{ 0x1 };
    constexpr int BIT_2{ 0x2 };
    constexpr int BIT_3{ 0x4 };
    constexpr int BIT_4{ 0x8 };
    constexpr int BIT_5{ 0x10 };
    constexpr int BIT_6{ 0x20 };
    constexpr int BIT_7{ 0x40 };
    constexpr int BIT_8{ 0x80 };
    constexpr int BIT_16{ 0x8000 };
    /// mark ambieguous cases, requiring asymtotic decider
    constexpr int MC_AMBIGUOUS{ 105 };

    /// <summary>
    /// Blocksize for CUDA kernels
    /// </summary>
    constexpr int MC_BLOCKSIZE{ 128 };
    // hash tables and element coloring
    constexpr int INVALID_INDEX = -1;
    /// <summary>
    /// There are only 24 colors, set color to an invalid value
    /// </summary>
    constexpr int INVALID_COLOR = 0x1F;
    /// <summary>
    /// Empty bucket in a hash table, where the keys are unsigned long long
    /// i.e. 64 bit unsigned integers
    /// </summary>
    constexpr unsigned long long EMPTY_BUCKET_64 = 0ull;
    /// <summary>
    /// Empty bucket in a hash table with keys which are signed 32bit integers
    /// </summary>
    constexpr int EMPTY_BUCKET_32 = -1;
    /// <summary>
    /// Quadrilateral with no particular valence pattern
    /// </summary>
    constexpr int P_NEUTRAL = 0;
    /// <summary>
    /// Elements with valence pattern 3X3Y, X >= 5 and Y >= 5
    /// Vertices incident to an element with valence pattern 3X3Y
    /// </summary>
    constexpr int P_3X3Y = 1;
    /// <summary>
    /// Elements with valence pattern 3333
    /// Vertices incident to an element with valence pattern 3333
    /// </summary>
    constexpr int P_3333 = 1;
    /// <summary>
    /// Neighbor of an element with valence pattern 3333
    /// The element will be removed to simplify the mesh.
    /// In case of neighbors all with valence pattern 3333, elements cannot
    /// be removed, e.g. a hexahedron
    /// </summary>
    constexpr int P_N_3333 = 2;
    /// <summary>
    /// Mark an element or vertex for removal
    /// </summary>
    constexpr int P_REMOVE = 3;
    /// <summary>
    /// Non manifold element, i.e. it shares an edge with more than two quadrilaterals
    /// These kind of elements cannot be simplified.
    /// </summary>
    constexpr int Q_NONMANIFOLD = 0xF;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // error handling
    //#define cudaCheckError() { \
    //	cudaError_t e=cudaGetLastError(); \
    //	if(e!=cudaSuccess) { \
    //		using d_out = utilities::DebugString; \
    //		std::ostringstream buf; \
    //		buf << "Cuda failure in file " << __FILE__ << ", line: " << __LINE__ << ", error: " << cudaGetErrorString(e) << "\n";  \
    //        d_out::print(buf.str()); \
    //		exit(0); \
    //	} \
    //}

    inline void cudaError_(const char* fileName, const size_t lineNumber)
    {
        cudaError_t e = cudaGetLastError();
        if (e != cudaSuccess)
        {
            std::cout << "Cuda failure in file " << fileName << ", line: " << lineNumber << ", error: " << cudaGetErrorString(e) << std::endl;
            exit(0);
        }
    }

#define cudaCheckError() cudaError_(__FILE__, __LINE__)


} // namespace
