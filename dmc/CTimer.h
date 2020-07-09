#pragma once

// C++ libs
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"

namespace p_mc {
	/// <summary>
	/// Use CUDA routines to measure computation time.
	/// </summary>
	struct CTimer {
		float e_milliseconds;
		cudaEvent_t c_start;
		cudaEvent_t c_stop;
		/// <summary>
		/// Constructor
		/// </summary>
		/// <returns></returns>
		CTimer() {
			cudaEventCreate(&c_start);
			cudaEventCreate(&c_stop);
		}
		/// <summary>
		/// Start timer
		/// </summary>
		/// <returns></returns>
		void __host__ start() {
			cudaEventRecord(c_start);
		}
		/// <summary>
		/// Stop timer
		/// </summary>
		/// <returns></returns>
		void __host__ stop() {
			cudaEventRecord(c_stop);
			cudaEventSynchronize(c_stop);
			cudaEventElapsedTime(&e_milliseconds, c_start, c_stop);
		}
		/// <summary>
		/// Print elapsed time between the call of start() and stop()
		/// </summary>
		/// <returns></returns>
		void __host__ print() {
			std::ostringstream buf;
			buf << std::setprecision(7) << " ... time in ms: " << e_milliseconds << std::endl;
			std::cout << buf.str() << std::endl;
		}
		/// <summary>
		/// Print elapsed time between the call of start() and stop(),
		/// print a user message
		/// </summary>
		/// <param name="m">user message</param>
		/// <returns></returns>
		void __host__ print(std::string& m) {
			std::ostringstream buf;
			buf << std::setprecision(7) << " ... " << m << " time in ms: " << e_milliseconds << std::endl;
			std::cout << buf.str() << std::endl;
		}
		/// <summary>
		/// Get elapsed time as string
		/// </summary>
		/// <param name="m"></param>
		/// <returns></returns>
		std::string getTime(const std::string& m) {
			std::ostringstream buf;
			buf << std::setprecision(7) << " ... " << m << " time in ms: " << e_milliseconds << std::endl;
			return buf.str();
		}

	};
} // namespace p_mc