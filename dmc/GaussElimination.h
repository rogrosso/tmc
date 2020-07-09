#pragma once

// cuda stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

namespace p_mc {
	/// <summary>
	/// Solve small linear systems using Gauss Elimnation
	/// </summary>
	struct GaussElimination {
		static const int N{ 3 };
		/// <summary>
		/// Convenience class to locally store a matrix. 
		/// Matrices are assumed to have size Nx(N+1) and store
		/// system matrix and rhd of equation
		/// </summary>
		struct Matrix {
			float m[N * (N + 1)]{ 0 };
			__device__ float& operator() (const int i, const int j) { return m[i * (N + 1) + j]; }
			__device__ const float& operator() (const int i, const int j) const { return m[i * (N + 1) + j]; }
		};
		/// <summary>
		/// Convenience class to locally store a vector. The
		/// vector will contain the solution of the linear system.
		/// </summary>
		struct Vector {
			float v[N]{ 0 };
			__device__ float& operator() (const int i) { return v[i]; }
			__device__ const float& operator() (const int i) const { return v[i]; }
		};
		/// <summary>
		/// Main method which performs the Gauss Elimination and sotr
		/// </summary>
		/// <param name="m"></param>
		/// <param name="x"></param>
		/// <returns></returns>
		__device__ bool gaussianElimination(Matrix& m, Vector& x)
		{
			bool singular_matrix = forwardElimination(m);

			/* if matrix is singular */
			if (singular_matrix)
			{
				return false;
			}

			/* get solution to system and print it using
			backward substitution */
			bool r = backSubstitution(m, x);
			return r;
		}

		// function for elementary operation of swapping two rows 
		__device__ void swap_row(Matrix& m, int i, int j)
		{
			for (int k = 0; k <= N; k++)
			{
				float t_ = m(i, k);
				m(i, k) = m(j, k);
				m(j, k) = t_;
			}
		}

		// function to reduce matrix to r.e.f. 
		__device__ bool forwardElimination(Matrix& m)
		{
			for (int k = 0; k < N; k++)
			{
				int i_max = k;
				int v_max = m(i_max, k);

				/* find greater amplitude for pivot if any */
				for (int i = k + 1; i < N; i++)
				{
					if (fabsf(m(i, k)) > v_max)
					{
						v_max = m(i, k);
						i_max = i;
					}
				}

				/* if a prinicipal diagonal element is zero,
				* it denotes that matrix is singular, and
				* will lead to a division-by-zero later. */
				if (m(k, i_max) == 0.0)
				{
					return true; // Matrix is singular 
				}

				/* Swap the greatest value row with current row */
				if (i_max != k)
					swap_row(m, k, i_max);


				for (int i = k + 1; i < N; i++)
				{
					/* factor f to set current row kth element to 0,
					* and subsequently remaining kth column to 0 */
					float f = m(i, k) / m(k, k);

					/* subtract fth multiple of corresponding kth
					row element*/
					for (int j = k + 1; j <= N; j++)
					{
						m(i, j) -= m(k, j) * f;
					}

					/* filling lower triangular matrix with zeros*/
					m(i, k) = 0.0;
				}
			}
			return false;
		}

		// function to calculate the values of the unknowns 
		__device__ bool backSubstitution(Matrix& m, Vector& x)
		{
			/* Start calculating from last equation up to the
			first */
			for (int i = N - 1; i >= 0; i--)
			{
				/* start with the RHS of the equation */
				x(i) = m(i, N);

				/* Initialize j to i+1 since matrix is upper
				triangular*/
				for (int j = i + 1; j < N; j++)
				{
					/* subtract all the lhs values
					* except the coefficient of the variable
					* whose value is being calculated */
					x(i) -= m(i, j) * x(j);
				}

				/* divide the RHS by the coefficient of the
				unknown being calculated */
				x(i) = x(i) / m(i, i);
			}

			return true;
		}
		
	};
} // namespace 