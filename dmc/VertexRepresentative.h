#pragma once

// CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"
#include "MCPolygon.h"

namespace p_mc {
    struct VertexRepresentative {
        // types
        using ulong = unsigned long long;
		using uint = unsigned int;
        enum class ProjectionDirection { W_PROJECTION = 1, V_PROJECTION = 2, U_PROJECTION = 3 };
        // helper class
        MCPolygon mcp;
        // methods
        /// <summary>
        /// Compute the gradient of the trilinear interpolant within a cell at the given local coordinates.
        /// </summary>
        /// <param name="f">scalar function values at the cell vertices</param>
        /// <param name="u">first local coordinate</param>
        /// <param name="v">second local coordinate</param>
        /// <param name="w">third local coordinate</param>
        /// <returns>gradient of trilinear interpolant at point (u,v,w)</returns>
        __device__ float3 gradient(const float f[8], const float3& p)
        {
            float3 grad;
            grad.x = (1 - p.z) * ((1 - p.y) * (f[1] - f[0]) + p.y * (f[3] - f[2])) + p.z * ((1 - p.y) * (f[5] - f[4]) + p.y * (f[7] - f[6]));
            grad.y = (1 - p.z) * ((1 - p.x) * (f[2] - f[0]) + p.x * (f[3] - f[1])) + p.z * ((1 - p.x) * (f[6] - f[4]) + p.x * (f[7] - f[5]));
            grad.z = (1 - p.y) * ((1 - p.x) * (f[4] - f[0]) + p.x * (f[5] - f[1])) + p.y * ((1 - p.x) * (f[6] - f[2]) + p.x * (f[7] - f[3]));
            return grad;
        }
		/// <summary>
		/// Compute the direction of projection with the largest surface area, this is
        /// the direction with the smallest component of the gradient which is normal
        /// to the surface.
		/// </summary>
		/// <param name="ni">gradient of the trilinear interpolant</param>
		/// <returns>direction of projection</returns>
		__device__ auto lProjection(const float3& ni)
		{
			const double dx = fabsf(ni.x);
			const double dy = fabsf(ni.y);
			const double dz = fabsf(ni.z);
			if (dz >= dx && dz >= dy) return ProjectionDirection::W_PROJECTION;
			if (dy >= dx && dy >= dz) return ProjectionDirection::V_PROJECTION;
			if (dx >= dy && dx >= dz) return ProjectionDirection::U_PROJECTION;
		};
		/// <summary>
		/// In this case the direction of projection is determined by the size of
        /// the bounding box.
		/// </summary>
		/// <param name="u">bounding box extension in u-direction</param>
		/// <param name="v">bounding box extension in v-direction</param>
		/// <param name="w">bounding box extension in w-direction</param>
		/// <returns>direction of projection</returns>
		__device__ auto minProjection(const float u, const float v, const float w)
		{
			if (w <= u && w <= v) return ProjectionDirection::W_PROJECTION;
			if (v <= u && v <= w) return ProjectionDirection::V_PROJECTION;
			if (u <= v && u <= w) return ProjectionDirection::U_PROJECTION;

		};
		/// <summary>
		/// Check if point in unit cell is within the bounding box of one of the other
		/// MC polygons intersecting the cell. This test is only necessary, if the MC polygon
		/// being processed has more than three vertices. The check has to be done if
		/// the neighbor bounding box has size three.
		/// </summary>
		/// <param name="u">local coordinate</param>
		/// <param name="v">local coordinate</param>
		/// <param name="w">local coordinate</param>
		/// <param name="cnt_">number of MC polygons intersecting the unit cell</param>
		/// <param name="c_">ulong encoding the MC polygons</param>
		/// <param name="t">Index of MC polygon being processed</param>
		/// <param name="bboxU">bounding box in u-direction</param>
		/// <param name="bboxB">bounding box in u-direction</param>
		/// <param name="bboxW">bounding box in u-direction</param>
		/// <returns>true in there is an intersection</returns>
		__device__ bool isInNeighborBBox(const float u, const float v, const float w, uint cnt_, ulong& c_, const int t,
			const float2 bboxU[4], const float2 bboxV[4], const float2 bboxW[4])
		{
			for (int i = 0; i < cnt_; i++)
			{
				if (t == i) continue;
				const int cnt_size = mcp.get_cnt_size(t, c_); //(c_ >> (4 * (i + 1))) & 0xF;
				//if (cnt_size == 3 && t != i)
				if (cnt_size == 3)
				{
					if (u < bboxU[i].x || bboxU[i].y < u)
						continue;
					if (v < bboxV[i].x || bboxV[i].y < v)
						continue;
					if (w < bboxW[i].x || bboxW[i].y < w)
						continue;
					return true;
				}
			}
			return false;
		}
		/// <summary>
		/// Sample iso-surface and choose best sample by the distance to the initial estimate of the
		/// vertex representative.
		/// </summary>
		/// <param name="r">projection case, extract indices to evaluate function at cell vertices</param>
		/// <param name="i0">iso-value</param>
		/// <param name="f">function values at cell vertices</param>
		/// <param name="cnt_">number of MC polygons</param>
		/// <param name="c_">MC polygons encoded in a unsigned long long</param>
		/// <param name="t">Index of MC polygon being processed by calling function</param>
		/// <param name="p">estimate of vertex representative</param>
		/// <param name="bboxU">bounding box in u-direction</param>
		/// <param name="bboxV">bounding box in v-direction</param>
		/// <param name="bboxW">bounding box in w-direction</param>
		/// <param name="nrSamples">number of samples to be used</param>
		/// <returns>true if a sample on the iso-surface was found</returns>
		__device__ bool projection(const uint r, const float i0, const float f[8],
			uint cnt_, ulong& c_, const int t, float3& p,
			const float2 bboxU[4], const float2 bboxV[4], const float2 bboxW[4],
			const int nrSamples)
		{
			const float du = (bboxU[t].y - bboxU[t].x) / (nrSamples - 1);
			const float dv = (bboxV[t].y - bboxV[t].x) / (nrSamples - 1);
			float minDistance = 100;
			float3 pt = make_float3(0,0, -1);
			const float eps{ 1e-5 };
			float wmin{ bboxW[t].x };
			float wmax{ bboxW[t].y };
			// consider tiny bounding box
			if (fabs(wmax - wmin) < eps) {
				wmin -= eps;
				wmax += eps;
			}

			const int cnt_size = mcp.get_cnt_size(t, c_);
			for (int i = 1; i < (nrSamples - 1); i++) {
				const float u = bboxU[t].x + i * du;
				for (int j = 1; j < (nrSamples - 1); j++) {
					const float v = bboxV[t].x + j * dv;
					const float g1 = (1 - v) * ((1 - u) * f[0] + u * f[(r >> 3) & 0x7]) + v * ((1 - u) * f[(r >> 6) & 0x7] + u * f[(r >> 9) & 0x7]);
					const float g2 = (1 - v) * ((1 - u) * f[(r >> 12) & 0x7] + u * f[(r >> 15) & 0x7]) + v * ((1 - u) * f[(r >> 18) & 0x7] + u * f[(r >> 21) & 0x7]);
					if (g1 == g2) continue;
					const float w = (i0 - g1) / (g2 - g1);
					if (wmin <= w && w <= wmax)
					{
						if (cnt_ == 1 || cnt_size == 3 || !isInNeighborBBox(u, v, w, cnt_, c_, t, bboxU, bboxV, bboxW))
						{
							const float d = (p.x - u) * (p.x - u) + (p.y - v) * (p.y - v) + (p.z - w) * (p.z - w);
							if (minDistance > d) {
								minDistance = d;
								pt.x = u;
								pt.y = v;
								pt.z = w;
							}
						}
					}
				}
			}
			if (pt.z > -1) {
				p.x = pt.x;
				p.y = pt.y;
				p.z = pt.z;
				return true;
			}
			else return false;
		}
        /// <summary>
        /// Compute point on iso-surface by evaluating the parametric representation of the iso-surface within
		/// the unit cell. Compute fist, a projection direction, compute corresponding 2d parametric representation
		/// and finally sample the iso-surface until a satisfactory vertex representative is found.
        /// </summary>
        /// <param name="i0">iso-value</param>
        /// <param name="f">function values at cell vertices</param>
        /// <param name="cnt_">number of MC polygons</param>
        /// <param name="c_">MC polygons encoded in an unsigned long long integer</param>
        /// <param name="t">MC polygon being processed by calling function</param>
        /// <param name="p">first estimate of vertex representative</param>
        /// <param name="bboxU">bounding box in u-direction (local coordinate)</param>
        /// <param name="bboxV">bounding box in u-direction (local coordinate)</param>
        /// <param name="bboxW">bounding box in u-direction (local coordinate)</param>
        /// <param name="count1">counter for checking, if surface projection failed at first try</param>
        /// <param name="count2">counter for checking if surface projection failed at second try, which means no
		///              vertex representative was found</param>
        /// <returns>void</returns>
        __device__ void representative(const float i0, const float f[8], uint cnt_, ulong& c_, int t, float3& p,
            float2 bboxU[4], float2 bboxV[4], float2 bboxW[4], uint& count1, uint& count2)
        {
			const uint w_proj = 16434824;
			const uint v_proj = 16362248;
			const uint u_proj = 16096528;
			const int nrSamples{ 9 };
            // compute case
            ProjectionDirection prj{ ProjectionDirection::W_PROJECTION };
            if (mcp.get_cnt_size(t, c_) >= 6) prj = lProjection(gradient(f, p));
            else prj = minProjection(bboxU[t].y - bboxU[t].x, bboxV[t].y - bboxV[t].x, bboxW[t].y - bboxW[t].x);
			float3 pt;
            switch (prj) {
            case ProjectionDirection::W_PROJECTION:
            {
				pt.x = p.x;
				pt.y = p.y;
				pt.z = p.z;
                if (!projection(w_proj, i0, f, cnt_, c_, t, pt, bboxU, bboxV, bboxW, nrSamples))
				{
					count1++;
					if (!projection(w_proj, i0, f, cnt_, c_, t, pt, bboxU, bboxV, bboxW, nrSamples * 3 + 1))
					{
						count2++;
					}
                }
				p.x = pt.x;
				p.y = pt.y;
				p.z = pt.z;
                break;
            }
            case ProjectionDirection::V_PROJECTION:
            {
				pt.x = p.x;
				pt.y = p.z;
				pt.z = p.y;
                if (!projection(v_proj, i0, f, cnt_, c_, t, pt, bboxU, bboxW, bboxV, nrSamples))
				{
					count1++;
					if (!projection(v_proj, i0, f, cnt_, c_, t, pt, bboxU, bboxW, bboxV, nrSamples * 3 + 1))
					{
						count2++;
					}
                }
				p.x = pt.x;
				p.y = pt.z;
				p.z = pt.y;
                break;
            }
            case ProjectionDirection::U_PROJECTION:
            {
				pt.x = p.y;
				pt.y = p.z;
				pt.z = p.x;
                if (!projection(u_proj, i0, f, cnt_, c_, t, pt, bboxV, bboxW, bboxU, nrSamples))
				{
					count1++;
					if (!projection(u_proj, i0, f, cnt_, c_, t, pt, bboxV, bboxW, bboxU, nrSamples * 3 + 1))
					{
						count2++;
					}
                }
				p.x = pt.z;
				p.y = pt.x;
				p.z = pt.y;
                break;
            }
            }
        }
    };
} // namespace p_mc
