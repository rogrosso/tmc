#pragma once

// Libs
#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include <fstream>
#include <list>
#include <vector>
#include <array>
#include <stack>
#include <map>
#include <utility>
#include <bitset>
#include <algorithm>
#include <random>
#include <chrono>

// project
#include "Vector.h"
#include "UniformGrid.h"

namespace cpp_mc {

	//
	//  DualMarchingCubes.h
	//
	//  Created by Roberto Grosso on 05.06.16.
	//  Copyright � 2016 Roberto Grosso. All rights reserved.
	//

	// constants
	constexpr int MC_AMBIGUOUS{ 105 };
	constexpr int INVALID_INDEX{ -1 };
	constexpr int INVALID_COLOR{ -1 };
	constexpr int MANIFOLD{ 0 };
	constexpr int NON_MANIFOLD{ 1 };
	constexpr int BIT_1{ 0x1 };
	constexpr int BIT_2{ 0x2 };
	constexpr int BIT_3{ 0x4 };
	constexpr int BIT_4{ 0x8 };
	constexpr int BIT_5{ 0x10 };
	constexpr int BIT_6{ 0x20 };
	constexpr int BIT_7{ 0x40 };
	constexpr int BIT_8{ 0x80 };
	constexpr int BIT_16{ 0x8000 };

		/** Computes a topologically correct and manifolg isosurface
		 *  from an uniform grid.
		 *  # Volume data
		 *  The implementation assumes the volume data is sotred in a uniform grid. The file containing the data
		 *  should be a binary file with the following format:
		 *  - three unsigned short indicating the dimensions in x, y and z
		 *  - three floats with the cell sizes in x, y and z
		 *  - the valume data as unsigned shorts.
		 *  # Element topology
		 *  The numbering of vertices is given by vertex coordinates of a unit reference hexahedron. Edges
		 *  are given by their end vertices. There are alos lists with the face-vertex and face-edge correspondences.
		 *  ## Vertices
		 *  - v0 = (0,0,0)
		 *  - v1 = (1,0,0)
		 *  - v2 = (1,0,0)
		 *  - v3 = (1,0,0)
		 *  - v4 = (1,0,0)
		 *  - v5 = (1,0,0)
		 *  - v6 = (1,0,0)
		 *  - v7 = (1,0,0)
		 *
		 *  ## Edges
		 *  - e0 = {0,1}
		 *  - e1 = {1,3}
		 *  - e2 = {2,3}
		 *  - e3 = {0,2}
		 *  - e4 = {4,5}
		 *  - e5 = {5,7}
		 *  - e6 = {6,7}
		 *  - e7 = {4,6}
		 *  - e8 = {0,4}
		 *  - e9 = {1,5}
		 *  - e10 = {3,7}
		 *  - e11 = {2,6}
		 *
		 *  ## Face-Edge correspondence
		 *  The order of the entries in the list corresponds to the orientation of the face
		 *  - f0 = {0,1,2,3}
		 *  - f1 = {4,5,6,7}
		 *  - f2 = {0,9,4,8}
		 *  - f3 = {2,10,6,11}
		 *  - f4 = {3,11,7,8}
		 *  - f5 = {1,10,5,9}
		 *
		 *  ## Face-Vertex correspondence
		 *  - f0 = {0,1,2,3}
		 *  - f1 = {4,5,6,7}
		 *  - f2 = {0,1,4,5}
		 *  - f3 = {2,3,6,7}
		 *  - f4 = {0,2,4,6}
		 *  - f5 = {1,3,5,7}
		 *
		 *  ## Orientation of contours
		 *  In order to construct contour which are oriented such that positive vertices are outside
		 *  we describes faces as been seeing from outside the hexahedron. In this case, the numbering
		 *  used to collect vertices for a face is as follows:
		 *  ### Oriented numbering of Face-Vertex correspondence
		 *  This correspondence is used to compute the intersection of the isosurface with the face
		 *  of the hexahedron and orient the segments such that positive vertices are outside the
		 *  contour
		 *  - f0 = {0,2,1,3}
		 *  - f1 = {5,7,4,6}
		 *  - f2 = {4,0,5,1}
		 *  - f3 = {2,6,3,7}
		 *  - f4 = {4,6,0,2}
		 *  - f5 = {1,3,5,7}
		 *
		 *   ### Oriented numbering of Face-Edge correspondence
		 *  - f0 = {3,2,1,0}
		 *  - f1 = {5,6,7,4}
		 *  - f2 = {8,0,0,4}
		 *  - f3 = {11,6,10,2}
		 *  - f4 = {7,11,3,8}
		 *  - f5 = {1,10,5,9}
		 *
		 *  # Remark
		 *  The function `slice()` implements the algorithm to compute the intersection of the iso-surface
         *  with the cell.
		 */
		class DualMarchingCubes {
		public:
			using uchar = unsigned char;
			using ushort = unsigned short;
			using uint = unsigned int;
			using ulong = unsigned long long;
			using Scalar = double;

			using Vertex = cpp_mc::Vector;
			using Point  = cpp_mc::Vector;
			using Normal = cpp_mc::Vector;
			using UGrid  = cpp_mc::UniformGrid;
			using Index  = UGrid::Index;


		public:
			/**
			 * Computes the duals marching cubes isosurface.
			 * @param[in] i0 iso-value
			 * @param[in] ugrid input uniform grid
			 * @param[out] vertices
             * @param[out] normals
             * @param[out] tris index list with m_triangles
             * @param[out] quads index list with quadrilaterals
			 */
			void dual_mc(const double i0, UGrid& ugrid,
				std::vector<Vertex>& vertices, std::vector<Normal>& normals,
				std::vector<int>& triangles, std::vector<int>& quads);

		private:
			/// <summary>
			/// Compute the intersection of the iso-surface with a cell
			/// </summary>
			/// <param name="i0">iso-value</param>
			/// <param name="i_case">MC case</param>
			/// <param name="i">i-index of cell</param>
			/// <param name="j">j-index of cell</param>
			/// <param name="k">k-index of cell</param>
			/// <param name="f">the eight values of the scalar at the cell vertices</param>
			/// <param name="ugrid">Uniform grid</param>
			/// <param name="v">vertices</param>
			/// <param name="n">normals</param>
			/// <param name="map">a hash map to reconstruct quadrilaterals</param>
			void slice(const double i0, const uint i_case, const int i, const int j, const int k, double f[8],
				UGrid& ugrid, std::vector<Vertex>& v, std::vector<Normal>& n, std::map<int, std::array<int, 5>>& map);
			/// <summary>
			/// Compute the MC polygon for the ambiguos cases
			/// </summary>
			/// <param name="i0">iso-values</param>
			/// <param name="f">scalar values at cell vertices</param>
			/// <param name="c_">mc polygon</param>
			/// <returns></returns>
			uint mc_polygon(const double i0, const double f[8], ulong& c_);
			/// <summary>
			/// Compute the MC polygon for the unambiguos case out of a modified MC lookup table.
			/// There might give up to four contours encoded in the unsigned long long c_
			/// </summary>
			/// <param name="i_case">MC case</param>
			/// <param name="c_">mc polygons</param>
			/// <returns></returns>
			uint mc_polygon(const int i_case, ulong& c_);
			/// <summary>
			/// Set size of MC polygon encoded in an unsigned long long c_. There might give
			/// up to four contours encoded in c_
			/// </summary>
			/// <param name="cnt">mc polygon id</param>
			/// <param name="size">new size</param>
			/// <param name="c_">MC polygons</param>
			void set_cnt_size(const int cnt, const int size, unsigned long long& c_) {
				// unset contour size
				c_ &= ~(0xF << 4 * cnt);
				c_ |= (size << 4 * cnt);
			}
			/// <summary>
			/// Get size of contour encoded in c_ (unsigned long long)
			/// </summary>
			/// <param name="cnt">mc polygon id</param>
			/// <param name="c_">MC polygons</param>
			/// <returns></returns>
			int get_cnt_size(const int cnt, unsigned long long& c_) {
				return static_cast<int>((c_ & (0xF << 4 * cnt)) >> 4 * cnt);
			}
			/// <summary>
			/// Set edge building a MC polygon
			/// </summary>
			/// <param name="cnt">mc polygon id</param>
			/// <param name="pos">position</param>
			/// <param name="val">edge id</param>
			/// <param name="c_">MC polygons</param>
			void set_c(const int cnt, const int pos, const int val, unsigned long long& c_) {
				const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
				const uint c_sz = c_ & mask[cnt];
				const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
				c_ &= ~(((unsigned long long)0xF) << e);
				c_ |= (((unsigned long long)val) << e);
			}
			/// <summary>
			/// Read edge from MC polygon
			/// </summary>
			/// <param name="cnt">polygon id</param>
			/// <param name="pos">edge position in polygon</param>
			/// <param name="c_">MC polygons</param>
			/// <returns></returns>
			int get_c(const int cnt, const int pos, unsigned long long c_) {
				const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
				const uint c_sz = (uint)(c_ & mask[cnt]);
				const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
				return static_cast<int>((c_ >> e) & 0xF);
			}
            /// <summary>
            /// Computes the intersection of the mc polygon with the edge and add the position
            /// to the vertex representative for this polygon. This way, the vertex representative
            /// can be computed as the mean value of the intersections of the mc polygon with the cell.
            /// The vertex representative is given in local coordinates.
            /// </summary>
            /// <param name="l_edge">vertices defining the edge</param>
            /// <param name="e">edge id</param>
            /// <param name="f">function values at cell vertices</param>
            /// <param name="i0">iso-value</param>
            /// <param name="u">vertex representative</param>
            /// <returns></returns>
			void getLocalCoordinates(const unsigned char l_edge, const int e, const double f[8], const double i0, Vertex& u)
			{
				const int v0 = (l_edge & 0xF);
				const int v1 = (l_edge >> 4) & 0xF;
				const float l = (i0 - f[v0]) / (f[v1] - f[v0]);
				unsigned long long int e_{ 75059404284193ULL };
				unsigned long long c_{ 38552806359568ULL };
				float val = (e_ & 1ULL << (4 * e)) >> (4 * e);
				float con = (c_ & 1ULL << (4 * e)) >> (4 * e);
				u[0] += l * val + con;
				val = (e_ & 1ULL << (4 * e + 1)) >> (4 * e + 1);
				con = (c_ & 1ULL << (4 * e + 1)) >> (4 * e + 1);
				u[1] += l * val + con;
				val = (e_ & 1ULL << (4 * e + 2)) >> (4 * e + 2);
				con = (c_ & 1ULL << (4 * e + 2)) >> (4 * e + 2);
				u[2] += l * val + con;
			}
            /// <summary>
            /// Move the vertex representative for a given mc polygon towards
            /// the iso-surface within the unit cell.
            /// </summary>
            /// <param name="f">scalar calues at cell vertices</param>
            /// <param name="i0">iso-value</param>
            /// <param name="p">vertex representative, which will be moved toward iso-surface</param>
            /// <returns></returns>
			void movePointToSurface(const double f[8], const double i0, Point& p)
			{
				const int max_iter = 20;
				double u1 = p[0], v1 = p[1], w1 = p[2];
				double val1 = trilinear(f, u1, v1, w1);
				double val2{ 0 };
				for (int iter = 0; iter < max_iter; iter++)
				{
					Vector grad = gradient(f, u1, v1, w1);
					if (val1 > i0)
					{
						grad.flip();
					}
					// normalize
					normalize(grad);
					const double step = 0.05f;
					double u2 = u1 + step * grad[0];
					double v2 = v1 + step * grad[1];
					double w2 = w1 + step * grad[2];

					// check if we are within the cell
					bool c = clamp(u2, v2, w2);
					//const float val2 = trilinear(f, u2, v2, w2);
					val2 = trilinear(f, u2, v2, w2);
					if ((val1 <= i0 && i0 < val2) || (val2 < i0 && i0 <= val1))
					{
						const float e = (i0 - val1) / (val2 - val1);
						p[0] = u1 + e * (u2 - u1);
						p[1] = v1 + e * (v2 - v1);
						p[2] = w1 + e * (w2 - w1);
						break;
					}
					if (c) break;
					// update
					u1 = u2;
					v1 = v2;
					w1 = w2;
					val1 = val2;
				}
			}
			/// <summary>Clamp local coordinates to unit cube</summary>
			/// <param name="u"></param>
			/// <param name="v"></param>
			/// <param name="w"></param>
			/// <returns>true is coordinated clamped</returns>
			bool clamp(double& u, double& v, double& w)
			{
				bool flag{ false };
				(u <= 0 ? (u = 0, flag = true) : u <= 1 ? u : (u = 1, flag = true));
				(v <= 0 ? (v = 0, flag = true) : v <= 1 ? v : (v = 1, flag = true));
				(w <= 0 ? (w = 0, flag = true) : w <= 1 ? w : (w = 1, flag = true));
				return flag;
			}
            /// <summary>Trilinear interpolation within unit cell</summary>
            /// <param name="f">scalar values at cell vertices</param>
			/// <param name="u">u-local coordinate</param>
			/// <param name="v">v-local coordinate</param>
            /// <param name="w">w-local coordinate</param>
			/// <returns>the interpolated value at the local coordinates (u,v,w)</returns>
			double trilinear(const double f[8], const double u, const double v, const double w)
			{
				return (1 - w) * ((1 - v) * ((1 - u) * f[0] + u * f[1]) + v * ((1 - u) * f[2] + u * f[3]))
					+ w * ((1 - v) * ((1 - u) * f[4] + u * f[5]) + v * ((1 - u) * f[6] + u * f[7]));
			}
            /// <summary>
            /// computes the gradient of the scalar field in terms of the
            /// local coordinates as it is defined whithin the unit cells
            /// </summary>
            /// <param name="f">scalar values at cell vertices</param>
			/// <param name="u">u-local coordinate</param>
			/// <param name="v">v-local coordinate</param>
            /// <param name="w">w-local coordinate</param>
			/// <returns>the value of the gradient at the local coordinates (u,v,w)</returns>
			Vector gradient(const double f[8], const double u, const double v, const double w)
			{
				Vector g;
				g[0] = (1 - w) * ((1 - v) * (f[1] - f[0]) + v * (f[3] - f[2])) + w * ((1 - v) * (f[5] - f[4]) + v * (f[7] - f[6]));
				g[1] = (1 - w) * ((1 - u) * (f[2] - f[0]) + u * (f[3] - f[1])) + w * ((1 - u) * (f[6] - f[4]) + u * (f[7] - f[5]));
				g[2] = (1 - v) * ((1 - u) * (f[4] - f[0]) + u * (f[5] - f[1])) + v * ((1 - u) * (f[6] - f[2]) + u * (f[7] - f[3]));
				return g;
			}
            /// <summary>
            /// Add vertex to hash map, required to reconstruct the quadrilaterals
            /// in the following processing step.
            /// </summary>
            /// <param name="key">key, unique index of edge in uniform grid</param>
			/// <param name="pos">position of vertex index within the quadrilateral</param>
            /// <param name="v_addr">vertex index in vertex buffer</param>
			/// <param name="m_">the hash map</param>
			/// <returns>true is success</returns>
			bool addVertex(const int key, const int pos, const int v_addr, std::map<int, std::array<int, 5>>& m_)
			{
				auto e = m_.find(key);
				if (e != m_.end())
				{
					// quadrilateral already exists
					e->second[pos] = v_addr;
				}
				else
				{
					// create quadrilateral
					std::array<int, 5> q{ INVALID_INDEX,INVALID_INDEX,INVALID_INDEX,INVALID_INDEX, INVALID_COLOR };
					q[pos] = v_addr;
					m_[key] = q;
				}
				return true;
			}
            /// <summary>
            /// Add vertex to hash map, required to reconstruct the quadrilaterals
            /// in the following processing step. This function adds a vertex index and colors
            /// for the given quadrilateral
            /// </summary>
            /// <param name="key">key, unique index of edge in uniform grid</param>
			/// <param name="pos">position of vertex index within the quadrilateral</param>
            /// <param name="v_addr">vertex index in vertex buffer</param>
            /// <param name="color">face color</param>
			/// <param name="m_">the hash map</param>
			/// <returns>true is success</returns>
			bool addVertex(const int key, const int pos, const int v_addr, const int color, std::map<int, std::array<int, 5>>& m_)
			{
				auto e = m_.find(key);
				if (e != m_.end())
				{
					// quadrilateral already exists
					e->second[pos] = v_addr;
					e->second[4] = color;
				}
				else
				{
					// create quadrilateral
					std::array<int, 5> q{ INVALID_INDEX,INVALID_INDEX,INVALID_INDEX,INVALID_INDEX, INVALID_COLOR };
					q[pos] = v_addr;
					q[4] = color;
					m_[key] = q;
				}
				return true;
			}
            /// <summary>
            /// collect quadrilaterals from hash map.
            /// </summary>
            /// <param name="quads">index array containing the quadrilaterals</param>
			/// <param name="colors">index array containing the colors of the faces</param>
            /// <param name="m_">hash map</param>
			/// <returns></returns>
			void collectQuadrilaterals(std::vector<int>& quads, std::vector<int>& colors, std::map<int, std::array<int, 5>> m_)
			{
				for (auto [k, q] : m_)
				{
					if (q[0] != INVALID_INDEX && q[1] != INVALID_INDEX && q[2] != INVALID_INDEX && q[3] != INVALID_INDEX)
					{
						quads.push_back(q[0]);
						quads.push_back(q[1]);
						quads.push_back(q[2]);
						quads.push_back(q[3]);
						colors.push_back(q[4]);
					}
				}
			}
            /// <summary>
            /// compute triangles from the list of quadrilaterals. Use the min-angle cirterion
            /// to compute the best quality triangles.
            /// </summary>
            /// <param name="tris">index array containing the triangles</param>
			/// <param name="colors">index array containing quadrilaterals</param>
            /// <param name="v">vertices required to compute the min-angle criterion</param>
			/// <returns></returns>
			void collectTriangles(std::vector<int>& tris, std::vector<int>& quads, std::vector<Vertex>& v)
			{
				const int nr_q = static_cast<int>(quads.size()) / 4;
				for (int i = 0; i < nr_q; i++)
				{
					const int v0 = quads[4 * i];
					const int v1 = quads[4 * i + 1];
					const int v2 = quads[4 * i + 2];
					const int v3 = quads[4 * i + 3];
					/*tris.push_back(v0);
					tris.push_back(v1);
					tris.push_back(v2);
					// second triangle
					tris.push_back(v0);
					tris.push_back(v2);
					tris.push_back(v3);
					*/
					double a1_ = minAngle(v[v0], v[v1], v[v2]);
					double a2_ = minAngle(v[v0], v[v2], v[v3]);
					const double b1_ = std::min(a1_, a2_);
					const double b2_ = std::max(a1_, a2_);
					a1_ = minAngle(v[v1], v[v3], v[v0]);
					a2_ = minAngle(v[v1], v[v2], v[v3]);
					const double c1_ = std::min(a1_, a2_);
					const double c2_ = std::max(a1_, a2_);

					if (b1_ < c1_ || (b1_ == c1_ && b2_ <= c2_))
					{
						// first triangle
						tris.push_back(v1);
						tris.push_back(v3);
						tris.push_back(v0);
						// second triangle
						tris.push_back(v1);
						tris.push_back(v2);
						tris.push_back(v3);
					}
					else
					{
						tris.push_back(v0);
						tris.push_back(v1);
						tris.push_back(v2);
						// second triangle
						tris.push_back(v0);
						tris.push_back(v2);
						tris.push_back(v3);
					}
				}
			}
            /// <summary>
            /// compute the min angle for a traingle.
            /// </summary>
            /// <param name="v0">vertex 0</param>
            /// <param name="v1">vertex 1</param>
            /// <param name="v2">vertex 2</param>
			/// <returns>min angle</returns>
			double minAngle(const Vertex v0, const Vertex v1, const Vertex v2)
			{
				const float a = distance(v0, v1); // std::sqrt((v1.x - v0.x) * (v1.x - v0.x) + (v1.y - v0.y) * (v1.y - v0.y) + (v1.z - v0.z) * (v1.z - v0.z));
				const float b = distance(v1, v2); // std::sqrt((v2.x - v1.x) * (v2.x - v1.x) + (v2.y - v1.y) * (v2.y - v1.y) + (v2.z - v1.z) * (v2.z - v1.z));
				const float c = distance(v2, v0); // std::sqrt((v0.x - v2.x) * (v0.x - v2.x) + (v0.y - v2.y) * (v0.y - v2.y) + (v0.z - v2.z) * (v0.z - v2.z));
				const float A = std::acos((b * b + c * c - a * a) / (2 * b * c));
				const float B = std::acos((a * a + c * c - b * b) / (2 * a * c));
				const float C = std::acos((b * b + a * a - c * c) / (2 * b * a));

				return std::min(std::min(A, B), C);
			}
			// Mesh simplification
			using Halfedge = std::array<int, 5>;
			void halfedges(const int nr_v, std::vector<int>& quads, std::vector<Halfedge>& he, std::vector<int>& he_v, std::vector<int>& he_f);
			std::array<int, 4> collectNeighbors(const int quad, std::vector<Halfedge>& he, std::vector<int>& he_f);
			void colorFaces(std::vector<Halfedge>& he, std::vector<int>& he_f, std::vector<int>& colors);
			void vertexValence(const int nr_v, std::vector<Halfedge>& he, std::vector<int>& vV_);
			bool isNonManifold(const int quad, std::vector<Halfedge>& he, std::vector<int>& he_f);
			void mark3X3Y(std::vector<int>& quads, std::vector<int>& vV, std::vector<bool>& p3X3Y);
			void mergeVertices3X3Y(std::vector<Vertex>& v, std::vector<Normal>& normals, std::vector<bool> p3X3Y, std::vector<int>& vV, std::vector<int>& colors,
				std::vector<Halfedge>& he, std::vector<int>& he_f, std::vector<std::pair<bool, int>>& vm_, std::vector<bool>& em_);
			void mergeVertices3X3Y(std::vector<Vertex>& v, std::vector<Normal>& normals, std::vector<bool> p3X3Y, std::vector<int>& vV,
				std::vector<Halfedge>& he, std::vector<int>& he_f, std::vector<std::pair<bool, int>>& vm_, std::vector<bool>& em_);
			void removeVertices3X3Y(std::vector<Vertex>& v, std::vector<Normal>& n, std::vector<std::pair<bool, int>>& vm_, std::vector<Vertex>& nv, std::vector<Normal>& nn);
			void removeQuadrilaterals3X3Y(std::vector<int>& q, std::vector<bool>& em, std::vector<std::pair<bool, int>>& vm, std::vector<int>& nq);
			void simplify3X3Y(std::vector<Vertex>& v, std::vector<Normal>& n, std::vector<int>& quads, std::vector<int>& colors);

		};

} // namespace homotopy
