#pragma once
#pragma once

// Libs
#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>

// project
#include "Vector.h"

namespace cpp_mc {

	class UniformGrid {
	public:
		using ushort = unsigned short;
		using Point = Vector;
		using Normal = Vector;
		using Index = std::array<int, 3>;
		using BBox = std::array<Point, 8>;
	public:
		void init(const std::string& filename);
		void init(const int nx, const int ny, const int nz);
		void init(const int nx, const int ny, const int nz, BBox& bb);
		void init(const int nx, const int ny, const int nz, BBox& bb, const double val);
		/// copy
		void copy(const UniformGrid& ug)
		{
			m_nx = ug.m_nx; //!< grid size in x-direction
			m_ny = ug.m_ny; //!< grid size in y-direction
			m_nz = ug.m_nz; //!< grid size in z-direction
			m_dx = ug.m_dx; //!< grid spacing in x-direction
			m_dy = ug.m_dy; //!< grid spacing in y-direction
			m_dz = ug.m_dz; //!< grid spacing in z-direction
			m_bbox = ug.m_bbox; //!< the bounding box of the ugrid.
			m_scalars = ug.m_scalars; //!< scalar values stored in the ugrid
			m_gradient = ug.m_gradient; //!< vector values stored in the ugrid
		}
		/// grid spacing in x-direction.
		double dx() { return m_dx; }
		double dx() const { return m_dx; }
		/// set grid spacing in x-direction
		void set_dx(const double d) { m_dx = d; }
		/// grid spacing in y-direction.
		double dy() { return m_dy; }
		double dy() const { return m_dy; }
		/// set grid spacing in y-direction
		void set_dy(const double d) { m_dy = d; }
		/// grid spacing in z-direction.
		double dz() { return m_dz; }
		double dz() const { return m_dz; }
		/// set grid spacing in z-direction
		void set_dz(const double d) { m_dz = d; }
		/// total number of grid points.
		int size() { return m_nx * m_ny * m_nz; }
		int size() const { return m_nx * m_ny * m_nz; }
		/// grid size in x-direction.
		int x_size() { return m_nx; }
		int x_size() const { return m_nx; }
		/// grid size in y-direction.
		int y_size() { return m_ny; }
		int y_size() const { return m_ny; }
		/// grid size in z-direction.
		int z_size() { return m_nz; }
		int z_size() const { return m_nz; }
		/// bounding box size
		std::array<Point, 8> bbox() { return m_bbox; }
		std::array<Point, 8> bbox() const { return m_bbox; }
		double minX() { return m_bbox[0][0]; }
		double minY() { return m_bbox[0][1]; }
		double minZ() { return m_bbox[0][2]; }
		double maxX() { return m_bbox[7][0]; }
		double maxY() { return m_bbox[7][1]; }
		double maxZ() { return m_bbox[7][2]; }

		/// Returns a point with the euclidean position of the vertex (i,j,k).
		/** This methods does not check if vertex is within the uniform grid
		*  It returns the vertex position for an infinite grid.
		*  @param[in] i cell index along x-coordinate
		*  @param[in] i cell index along x-coordinate
		*  @param[in] i cell index along x-coordinate
		*  @return _Point_ with coordinates of vertex {i,j,k}
		*/
		Point point(const int i, const int j, const int k) { return { m_bbox[0][0] + i * m_dx, m_bbox[0][1] + j * m_dy, m_bbox[0][2] + k * m_dz }; }
		Point point(const int i, const int j, const int k) const { return { m_bbox[0][0] + i * m_dx, m_bbox[0][1] + j * m_dy, m_bbox[0][2] + k * m_dz }; }
		Point point(const Index i) { return { m_bbox[0][0] + i[0] * m_dx, m_bbox[0][1] + i[1] * m_dy, m_bbox[0][2] + i[2] * m_dz }; }
		Point point(const Index i) const { return { m_bbox[0][0] + i[0] * m_dx, m_bbox[0][1] + i[1] * m_dy, m_bbox[0][2] + i[2] * m_dz }; }
		Point point(const int gl_index)
		{
			Index i = local_index(gl_index);
			return { m_bbox[0][0] + i[0] * m_dx, m_bbox[0][1] + i[1] * m_dy, m_bbox[0][2] + i[2] * m_dz };
		}
		/// Set all scalar to an input value
		void setScalars(const double val) { std::fill(m_scalars.begin(), m_scalars.end(), val); }
		/// Set the scalar value at grid node using node's global index.
		void scalar(const int gindex, const double val) { m_scalars[gindex] = val; }
		/// set the scalar value at grid node specified by indices (i,j,k).
		void scalar(const int i, const int j, const int k, const double val) { m_scalars[global_index(i, j, k)] = val; }
		/// set the scalar value at grid node specified by Index i.
		void scalar(const Index i, const double val) { m_scalars[global_index(i)] = val; }
		/// returns scalar value at grid node using global node index.
		double scalar(const int gindex) { return m_scalars[gindex]; }
		double scalar(const int gindex) const { return m_scalars[gindex]; }
		/// returns scalar value at grid node specified by indices (i,k,j).
		double scalar(const int i, const int j, const int k) { return m_scalars[global_index(i, j, k)]; }
		double scalar(const int i, const int j, const int k) const { return m_scalars[global_index(i, j, k)]; }
		/// returns scalar value at grid node specified by index i.
		double scalar(const Index i) { return m_scalars[global_index(i)]; }
		double scalar(const Index i) const { return m_scalars[global_index(i)]; }

		/// returns the normal vector at grid's node using node's global index.
		Vector gradient(const int gindex) { return m_gradient[gindex]; }
		const Vector gradient(const int gindex) const { return m_gradient[gindex]; }
		/// returns the normal vector at grid's node specified using i, j and k indices.
		Vector gradient(const int i, const int j, const int k) { return m_gradient[global_index(i, j, k)]; }
		const Vector gradient(const int i, const int j, const int k) const { return m_gradient[global_index(i, j, k)]; }
		Vector gradient(const Index i) { return m_gradient[global_index(i)]; }
		const Vector gradient(const Index i) const { return m_gradient[global_index(i)]; }
		/// invert normals
		void flip_gradient();
		/// check if a point is strictly inside the ugrid.
		bool in_bbox(const Point& p) {
			return  ((p[0] > m_bbox[0][0] && p[0] < m_bbox[1][0]) &&
				(p[1] > m_bbox[0][1] && p[1] < m_bbox[2][1]) &&
				(p[2] > m_bbox[0][2] && p[2] < m_bbox[4][2]));
		}
		/// compute index i,j,k of cell containing the point.
		Index cell_index(Point p)
		{
			Index id;
			id[0] = static_cast<int>((p[0] - m_bbox[0][0]) / m_dx);
			id[1] = static_cast<int>((p[1] - m_bbox[0][1]) / m_dy);
			id[2] = static_cast<int>((p[2] - m_bbox[0][2]) / m_dz);
			return id;
		}
		void cell(Point p, Point v[8])
		{
			const int i0 = static_cast<int>((p[0] - m_bbox[0][0]) / m_dx);
			const int j0 = static_cast<int>((p[1] - m_bbox[0][1]) / m_dy);
			const int k0 = static_cast<int>((p[2] - m_bbox[0][2]) / m_dz);
			v[0] = point(i0, j0, k0);
			v[1] = point(i0 + 1, j0, k0);
			v[2] = point(i0, j0 + 1, k0);
			v[3] = point(i0 + 1, j0 + 1, k0);
			v[4] = point(i0, j0, k0 + 1);
			v[5] = point(i0 + 1, j0, k0 + 1);
			v[6] = point(i0, j0 + 1, k0 + 1);
			v[7] = point(i0 + 1, j0 + 1, k0 + 1);
		}
		/// compute global index.
		int global_index(const Point p)
		{
			const int i = static_cast<int>((p[0] - m_bbox[0][0]) / m_dx);
			const int j = static_cast<int>((p[1] - m_bbox[0][1]) / m_dy);
			const int k = static_cast<int>((p[2] - m_bbox[0][2]) / m_dz);
			return k * m_ny*m_nx + j * m_nx + i;
		}
		int global_index(const int i, const int j, const int k) { return (k*m_ny*m_nx + j * m_nx + i); }
		int global_index(const int i, const int j, const int k) const { return (k*m_ny*m_nx + j * m_nx + i); }

		int global_index(const Index i) { return (i[2]*m_ny*m_nx + i[1] * m_nx + i[0]); }
		int global_index(const Index i) const { return (i[2] * m_ny*m_nx + i[1] * m_nx + i[0]); }
		Index local_index(const int g_index) { return Index{ g_index % m_nx,(g_index / m_nx) % m_ny, g_index / (m_nx*m_ny) }; }
		/// interpolate scalar
		//double interpolate_scalar(const Point& p);
		/// interplate normal
		//bool interpolate_normal(const Point& p, Normal& n);
		/// returns maximum scalar value in the grid
		double max_scalar() { return *std::max_element(m_scalars.begin(), m_scalars.end()); }
		double max_scalar() const { return *std::max_element(m_scalars.begin(), m_scalars.end()); }
		/// compute minimum value of scalar data stored on the grid
		double min_scalar() { return *std::min_element(m_scalars.begin(), m_scalars.end()); }
		double min_scalar() const { return *std::min_element(m_scalars.begin(), m_scalars.end()); }
		void estimateGradient();
		void normal(Normal& n, const int i, const int j, const int k, const double u, const double v, const double w)
		{
			const Vector n0 = gradient(i, j, k);
			const Vector n1 = gradient(i+1, j, k);
			const Vector n2 = gradient(i, j+1, k);
			const Vector n3 = gradient(i+1, j+1, k);
			const Vector n4 = gradient(i, j, k+1);
			const Vector n5 = gradient(i+1, j, k+1);
			const Vector n6 = gradient(i, j+1, k+1);
			const Vector n7 = gradient(i+1, j+1, k+1);
			n = (1 - w) * ((1 - v) * ((1 - u) * n0 + u * n1) + v * ((1 - u) * n2 + u * n3))
				+ w * ((1 - v) * ((1 - u) * n4 + u * n5) + v * ((1 - u) * n6 + u * n7));
			n.normalize();
		}
		void position(Point& p, const int i, const int j, const int k, const double u, const double v, const double w)
		{
			p[0] = m_bbox[0][0] + (i + u) * m_dx;
			p[1] = m_bbox[0][1] + (j + v) * m_dy;
			p[2] = m_bbox[0][2] + (k + w) * m_dz;
		}
		void writeVolume(std::string filename)
		{
			const size_t bytes{ m_scalars.size() };
			auto myfile = std::fstream(filename, std::ios::out | std::ios::binary);
			myfile.write((char*)&m_nx, sizeof(int));
			myfile.write((char*)&m_ny, sizeof(int));
			myfile.write((char*)&m_nz, sizeof(int));
			myfile.write((char*)&m_scalars[0], sizeof(double) * bytes);
			myfile.close();
		}

	private:
		int m_nx{ 0 }; //!< grid size in x-direction
		int m_ny{ 0 }; //!< grid size in y-direction
		int m_nz{ 0 }; //!< grid size in z-direction
		double m_dx{ 0 }; //!< grid spacing in x-direction
		double m_dy{ 0 }; //!< grid spacing in y-direction
		double m_dz{ 0 }; //!< grid spacing in z-direction
		std::array<Point, 8> m_bbox; //!< the bounding box of the ugrid.
		std::vector<double> m_scalars; //!< scalar values stored in the ugrid
		std::vector<Normal> m_gradient; //!< estimated gradient of scalar field

	public:
		/// print maximum scalar value stored in the ugrid.
		void print_max_scalar() { std::cout << "min scalar: " << *std::max_element(m_scalars.begin(), m_scalars.end()) << std::endl; }
		/// print minimum scalar valued stored in the ugrid.
		void print_min_scalar() { std::cout << "max scalar: " << *std::min_element(m_scalars.begin(), m_scalars.end()) << std::endl; }
		/// handle indices at grid boundaries using periodic conditions.
		int mod(int n, int m) { return (n >= 0) ? (n%m) : (m - (-n % m)) % m; }
		int mod(int n, int m) const { return (n >= 0) ? (n%m) : (m - (-n % m)) % m; }
		/// trilinar interpolation given local coordinates and vertex values
		//double trilinear(const double u, const double v, const double w, const std::array<double, 8>& F);
	};

} // homotopy
