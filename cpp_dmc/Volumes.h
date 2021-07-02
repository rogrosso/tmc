#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

#include "UniformGrid.h"

namespace dmc {
	/**
	 * Use boolean operation to construct implicite surfaces
	 * Given two sets A, B
	 * F(intersection(A,B)) = MAX(A,B)
	 * F(union(A,B)) = MIN(A,B)
	 * F(subtraction(A,B)) = MAX(A,-B)
	 */
	class Volumes {
	public:
		using uchar = unsigned char;
		using ushort = unsigned short;
		using uint = unsigned int;
		using Scalar = double;

		using Vertex = dmc::Vector;
		using Point  = dmc::Vector;
		using Normal = dmc::Vector;
		using UGrid  = dmc::UniformGrid;
		using Index  = UGrid::Index;
		using BBox = UGrid::BBox;
	    // Surfaces generated implicitly
		enum class Surface {
			Sphere,
			Torus,
			TwoHoledTorus,
			MonkeySaddle,
			GenusTwo,
			iWP,
			Neovius,
			SteinerRoman,
			Kummer,
			Tetrahedron
		};
	    // read volume from file
		template<typename T>
        void readFromFile(const std::string& i_file, UGrid& ugrid)
        {
            std::ifstream ifile;
            ifile.open(i_file, std::ios::binary);
            if (!ifile.is_open()) {
                exit(1);
            }
            int nx, ny, nz;
            float dx, dy, dz;
            ifile.read(reinterpret_cast<char*>(&nx), sizeof(int));
            ifile.read(reinterpret_cast<char*>(&ny), sizeof(int));
            ifile.read(reinterpret_cast<char*>(&nz), sizeof(int));
            ifile.read(reinterpret_cast<char*>(&dx), sizeof(float));
            ifile.read(reinterpret_cast<char*>(&dy), sizeof(float));
            ifile.read(reinterpret_cast<char*>(&dz), sizeof(float));

            double xmax = static_cast<double>(dx * (nx - 1));
            double ymax = static_cast<double>(dy * (ny - 1));
            double zmax = static_cast<double>(dz * (nz - 1));
            BBox bbox;
            bbox[0] = Point{ 0, 0, 0 };
            bbox[1] = Point{ xmax, 0, 0 };
            bbox[2] = Point{ 0, ymax, 0 };
            bbox[3] = Point{ xmax, ymax, 0 };
            bbox[4] = Point{ 0, 0, zmax };
            bbox[5] = Point{ xmax, 0, zmax };
            bbox[6] = Point{ 0, ymax, zmax };
            bbox[7] = Point{ xmax, ymax, zmax };
            ugrid.init(nx, ny, nz, bbox, 0);
            ugrid.set_dx(dx);
            ugrid.set_dy(dy);
            ugrid.set_dz(dz);

            size_t size_ = static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
            std::vector<double> v_data(size_);
            //ushort* t_buff = new ushort[size_];
            std::vector<T> t_buff(size_);
            ifile.read(reinterpret_cast<char*>(&t_buff[0]), size_ * sizeof(T));
            ifile.close();

#pragma omp parallel for
            for (int index = 0; index < static_cast<int>(size_); index++)
            {
                ugrid.scalar(index, static_cast<double>(t_buff[index]));
            }
            /*for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    for (int i = 0; i < nx; i++)
                    {
                        ugrid.scalar(i, j, k, static_cast<double>(t_buff[k*ny*nx + j*nx + i]));
                    }
                }
            }*/
            // compute gradient for shading purpose
            ugrid.estimateGradient();
            ugrid.flip_gradient();
        }
        // generate volume to compute implicit surface
        template<Surface T>
		void scalar(UGrid& ugrid, const int nx, const int ny, const int nz)
        {
            // center volume in [-1,1]^3
            initUGrid(ugrid, nx, ny, nz);
            const double minX = ugrid.minX();
            const double minY = ugrid.minX();
            const double minZ = ugrid.minX();
            const double dx = ugrid.dx();
            const double dy = ugrid.dy();
            const double dz = ugrid.dz();
            //double x = minX;

            const int size_ = nx * ny * nz;
#pragma omp parallel for
            for (int  g_index = 0; g_index < size_; g_index++)
            {
                const int i = g_index % nx;
                const int j = (g_index / nx) % ny;
                const int k = g_index / (nx * ny);
                const double x = minX + i * dx;
                const double y = minY + j * dy;
                const double z = minZ + k * dz;
                //ugrid.scalar(i, j, k, x * x + y * y + z * z);
                ugrid.scalar(i, j, k, surface<T>(x, y, z));
                Normal g;
                g[0] = (surface<T>(x + dx, y, z) - surface<T>(x - dx, y, z)) / (2 * dx);
                g[1] = (surface<T>(x, y + dy, z) - surface<T>(x, y - dy, z)) / (2 * dy);
                g[2] = (surface<T>(x, y, z + dz) - surface<T>(x, y, z - dz)) / (2 * dz);
                ugrid.gradient(i, j, k, g);
            }
            // compute gradient for shading purpose
            //ugrid.estimateGradient();
            //ugrid.flip_gradient();
        }
    private:
		template<Surface T>
		double surface(const double x, const double y, const double z)
		{
			return x * x + y * y + z * z;
		}
		template<>
		double surface<Surface::Sphere>(const double x, const double y, const double z)
		{
			return x * x + y * y + z * z;
		}
		template<>
		double surface<Surface::Torus>(const double x, const double y, const double z)
		{
			const double R = 0.6 * 0.6;
			const double r = 0.3 * 0.3;
			double val = (x * x + y * y + z * z + R - r);
			val = val * val;
			val = val - 4 * R * (x * x + y * y);
			return val;
		}
		template<>
		double surface<Surface::TwoHoledTorus>(const double x, const double y, const double z)
		{
			// center one torus at (-1/2,0,0), the other at (1/2,0,0)
			const double R = square(0.4);
			const double r = square(0.2);
			const double x1 = x + 0.4;
			const double x2 = x - 0.4;
			double val1 = square((square(x1) + square(y) + square(z) + R - r));
			val1 = val1 - 4 * R * (square(x1) + square(y));
			double val2 = square((square(x2) + square(y) + square(z) + R - r));
			val2 = val2 - 4 * R * (square(x2) + square(y));
			return std::min(val1, val2);
		}
		template<>
		double surface<Surface::MonkeySaddle>(const double x_, const double y_, const double z_)
		{
			const double alpha = 0.5;
			const double x = alpha * x_;
			const double y = alpha * y_;
			const double z = alpha * z_;
			return z - x * x * x - 3 * x * y * y;
		}
		template<>
		double surface<Surface::GenusTwo>(const double x_, const double y_, const double z_)
		{
			double alpha = 1.0;
			double x = (x_ + 1.0) / 2.0;
			double y = (y_ + 1.0) / 2.0;
			double z = (z_ + 1.0) / 2.0;
			x = alpha * (4 * x - 2);
			y = alpha * (4 * y - 2);
			z = alpha * (4 * z - 2);
			double val = 2 * y * (y * y - 3 * x * x) * (1 - z * z) + (x * x + y * y) * (x * x + y * y) - (9 * z * z - 1) * (1 - z * z);
			return val;
		}
		template<>
		double surface<Surface::iWP>(const double x_, const double y_, const double z_)
		{
			const double alpha = 5.01;
			//const float alpha = 1.01;
			const double x = alpha * (x_ + 1) * pi;
			const double y = alpha * (y_ + 1) * pi;
			const double z = alpha * (z_ + 1) * pi;
			return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z); // iso-value = 0
		}
		template<>
		double surface<Surface::Neovius>(const double x_, const double y_, const double z_)
		{
			const double alpha = 2;
			const double x = alpha * (x_ + 1) * pi;
			const double y = alpha * (y_ + 1) * pi;
			const double z = alpha * (z_ + 1) * pi;
			return 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z); // iso_value = 0.0
		}
		template<>
		double surface<Surface::SteinerRoman>(const double x_, const double y_, const double z_)
		{
			const double alpha = 1.f;
			//const float r = 1.5f;
			const double x = alpha * x_;
			const double y = alpha * y_;
			const double z = alpha * z_;
			auto sq = [](const double v) { return v * v;  };
			return sq(x * x + y * y + z * z - 1.0f) - (sq(z - 1) - 2.0f * x * x) * (sq(z + 1) - 2 * y * y);
			//return sq(x * y) + sq(x * z) + sq(y * z) - r * x * y * z;
		}
		template<>
		double surface<Surface::Kummer>(const double x_, const double y_, const double z_)
		{
			const double alpha{ 2 };
			const double x{ alpha * x_ };
			const double y{ alpha * y_ };
			const double z{ alpha * z_ };
			const double mu{ 1.3 };
			const double lambda{ (3 * mu * mu - 1) / (3 - mu * mu) };
			const double w2{ std::sqrt(2) };
			const double p = 1 - z - w2 * x;
			const double q = 1 - z + w2 * x;
			const double r = 1 + z + w2 * y;
			const double s = 1 + z - w2 * y;
			double v = (x * x + y * y + z * z - mu * mu);
			v = v * v;
			v = v - lambda * p * q * r * s;
			return v;
		}
		template<>
		double surface<Surface::Tetrahedron>(const double x_, const double y_, const double z_)
		{
			// set outside 0, inside 1
			double val{ 0 };
			Vertex v0{ -0.8,-0.8,-0.8 };
			Vertex v1{  0.8,-0.8,-0.8 };
			Vertex v2{  0.8, 0.8,-0.8 };
			Vertex v3{  0.0, 0.0, 0.8 };
			Vertex p{ x_,y_,z_ };

			double d0 = distancePointTriangle(v0, v2, v1, p);
			double d1 = distancePointTriangle(v0, v1, v3, p);
			double d2 = distancePointTriangle(v1, v2, v3, p);
			double d3 = distancePointTriangle(v0, v3, v2, p);
			double d = std::min(std::fabs(d0), std::fabs(d1));
			d = std::min(d, std::fabs(d2));
			d = std::min(d, std::fabs(d3));
			if (d == std::fabs(d0)) return d0;
			else if (d == std::fabs(d1)) return d1;
			else if (d == std::fabs(d2)) return d2;
			else if (d == std::fabs(d3)) return d3;
			return d;

		}


	private:
		void initUGrid(UGrid& ugrid, const int nx, const int ny, const int nz)
		{
			BBox bb;
			bb[0] = { -1,-1,-1 };
			bb[1] = {  1,-1,-1 };
			bb[2] = { -1, 1,-1 };
			bb[3] = {  1, 1,-1 };
			bb[4] = { -1,-1, 1 };
			bb[5] = {  1,-1, 1 };
			bb[6] = { -1, 1, 1 };
			bb[7] = {  1, 1, 1 };
			ugrid.init(nx, ny, nz, bb);
		}
		double square(const double x) { return x * x; }
		double pi{ 3.14159265358979323846 };
		// distance point to edge
		double distancePointEdge(const Vertex& v0, const Vertex& v1, const Vertex& p)
		{
			Vertex n = v1 - v0;
			n.normalize();
			const double l = dot(p - v0, n);
			if (l <= 0 || l >= distance(v0,v1))
			{
				// point is outside, measure distance to vertex
				const double l0 = distance(p, v0);
				const double l1 = distance(p, v1);
				if (l0 < l1) return l0;
				else return l1;
			}
			else
			{
				return distance(p, v0 + l * n);
			}
		}
		// distace point to triangle
		double distancePointTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2, const Vertex& p)
		{
			// project point to triangle's plane
			Vertex n = cross(v1 - v0, v2 - v0);
			n.normalize();
			const double l = -dot(p - v0, n);
			Vertex q = p + l * n;
			// compute generalized barycentric coords of q
			const double d0 = dot(n, cross(v0 - q, v1 - q));
			const double d1 = dot(n, cross(v1 - q, v2 - q));
			const double d2 = dot(n, cross(v2 - q, v0 - q));
			if (d0 > 0 && d1 > 0 && d2 > 0)
			{
				// point is within triangle
				return dot(n,p-v0);
			}
			else {
				// point is outside triangle, measure distance to edge
				const double e0 = distancePointEdge(v0, v1, p);
				const double e1 = distancePointEdge(v1, v2, p);
				const double e2 = distancePointEdge(v2, v0, p);
				return std::min(std::min(e0, e1), e2);
			}
		}
	};
} // namespace homotopy
