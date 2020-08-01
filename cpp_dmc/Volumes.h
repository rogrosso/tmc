#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

#include "UniformGrid.h"

namespace cpp_mc {
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

		using Vertex = cpp_mc::Vector;
		using Point  = cpp_mc::Vector;
		using Normal = cpp_mc::Vector;
		using UGrid  = cpp_mc::UniformGrid;
		using Index  = UGrid::Index;
		using BBox = UGrid::BBox;
	public:
        // surface cases that can be computed with this class
		enum class Surface { Sphere, Torus, TwoHoledTorus, MonkeySaddle, GenusTwo, iWP, Neovius, SternerRoman };
	public:
        // read uniform grid from file
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
			ugrid.set_dx(dz);
			ugrid.set_dx(dy);

            size_t size_ = static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
            std::vector<double> v_data(size_);
            //ushort* t_buff = new ushort[size_];
			std::vector<T> t_buff(size_);
            ifile.read(reinterpret_cast<char*>(&t_buff[0]), size_ * sizeof(T));
            ifile.close();
            for (int k = 0; k < nz; k++)
            {
                for (int j = 0; j < ny; j++)
                {
                    {
#pragma omp parallel for
                        for (int i = 0; i < nx; i++)
                        {
                            ugrid.scalar(i, j, k, static_cast<double>(t_buff[k*ny*nx + j*nx + i]));
                        }
                    }
                }
            }
            // compute gradient for shading purpose
            ugrid.estimateGradient();
            ugrid.flip_gradient();
        }
        // computes the scalar values of the implicit function
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
			double x = minX;
			for (int i = 0; i < ugrid.x_size(); i++)
			{
				double y = minY;
				for (int j = 0; j < ugrid.y_size(); j++)
				{
					{
						double z = minZ;
						for (int k = 0; k < ugrid.z_size(); k++)
						{
							//ugrid.scalar(i, j, k, x * x + y * y + z * z);
							ugrid.scalar(i, j, k, surface<T>(x, y, z));
							z += dz;
						}
					}
					y += dy;
				}
				x += dx;
			}
			// compute gradient for shading purpose
			ugrid.estimateGradient();
			ugrid.flip_gradient();
		};
		template<Surface T>
		double surface(const double x, const double y, const double z)
		{
			return x * x + y * y + z * z;
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
	};

    template<>
    double Volumes::surface<Volumes::Surface::Sphere>(const double x, const double y, const double z)
    {
        return x * x + y * y + z * z;
    }
    template<>
    double Volumes::surface<Volumes::Surface::Torus>(const double x, const double y, const double z)
    {
        const double R = 0.6 * 0.6;
        const double r = 0.3 * 0.3;
        double val = (x * x + y * y + z * z + R - r);
        val = val * val;
        val = val - 4 * R * (x * x + y * y);
        return val;
    }
    template<>
    double Volumes::surface<Volumes::Surface::TwoHoledTorus>(const double x, const double y, const double z)
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
    double Volumes::surface<Volumes::Surface::MonkeySaddle>(const double x_, const double y_, const double z_)
    {
        const double alpha = 0.5;
        const double x = alpha * x_;
        const double y = alpha * y_;
        const double z = alpha * z_;
        return z - x * x * x - 3 * x * y * y;
    }
    template<>
    double Volumes::surface<Volumes::Surface::GenusTwo>(const double x_, const double y_, const double z_)
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
    double Volumes::surface<Volumes::Surface::iWP>(const double x_, const double y_, const double z_)
    {
        const float alpha = 5.01;
        //const float alpha = 1.01;
        const float x = alpha * (x_ + 1) * pi;
        const float y = alpha * (y_ + 1) * pi;
        const float z = alpha * (z_ + 1) * pi;
        return cos(x) * cos(y) + cos(y) * cos(z) + cos(z) * cos(x) - cos(x) * cos(y) * cos(z); // iso-value = 0
    }
    template<>
    double Volumes::surface<Volumes::Surface::Neovius>(const double x_, const double y_, const double z_)
    {
        const float alpha = 1;
        const float x = alpha * (x_ + 1) * pi;
        const float y = alpha * (y_ + 1) * pi;
        const float z = alpha * (z_ + 1) * pi;
        return 3 * (cos(x) + cos(y) + cos(z)) + 4 * cos(x) * cos(y) * cos(z); // iso_value = 0.0
    }
    template<>
    double Volumes::surface<Volumes::Surface::SternerRoman>(const double x_, const double y_, const double z_)
    {
        const float alpha = 1.5f;
        const float x = alpha * x_;
        const float y = alpha * y_;
        const float z = alpha * z_;
        auto sq = [](const double v) { return v * v;  };
        return sq(x * x + y * y + z * z - 1.0f) - (sq(z - 1) - 2.0f * x * x) * (sq(z + 1) - 2 * y * y);
    }

} // namespace cpp_mc
