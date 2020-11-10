#pragma once
// Libs
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>

// utilities
#include "Timer.h"
// project files
#include "DualMarchingCubes.h"
#include "Volumes.h"
#include "UniformGrid.h"


namespace cpp_mc {
	class ImplicitSurface
	{
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
        // compute iso-surface using the DMC algorothm
		void dmc()
		{
			bool dType{ true };
			double i0{ 0 };
			std::cout << " ... generate volume data\n";
            m_nx = 256;
			m_ny = 256;
			m_nz = 256;
            // choose one surface case
			//m_volumes.scalar<Surface::GenusTwo>(m_ugrid, m_nx, m_ny, m_nz);
            m_volumes.scalar<Surface::iWP>(m_ugrid, m_nx, m_ny, m_nz);
            //m_volumes.scalar<Surface::MonkeySaddle>(m_ugrid, m_nx, m_ny, m_nz);
            //m_volumes.scalar<Surface::SternerRoman>(m_ugrid, m_nx, m_ny, m_nz);
            //m_volumes.scalar<Surface::Neovius>(m_ugrid, m_nx, m_ny, m_nz);
            // set the iso-value to be zero
			i0 = 0.0;
            // start computation
			std::cout << " ... compute DMC\n";
			Timer timer;
			timer.Start();
			std::vector<Vertex> v;
			std::vector<Normal> n;
			std::vector<int> tris;
			std::vector<int> quads;
			s_mc.dual_mc(i0, m_ugrid, v, n, tris, quads);
			timer.Stop();
			std::cout << " ... DMC: " << timer.GetMilisecondsElapsed() << " ms" << std::endl;
			std::cout << " ... Nr. vertices: " << v.size() << std::endl;
			std::cout << " ... Nr. quadrilaterals: " << (quads.size() / 4) << std::endl;
			std::cout << " ... write mesh to obj file\n";
			timer.Start();
			writeOBJ(v,n,tris,quads);
			timer.Stop();
			std::cout << " ... Mesh: " << timer.GetMilisecondsElapsed() << " ms" << std::endl;
		}

	private:
		UGrid m_ugrid;
		Volumes m_volumes;
		DualMarchingCubes s_mc;
		int m_nx{ 0 };
		int m_ny{ 0 };
		int m_nz{ 0 };
		double m_dx{ 0 };
		double m_dy{ 0 };
		double m_dz{ 0 };
		Timer timer;

	private:
		// write mesh into obj file
		void writeOBJ(std::vector<Vertex>& v_, std::vector<Normal>& n_, std::vector<int>& tris, std::vector<int>& quads)
		{
            std::ofstream objF;
            objF.open("./dualMC_result.obj");
            if (!objF.is_open()) {
                std::cout << "ERROR: can't open output file " << std::endl;
            }
            else {
                const int nr_v = static_cast<int>(v_.size());
                objF << "#Dual Marching Cubes\n";
                for (auto v : v_)
                {
                    objF << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
                }
                for (auto n : n_)
                {
                    objF << "vn " << n[0] << " " << n[1] << " " << n[2] << std::endl;
                }
                const int nr_q{ static_cast<int>(quads.size())/4};
                for (int f = 0; f < nr_q; f++)
                {
                    const int v0 = quads[4 * f];
                    const int v1 = quads[4 * f + 1];
                    const int v2 = quads[4 * f + 2];
                    const int v3 = quads[4 * f + 3];
                    //objF << "f " << (t[0] + 1) << "//" << (t[0] + 1) << " " << (t[1] + 1) << "//" << (t[1] + 1) << " " << (t[2] + 1) << "//" << (t[2] + 1) << std::endl;
                    objF << "f " << (v0 + 1) << "//" << (v0 + 1) << " " << (v1 + 1) << "//" << (v1 + 1) << " " << (v2 + 1) << "//" << (v2 + 1) << " " << (v3 + 1) << "//" << (v3 + 1) << std::endl;
                }
                objF.close();
            }
		}
	};
} // namespace homotopy
