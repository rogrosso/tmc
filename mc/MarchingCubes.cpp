//  MarchingCubes.cpp
//  Quitte
//
//  Created by Roberto Grosso on 05.06.16.
//  Copyright © 2016 Roberto Grosso. All rights reserved.
//

#include "MarchingCubes.h"


void
tmc::MarchingCubes::operator() (const double i0, const std::string& i_file,const bool c_flag, const bool obj_flag, const std::string& o_objF,const bool off_flag,const std::string& o_offF)
{
    std::cout << " ... reading data \n";

	/*std::FILE* f{ nullptr };
	errno_t status = fopen_s(&f, i_file.c_str(), "rb");
	if (status != 0) {
		std::cerr << "ERROR: can't open file " << i_file.c_str() << std::endl;
		exit(1);
	}*/
	std::FILE* f = fopen(i_file.c_str(), "rb");
	if (f == nullptr) {
		std::cout << "ERROR: can't open input file " << i_file << std::endl;
		exit(1);
	}
    ushort x_size;
    ushort y_size;
    ushort z_size;
    std::fread(&x_size, sizeof(ushort), 1, f);
    std::fread(&y_size, sizeof(ushort), 1, f);
    std::fread(&z_size, sizeof(ushort), 1, f);
    float dx;
    float dy;
    float dz;
    std::fread(&dx, sizeof(float), 1, f);
    std::fread(&dy, sizeof(float), 1, f);
    std::fread(&dz, sizeof(float), 1, f);
    m_nx = x_size;
    m_ny = y_size;
    m_nz = z_size;
    m_dx = double(dx);
    m_dy = double(dy);
    m_dz = double(dz);
    double xmax = m_dx * (m_nx - 1.);
    double ymax = m_dy * (m_ny - 1.);
    double zmax = m_dz * (m_nz - 1.);
    std::array<Point, 8> bb;
	bb[0] = Point{ { 0, 0, 0 } };
	bb[1] = Point{ { xmax, 0, 0 } };
	bb[2] = Point{ { 0, ymax, 0 } };
	bb[3] = Point{ { xmax, ymax, 0 } };
	bb[4] = Point{ { 0, 0, zmax } };
	bb[5] = Point{ { xmax, 0, zmax } };
	bb[6] = Point{ { 0, ymax, zmax } };
	bb[7] = Point{ { xmax, ymax, zmax } };

    m_ugrid.init(m_nx, m_ny, m_nz,bb);
    int m_size = m_nx * m_ny * m_nz;
    unsigned short* t_buff = new unsigned short[x_size*y_size*z_size];
    std::fread(&t_buff[0],sizeof(unsigned short),x_size*y_size*z_size,f);
    for (int i = 0; i < m_size; i++) {
        m_ugrid.scalar(i,(double)t_buff[i]);
    }
    std::fclose(f);
    delete [] t_buff;



    // compute gradient for normals
    m_ugrid.gradient();
   
    // invert gradient
    //m_ugrid.invert_normals();
   
    // compute isosurface
    std::cout << " ... computing isosurface\n";
	//const double i0 = 162.61;
	//const double i0 = 1003.4;
	//const double i0 = 1233.6;
    t_mc(i0);


    const int nr_v = (int)m_points.size();
    const int nr_t = (int)m_triangles.size();

    std::cout << "tot. nr. of triangles: " << nr_t << std::endl;
    std::cout << "tot. nr. of vertices:  " << nr_v << std::endl;


    // Check mesh topology and correct triangle orientation if necessary
	if (c_flag) {
		std::cout << " ... check and correct triangle orientation\n";
		connectivity();
	}
	

    // write obj
	if (obj_flag) {
		std::cout << " ... write obj file\n";
		std::ofstream objF;
		objF.open(o_objF.c_str());
		if (objF.is_open()) {
			objF << "# Topologically correct and manifold isosurface\n";
			for (int i = 0; i < nr_v; i++) {
				objF << "v " << std::setprecision(7) << std::fixed << m_points[i][0] << " " << m_points[i][1] << " " << m_points[i][2] << std::endl;
			}
			for (int n = 0; n < nr_v; n++) {
				objF << "vn " << std::setprecision(7) << std::fixed << m_pnorms[n][0] << " " << m_pnorms[n][1] << " " << m_pnorms[n][2] << std::endl;
			}
			for (auto t : m_triangles) {
				objF << "f " << (t.v[0] + 1) << "//" << (t.v[0] + 1) << " " << (t.v[1] + 1) << "//" << (t.v[1] + 1) << " " << (t.v[2] + 1) << "//" << (t.v[1] + 1) << std::endl;
			}
			objF.close();
		}
		else {
			std::cout << "ERROR: can't open output file " << o_objF << std::endl;
		}
	}
	

    // write obj output file
	if (off_flag) {
		std::cout << " ... write OFF file: " << o_objF << std::endl;
		std::ofstream offF;
		offF.open(o_offF.c_str());
		if (!offF.is_open()) {
			std::cout << "ERROR: can't open output file " << o_objF << std::endl;
		}
		else {
			offF << "OFF\n";
			offF << nr_v << " " << nr_t << " " << 0 << std::endl;
			for (int i = 0; i < nr_v; i++) {
				offF << m_points[i][0] << " " << m_points[i][1] << " " << m_points[i][2] << std::endl;
			}
			for (auto t : m_triangles) {
				offF << "3 " << t.v[0] << " " << t.v[1] << " " << t.v[2] << std::endl;
			}
			offF.close();
		}
	}
	return;
}



//*******************************************************************************************************************************************
//  IMPLEMENTATION t_mc
//*******************************************************************************************************************************************
void
tmc::MarchingCubes::t_mc(const double i0)
{
    // edges are uniquely characterized by the two end vertices, which have a unique vertex id
    // the end vertices of the edge are computed in the cell by giving the indices (i,j,k).
    // These indices are obtained from the cell index by adding 0 or 1 to i, j or k respectively
    // Example: edge 0: (i,j,k) - (i+1,j,k)
    //          edge 1: (i+1,j,k) - (i+1,j+1,k)
    // The first 3 indices are for the first vertex and the second 3 for the second vertex.
    // there are 12 edges, assign to each vertex three edges, the global edge numbering
    // consist of 3*global_vertex_id + edge_offset.
    const int global_edge_id[][4] = {{0,0,0,0},{1,0,0,1},{0,1,0,0},{0,0,0,1},
        {0,0,1,0},{1,0,1,1},{0,1,1,0},{0,0,1,1},
        {0,0,0,2},{1,0,0,2},{1,1,0,2},{0,1,0,2}};
    // the end vertices of an edge
    int l_edges[12][2] = {{0,1}, {1,3}, {2,3}, {0,2},
        {4,5}, {5,7}, {6,7}, {4,6},
        {0,4}, {1,5}, {3,7}, {2,6}};
    // compute sizes
    const int nx = m_ugrid.x_size();
    const int ny = m_ugrid.y_size();
    const int nz = m_ugrid.z_size();

    // we need to compute up to 3 vertices at the interior of a cell, therefore
    // the cell shift factor is set to 3+3 = 6, i.e. 3 edges assigned to a cell for global numberig
    // and 3 vertices in the interior of the cell
    m_cell_shift_factor = 6;

    // there can be at most 12 intersections
    std::vector<Vertex> vertices(12);
    std::vector<Point>  ip(12);
    std::vector<Normal> in(12);

    timer.start();
    int i_case_count = 0;
    // marching cubes
    for (int k = 0; k < (nz-1); k++) {
        m_kindex = k;
        for (int j = 0; j < (ny-1); j++) {
            m_jindex = j;
            for (int i = 0; i < (nx-1); i++) {
                m_iindex = i;
                // slice hex
                // collect function values and build index
                double u[8];
                Point  p[8];
                Normal n[8];
                int vi{0};
                std::bitset<8> index = 0;
                for (int kl = 0; kl <= 1; kl++) {
                    for (int jl = 0; jl <= 1; jl++) {
                        for (int il = 0; il <= 1; il++) {
                            // collect scalar values and computex index
                            p[vi] = m_ugrid.point(i+il,j+jl,k+kl);
                            u[vi] = m_ugrid.scalar(i+il,j+jl,k+kl);
                            if (u[vi] >= i0) {
                                //index.set(VertexMapping[vi]);
                                index.set(vi);

                            }
                            // probably better get normals here
                            n[vi] = m_ugrid.normal(i+il, j+jl, k+kl);
                            // next cell vertex
                            vi++;
                        }
                    }
                }

                // collect edges from table and
                // interpolate triangle vertex positon
                int i_case = int(index.to_ullong());
                //int tcm = (int)t_entry[i_case];
                int tcm = (int)t_ambig[i_case];
                if (tcm == 105) {
                    i_case_count++;
                    //t_slice(i,j,k,i0,u,p,n,i_case);
                    p_slice(i,j,k,i0,u,p,n,i_case);
                } else {
                    // compute for this case the vertices
                    ushort flag = 1;
                    for (int eg = 0; eg < 12; eg++) {
                        if (flag & e_pattern[i_case]) {
                            // the edge global index is given by the vertex global index + the edge offset
                            const int ix = i + global_edge_id[eg][0];
                            const int iy = j + global_edge_id[eg][1];
                            const int iz = k + global_edge_id[eg][2];
                            vertices[eg].g_edg = uint(m_cell_shift_factor*m_ugrid.global_index(ix, iy, iz) + global_edge_id[eg][3]);
                            // generate vertex here, do not care at this point if vertex already exist
                            int* vert = l_edges[eg];
                            // interpolation weight
                            const int v0 = vert[0];
                            const int v1 = vert[1];
                            double l = (i0 - u[v0]) / (u[v1] - u[v0]);
                            // interpolate vertex
                            ip[eg][0] = (1 - l)*p[v0][0] + l*p[v1][0];
                            ip[eg][1] = (1 - l)*p[v0][1] + l*p[v1][1];
                            ip[eg][2] = (1 - l)*p[v0][2] + l*p[v1][2];

                            // interpolate normal
                            in[eg][0] = (1 - l)*n[v0][0] + l*n[v1][0];
                            in[eg][1] = (1 - l)*n[v0][1] + l*n[v1][1];
                            in[eg][2] = (1 - l)*n[v0][2] + l*n[v1][2];
                            const double nlength = std::sqrt(in[eg][0]*in[eg][0]+in[eg][1]*in[eg][1]+in[eg][2]*in[eg][2]);
                            in[eg][0] = in[eg][0] / nlength;
                            in[eg][1] = in[eg][1] / nlength;
                            in[eg][2] = in[eg][2] / nlength;

                            // set vertex index
                            auto s_index = m_vertices.find(vertices[eg].g_edg);
                            if (s_index == m_vertices.end()) {
                                // index not found!
                                const int g_idx = (int)m_points.size();
                                m_vertices[vertices[eg].g_edg] = g_idx;
                                vertices[eg].g_idx = g_idx;
                                m_points.push_back(ip[eg]);
                                m_pnorms.push_back(in[eg]);
                            } else {
                                vertices[eg].g_idx = s_index->second; // this is vertex global index g_idx
                            }
                        }
                        flag <<=1;
                    }

                    // construct triangles
                    for (int t = 0; t < 16; t += 3) {
                        const int t_index = i_case * 16 + t;
                        //if (e_tris_list[t_index] == 0x7f)
                        if (t_pattern[t_index] == -1)
                            break;
                        Triangle tri;
                        //                        const int eg0 = e_tris_list[t_index];
                        //                        const int eg1 = e_tris_list[t_index+1];
                        //                        const int eg2 = e_tris_list[t_index+2];
                        const int eg0 = t_pattern[t_index];
                        const int eg1 = t_pattern[t_index+1];
                        const int eg2 = t_pattern[t_index+2];
						// change MC triangle orientation
						// positive vertices must be outside the surface
						tri.v[0] = (int)vertices[eg2].g_idx;
						tri.v[1] = (int)vertices[eg1].g_idx;
						tri.v[2] = (int)vertices[eg0].g_idx;

                        // insert new triangle in list
                        m_triangles.push_back(tri);
                    }
                }
            }
        }
    }

    timer.stop();
    timer.print();

    std::cout << "tot. nr. of vertices:        " << m_points.size()    << std::endl;
    std::cout << "tot. nr. of triangles:       " << m_triangles.size() << std::endl;
    std::cout << "tot. nr. of specials cases:  " << i_case_count << std::endl;
    std::cout << "tot. nr cases 3: " << m_ccases_3 << std::endl;
    std::cout << "tot. nr cases 4: " << m_ccases_4 << std::endl;
    std::cout << "tot. nr cases 5: " << m_ccases_5 << std::endl;
    std::cout << "tot. nr cases 6: " << m_ccases_6 << std::endl;
    std::cout << "tot. nr cases 6a: " << m_ccases_6a << std::endl;
    std::cout << "tot. nr cases 7: " << m_ccases_7 << std::endl;
    std::cout << "tot. nr cases 8: " << m_ccases_8 << std::endl;
    std::cout << "tot. nr cases 9: " << m_ccases_9 << std::endl;
    std::cout << "tot. nr cases tunnel: " << m_ccases_tunnel << std::endl;
    std::cout << "tot. nr cases vert12: " << m_ccases_12cont << std::endl;
}


//*******************************************************************************************************************************************
//  IMPLEMENTATION p_slice
//*******************************************************************************************************************************************
void
tmc::MarchingCubes::p_slice(const int i_index, const int j_index, const int k_index, const double i0, double* F, Point* p, Normal* n, const int i_case)
{
	// there are 12 edges, assign to each vertex three edges, the global edge numbering
	// consist of 3*global_vertex_id + edge_offset.
	const unsigned long long gei_pattern_ = 670526590282893600ull;
	const uint axis_mask = 1;
	const uint offs_mask = 3;

	// code edge end vertices for each of the 12 edges
	const unsigned char l_edges_[12] = { 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
	auto get_edge_vertex = [](const int e, unsigned int& v0, unsigned int& v1, const unsigned char l_edges_[12]) {
		v0 = (unsigned int)(l_edges_[e] & 0xF);
		v1 = (unsigned int)(l_edges_[e] >> 4) & 0xF;
	};

	// A hexahedron has twelve edges, save the intersection of the isosurface with the edge
	// save global edge and global vertex index of isosurface
	std::vector<Vertex> vertices(12);
	// save coordinates of intersection points in 3D
	std::vector<Point>  ip(12);
	// save normals of intersection points from scalar data
	std::vector<Normal> in(12);
	// save loca coordinate along the edge of intersection point
	std::vector<double> ecoord{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// collect vertices
	ushort   flag{ 1 };
	for (int eg = 0; eg < 12; eg++) {
		vertices[eg].g_edg = -1;
		vertices[eg].g_idx = -1;
		if (flag & e_pattern[i_case]) {
			// the edge global index is given by the vertex global index + the edge offset
			uint shift = 5 * eg;
			const int ix = i_index + (int)((gei_pattern_ >> shift) & 1); // global_edge_id[eg][0];
			const int iy = j_index + (int)((gei_pattern_ >> (shift + 1)) & 1); // global_edge_id[eg][1];
			const int iz = k_index + (int)((gei_pattern_ >> (shift + 2)) & 1); // global_edge_id[eg][2];
			const int off_val = (int)((gei_pattern_ >> (shift + 3)) & 3);

			vertices[eg].g_edg = int(m_cell_shift_factor*m_ugrid.global_index(ix, iy, iz) + off_val);

			// generate vertex here, do not care at this point if vertex already exist
			//int* vert = l_edges[eg];
			// interpolation weight
			//uint v0 = (l_edges_[eg]&0xF);
			//uint v1 = (l_edges_[eg]>>4)&0xF;
			//            int v0 = vert[0];
			//            int v1 = vert[1];
			uint v0, v1;
			get_edge_vertex(eg, v0, v1, l_edges_);

			double l = (i0 - F[v0]) / (F[v1] - F[v0]);
			ecoord[eg] = l;
			// interpolate vertex
			ip[eg][0] = (1 - l)*p[v0][0] + l*p[v1][0];
			ip[eg][1] = (1 - l)*p[v0][1] + l*p[v1][1];
			ip[eg][2] = (1 - l)*p[v0][2] + l*p[v1][2];

			// interpolate normal
			in[eg][0] = (1 - l)*n[v0][0] + l*n[v1][0];
			in[eg][1] = (1 - l)*n[v0][1] + l*n[v1][1];
			in[eg][2] = (1 - l)*n[v0][2] + l*n[v1][2];

			const double n_size = std::sqrt(in[eg][0] * in[eg][0] + in[eg][1] * in[eg][1] + in[eg][2] * in[eg][2]);
			in[eg][0] /= n_size;
			in[eg][1] /= n_size;
			in[eg][2] /= n_size;

			// set vertex in map
			// set vertex index
			auto s_index = m_vertices.find(vertices[eg].g_edg);
			if (s_index == m_vertices.end()) {
				const int g_idx = (int)m_points.size();
				vertices[eg].g_idx = g_idx;
				m_vertices[vertices[eg].g_edg] = g_idx;
				m_points.push_back(ip[eg]);
				m_pnorms.push_back(in[eg]);
			}
			else {
				vertices[eg].g_idx = s_index->second;
			}
		}
		/*else {
		e_set[eg] = false;
		}*/
		//next edge
		flag <<= 1;
	}

	// compute oriented contours
	// A countour consists of segment at the faces connecting the intersection of the
	// iso-surface with the edges. For each edge we store the edge to which the segment
	// is outgoing and the edge from which the segment in comming. Therefore a contour
	// cab be reconstructed by connecting the edges in the direccion of the outgoing.
	// The contour is oriented in such a way, that the positive vertices are outside.
	// 1. build segments
	// 2. connect segments
	// build up segments
	// set segments map
	unsigned char segm_[12] = { 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
	auto set_segm = [](const int e, const int pos, const int val, unsigned char segm_[12]) {
		if (pos == 0) {
			segm_[e] &= 0xF0;
			segm_[e] |= (unsigned char)val & 0xF;
		}
		else if (pos == 1) {
			segm_[e] &= 0xF;
			segm_[e] |= val << 4;
		}
	};
	auto get_segm = [](const int e, const int pos, unsigned char segm_[12]) {
		if (pos == 0)
			return (int)(segm_[e] & 0xF);
		else
			return (int)((segm_[e] >> 4) & 0xF);
	};
	auto is_segm_set = [](const int e, unsigned char segm_[12]) {
		return (segm_[e] != 0xFF);
	};
	auto unset_segm = [](const int e, unsigned char segm_[12]) {
		segm_[e] = 0xFF;
	};
	// In order to compute oriented segments, the hexahedron has to be flatten.
	// The insides of the faces of the hexahedron have to be all at the same
	// side of the flattend hexa. This requires changing the order of the
	// edges when reading from the faces
	// code edges at face
	//unsigned short face_e_[6] = { 12816, 30292, 33936, 46754, 34739, 38305 };
	std::array<unsigned short, 6> e_face_{ { 291, 18277, 18696, 10859, 33719, 38305 } };
	// code vertices at face
	//unsigned short face_v_[6] = { 12816, 30292, 21520, 30258, 25632, 30001 };
	std::array<unsigned short, 6> v_face_{ { 12576, 25717, 5380, 29538, 8292, 30001 } };

	// reading edge from face
	auto get_face_e = [e_face_](const int f, const int e) { return ((e_face_[f] >> (4 * e)) & 0xF); };
	auto get_face_v = [v_face_](const int f, const int e) { return ((v_face_[f] >> (4 * e)) & 0xF); };
	// compute oriented segments using the isoline scheme at the faces
	const unsigned int BIT_1 = 1;
	const unsigned int BIT_2 = 2;
	const unsigned int BIT_3 = 4;
	const unsigned int BIT_4 = 8;
	auto asymptotic_decider = [](const double f0, const double f1, const double f2, const double f3) {
		return (f0*f3 - f1*f2) / (f0 + f3 - f1 - f2);
	};
	std::vector<bool> f_flag(6, false);
	for (int f = 0; f < 6; f++) {
		// classify face
		unsigned int f_case{ 0 };
		uint v0 = get_face_v(f, 0);
		uint v1 = get_face_v(f, 1);
		uint v2 = get_face_v(f, 2);
		uint v3 = get_face_v(f, 3);
		uint e0 = get_face_e(f, 0);
		uint e1 = get_face_e(f, 1);
		uint e2 = get_face_e(f, 2);
		uint e3 = get_face_e(f, 3);
		double f0 = F[v0];
		double f1 = F[v1];
		double f2 = F[v2];
		double f3 = F[v3];
		if (f0 >= i0)
			f_case |= BIT_1;
		if (f1 >= i0)
			f_case |= BIT_2;
		if (f2 >= i0)
			f_case |= BIT_3;
		if (f3 >= i0)
			f_case |= BIT_4;
		switch (f_case)
		{
		case 1:
			set_segm(e0, 0, e3, segm_);
			set_segm(e3, 1, e0, segm_);
			break;
		case 2:
			set_segm(e1, 0, e0, segm_);
			set_segm(e0, 1, e1, segm_);
			break;
		case 3:
			set_segm(e1, 0, e3, segm_);
			set_segm(e3, 1, e1, segm_);
			break;
		case 4:
			set_segm(e3, 0, e2, segm_);
			set_segm(e2, 1, e3, segm_);
			break;
		case 5:
			set_segm(e0, 0, e2, segm_);
			set_segm(e2, 1, e0, segm_);
			break;
		case 6:
		{
			const double val = asymptotic_decider(f0, f1, f2, f3);
			if (val > i0) {
				set_segm(e3, 0, e0, segm_);
				set_segm(e0, 1, e3, segm_);
				set_segm(e1, 0, e2, segm_);
				set_segm(e2, 1, e1, segm_);
			}
			else if (val < i0) {
				set_segm(e1, 0, e0, segm_);
				set_segm(e0, 1, e1, segm_);
				set_segm(e3, 0, e2, segm_);
				set_segm(e2, 1, e3, segm_);
			}
			else {
				f_flag[f] = true;
				// singular case val == i0, there are no asymptotes
				// check if there is a reasonable triangulation of the face
				unsigned short e_flag = e_flag = 0x218;
				unsigned short bit_1 = 0x1;
				unsigned short bit_2 = 0x2;
				double ec0 = ecoord[e0];
				double ec1 = ecoord[e1];
				double ec2 = ecoord[e2];
				double ec3 = ecoord[e3];
				if ((e_flag >> (f * 2)) & bit_1) {
					ec0 = 1 - ec0;
					ec2 = 1 - ec2;
				}
				if ( (e_flag >> (f * 2)) & bit_2) {
					ec1 = 1 - ec1;
					ec3 = 1 - ec3;
				}
				if (ec1 < ec3 && ec0 > ec2) {
					set_segm(e1, 0, e0, segm_);
					set_segm(e0, 1, e1, segm_);
					set_segm(e3, 0, e2, segm_);
					set_segm(e2, 1, e3, segm_);
				}
				else if (ec1 > ec3 && ec0 < ec2) {
					set_segm(e3, 0, e0, segm_);
					set_segm(e0, 1, e3, segm_);
					set_segm(e1, 0, e2, segm_);
					set_segm(e2, 1, e1, segm_);
				}
				else {
					std::cerr << "ERROR: can't correctly triangulate cell's face\n";
					return;
				}
			}
		}
		break;
		case 7:
			set_segm(e1, 0, e2, segm_);
			set_segm(e2, 1, e1, segm_);
			break;
		case 8:
			set_segm(e2, 0, e1, segm_);
			set_segm(e1, 1, e2, segm_);
			break;
		case 9:
		{
			const double val = asymptotic_decider(f0, f1, f2, f3);
			if (val > i0){
				set_segm(e0, 0, e1, segm_);
				set_segm(e1, 1, e0, segm_);
				set_segm(e2, 0, e3, segm_);
				set_segm(e3, 1, e2, segm_);
			}
			else if (val < i0) {
				set_segm(e0, 0, e3, segm_);
				set_segm(e3, 1, e0, segm_);
				set_segm(e2, 0, e1, segm_);
				set_segm(e1, 1, e2, segm_);
			}
			else {
				f_flag[f] = true;
				// singular case val == i0, there are no asymptotes
				// check if there is a reasonable triangulation of the face
				unsigned short e_flag = e_flag = 0x218;
				unsigned short bit_1 = 0x1;
				unsigned short bit_2 = 0x2;
				double ec0 = ecoord[e0];
				double ec1 = ecoord[e1];
				double ec2 = ecoord[e2];
				double ec3 = ecoord[e3];
				if ((e_flag >> (f * 2)) & bit_1) {
					ec0 = 1 - ec0;
					ec2 = 1 - ec2;
				}
				if ((e_flag >> (f * 2)) & bit_2) {
					ec1 = 1 - ec1;
					ec3 = 1 - ec3;
				}
				if (ec1 < ec3 && ec0 > ec2) {
					set_segm(e0, 0, e1, segm_);
					set_segm(e1, 1, e0, segm_);
					set_segm(e2, 0, e3, segm_);
					set_segm(e3, 1, e2, segm_);
				}
				else if (ec1 > ec3 && ec0 < ec2) {
					set_segm(e0, 0, e3, segm_);
					set_segm(e3, 1, e0, segm_);
					set_segm(e2, 0, e1, segm_);
					set_segm(e1, 1, e2, segm_);
				}
				else {
					std::cerr << "ERROR: can't correctly triangulate cell's face\n";
					return;
				}
			}
		}
		break;
		case 10:
			set_segm(e2, 0, e0, segm_);
			set_segm(e0, 1, e2, segm_);

			break;
		case 11:
			set_segm(e2, 0, e3, segm_);
			set_segm(e3, 1, e2, segm_);

			break;
		case 12:
			set_segm(e3, 0, e1, segm_);
			set_segm(e1, 1, e3, segm_);

			break;
		case 13:
			set_segm(e0, 0, e1, segm_);
			set_segm(e1, 1, e0, segm_);

			break;
		case 14:
			set_segm(e3, 0, e0, segm_);
			set_segm(e0, 1, e3, segm_);
			break;
		default:
			break;
		}
	}

	// connect oriented segments into oriented contours
	// closed contours are coded in 64 bit unsigned long long
	// 1) Each entry has 4 bits
	// 2) The first 4 entries are reserved for the size of the contours
	// 3) The next 12 entries are the indices of the edges constituting the contorus
	//    The indices are numbers from 0 to 12
	unsigned long long c_ = 0xFFFFFFFFFFFF0000;
	// in the 4 first bits store size of contours
	auto get_cnt_size = [](const int cnt, unsigned long long &c_) {
		return (size_t)((c_ & (0xF << 4 * cnt)) >> 4 * cnt);
	};
	auto set_cnt_size = [](const int cnt, const int size, unsigned long long &c_) {
		// unset contour size
		c_ &= ~(0xF << 4 * cnt);
		c_ |= (size << 4 * cnt);
	};
	// set corresponging edge
	auto set_c = [](const int cnt, const int pos, const int val, unsigned long long &c_) {
		const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
		const uint c_sz = c_ & mask[cnt];
		const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
		c_ &= ~(((unsigned long long)0xF) << e);
		c_ |= (((unsigned long long)val) << e);
	};
	// read edge from contour
	auto get_c = [](const int cnt, const int pos, unsigned long long c_) {
		const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
		const uint c_sz = (uint)(c_ & mask[cnt]);
		const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
		return (int)((c_ >> e) & 0xF);
	};


	// connect oriented contours
	uint cnt_{ 0 };
	for (uint e = 0; e < 12; e++) {
		if (is_segm_set(e, segm_)) {
			uint eTo = get_segm(e, 0, segm_);
			uint eIn = get_segm(e, 1, segm_);
			uint eStart = e;
			uint pos = 0;
			set_c(cnt_, pos, eStart, c_);
			while (eTo != eStart) {
				pos = pos + 1;
				set_c(cnt_, pos, eTo, c_);
				eIn = eTo;
				eTo = get_segm(eIn, 0, segm_);
				unset_segm(eIn, segm_);
			}
			// set contour length
			set_cnt_size(cnt_, pos + 1, c_);
			// update number of contours
			cnt_ = cnt_ + 1;
		}
	}

	// compute intersection of opposite faces
	// It is enough to compute a pair of solutions for one face
	// The other solutions are obtained by evaluating the equations
	// for the common variable
	double ui[2]{};
	double vi[2]{};
	double wi[2]{};
	unsigned char q_sol{ 0 };
	const double a = (F[0] - F[1])*(-F[6] + F[7] + F[4] - F[5]) - (F[4] - F[5])*(-F[2] + F[3] + F[0] - F[1]);
	const double b = (i0 - F[0])*(-F[6] + F[7] + F[4] - F[5]) + (F[0] - F[1])*(F[6] - F[4]) - (i0 - F[4])*(-F[2] + F[3] + F[0] - F[1]) - (F[4] - F[5])*(F[2] - F[0]);
	const double c = (i0 - F[0])*(F[6] - F[4]) - (i0 - F[4])*(F[2] - F[0]);;
	double d = b*b - 4 * a*c;
	if (d > 0) {
		d = std::sqrt(d);
		// compute u-coord of solutions
		ui[0] = (-b - d) / (2 * a);
		ui[1] = (-b + d) / (2 * a);
		// compute v-coord of solutions
		double g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
		double g2 = F[2] * (1 - ui[0]) + F[3] * ui[0];
		vi[0] = (i0 - g1) / (g2 - g1);
		if (std::isnan(vi[0]) || std::isinf(vi[0])) vi[0] = -1.f;
		g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
		g2 = F[2] * (1 - ui[1]) + F[3] * ui[1];
		vi[1] = (i0 - g1) / (g2 - g1);
		if (std::isnan(vi[1]) || std::isinf(vi[1])) vi[1] = -1.f;
		// compute w-coordinates of solutions
		g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
		g2 = F[4] * (1 - ui[0]) + F[5] * ui[0];
		wi[0] = (i0 - g1) / (g2 - g1);
		if (std::isnan(wi[0]) || std::isinf(wi[0])) wi[0] = -1.f;
		g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
		g2 = F[4] * (1 - ui[1]) + F[5] * ui[1];
		wi[1] = (i0 - g1) / (g2 - g1);
		if (std::isnan(wi[1]) || std::isinf(wi[1])) wi[1] = -1.f;
		// correct values for roots of quadratic equations
		// in case the asymptotic decider has failed
		if (f_flag[0] == true) { // face 1, w = 0;
			if (wi[0] < wi[1]) wi[0] = 0;
			else wi[1] = 0;
		}
		if (f_flag[1] == true) { // face 2, w = 1
			if (wi[0] > wi[1]) wi[1] = 1;
			else wi[1] = 1;
		}
		if (f_flag[2] == true) { // face 3, v = 0
			if (vi[0] < vi[1]) vi[0] = 0;
			else vi[1] = 0;
		}
		if (f_flag[3] == true) { // face 4, v = 1
			if (vi[0] > vi[1]) vi[0] = 1;
			else vi[1] = 1;
		}
		if (f_flag[4] == true) { // face 5, u = 0
			if (ui[0] < ui[1]) ui[0] = 0;
			else ui[1] = 0;
		}
		if (f_flag[5] == true) { // face 6, u = 1
			if (ui[0] > ui[1]) ui[0] = 1;
			else ui[1] = 1;
		}

		// check solution intervals
		if (0 < ui[0] && ui[0] < 1) {
			q_sol |= 1;
		}
		if (0 < ui[1] && ui[1] < 1) { 
			q_sol |= 2;
		}
		if (0 < vi[0] && vi[0] < 1) { 
			q_sol |= 4;
		}
		if (0 < vi[1] && vi[1] < 1) { 
			q_sol |= 8;
		}
		if (0 < wi[0] && wi[0] < 1) { 
			q_sol |= 16;
		} 
		if (0 < wi[1] && wi[1] < 1) { 
			q_sol |= 32;
		}
	}

	//
	// count the number of set bits
	auto numberOfSetBits = [](const unsigned char n) {
		// C or C++: use uint32_t
		uint b = (uint)n;
		b = b - ((b >> 1) & 0x55555555);
		b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
		return (((b + (b >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	};
	// compute the number of solutions to the quadratic equation for a given face
	auto nrQSolFace = [](const uint f, const unsigned char n)  {
		uint nr{ 0 };
		switch (f) {
		case 0:
			if ((n & 0x5) == 0x5)
				nr = nr + 1;
			if ((n & 0xA) == 0xA)
				nr = nr + 1;
			break;
		case 1:
			if ((n & 0x11) == 0x11) nr = nr + 1;
			if ((n & 0x22) == 0x22) nr = nr + 1;
			break;
		case 2:
			if ((n & 0x18) == 0x18) nr = nr + 1;
			if ((n & 0x24) == 0x24) nr = nr + 1;
			break;
		}
		return nr;
	};



	// triangulate contours
	// if all bits are set, then there are three pairs of nontrivial solutions
	// to the quadratic equations. In this case, there is a tunnel or a contour
	// with 12 vertices. If there are three contours, then there is a tunnel and
	// one of the contorus with only three vertices is not part of it.
	if (numberOfSetBits(q_sol) == 6) {
		m_ccases_tunnel += 1;
		// there are at most three contours
		// Possible cases:
		//  1) a single contour with 12 vertices
		//  2) two contours which build a tunnel
		//  3) three contours, one has only 3 vertices and does not belong to the tunnel

		// construct the six vertices of the inner hexagon
		double hvt[6][3];
		hvt[0][0] = ui[0]; hvt[0][1] = vi[0]; hvt[0][2] = wi[0];
		hvt[1][0] = ui[0]; hvt[1][1] = vi[0]; hvt[1][2] = wi[1];
		hvt[2][0] = ui[1]; hvt[2][1] = vi[0]; hvt[2][2] = wi[1];
		hvt[3][0] = ui[1]; hvt[3][1] = vi[1]; hvt[3][2] = wi[1];
		hvt[4][0] = ui[1]; hvt[4][1] = vi[1]; hvt[4][2] = wi[0];
		hvt[5][0] = ui[0]; hvt[5][1] = vi[1]; hvt[5][2] = wi[0];

		// construct vertices at intersections with the edges
		auto e_vert = [&ecoord](const int e, const int i) {
			const unsigned int l_coord[3]{ 1324855, 5299420, 16733440 };
			unsigned char flag = (l_coord[i] >> (2 * e)) & 3;
			if (flag == 3)
				return ecoord[e];
			else
				return (double)(flag);

		};

		// if there are three contours, then there is a tunnel and one
		// of the contours is not part of it.
		unsigned char _not_tunnel = 0xF;
		if (cnt_ == 3) {
			// loop over the contorus
			// triangulate the contour which is not part of
			// the tunnel
			const double uc_min = (ui[0] < ui[1]) ? ui[0] : ui[1];
			const double uc_max = (ui[0] < ui[1]) ? ui[1] : ui[0];
			for (int t = 0; t < (int)cnt_; t++) {
				if (get_cnt_size(t, c_) == 3) {
					double umin = 2;
					double umax = -2;
					uint e0 = get_c(t, 0, c_);
					uint e1 = get_c(t, 1, c_);
					uint e2 = get_c(t, 2, c_);
					const double u_e0 = e_vert(e0, 0);
					const double u_e1 = e_vert(e1, 0);
					const double u_e2 = e_vert(e2, 0);
					umin = (u_e0 < umin) ? u_e0 : umin;
					umin = (u_e1 < umin) ? u_e1 : umin;
					umin = (u_e2 < umin) ? u_e2 : umin;
					umax = (u_e0 > umax) ? u_e0 : umax;
					umax = (u_e1 > umax) ? u_e1 : umax;
					umax = (u_e2 > umax) ? u_e1 : umax;
					if (uc_min > umax || uc_max < umin) {
						// this contour is not part of the tunnel
						_not_tunnel = t;
						Triangle tr;
						tr.v[0] = vertices[e0].g_idx;
						tr.v[1] = vertices[e1].g_idx;
						tr.v[2] = vertices[e2].g_idx;
						m_triangles.push_back(tr);
						m_ccases_3++;
					}
				}
			}
		}

		// compute vertices of inner hexagon, save new vertices in list and compute and keep
		// global vertice index to build triangle connectivity later on.
		uint tg_idx[6];
		Point p_ih[6]; // remember the six points of inner hexagon for debugging
		for (int i = 0; i < 6; i++) {
			Point  hp;
			Normal hn;
			const double u = hvt[i][0]; const double v = hvt[i][1]; const double w = hvt[i][2];
			hp[0] = (1 - w)*((1 - v)*(p[0][0] + u*(p[1][0] - p[0][0])) + v*(p[2][0] + u*(p[3][0] - p[2][0]))) + w*((1 - v)*(p[4][0] + u*(p[5][0] - p[4][0])) + v*(p[6][0] + u*(p[7][0] - p[6][0])));
			hp[1] = (1 - w)*((1 - v)*(p[0][1] + u*(p[1][1] - p[0][1])) + v*(p[2][1] + u*(p[3][1] - p[2][1]))) + w*((1 - v)*(p[4][1] + u*(p[5][1] - p[4][1])) + v*(p[6][1] + u*(p[7][1] - p[6][1])));
			hp[2] = (1 - w)*((1 - v)*(p[0][2] + u*(p[1][2] - p[0][2])) + v*(p[2][2] + u*(p[3][2] - p[2][2]))) + w*((1 - v)*(p[4][2] + u*(p[5][2] - p[4][2])) + v*(p[6][2] + u*(p[7][2] - p[6][2])));
			hn[0] = (1 - w)*((1 - v)*(n[0][0] + u*(n[1][0] - n[0][0])) + v*(n[2][0] + u*(n[3][0] - n[2][0]))) + w*((1 - v)*(n[4][0] + u*(n[5][0] - n[4][0])) + v*(n[6][0] + u*(n[7][0] - n[6][0])));
			hn[1] = (1 - w)*((1 - v)*(n[0][1] + u*(n[1][1] - n[0][1])) + v*(n[2][1] + u*(n[3][1] - n[2][1]))) + w*((1 - v)*(n[4][1] + u*(n[5][1] - n[4][1])) + v*(n[6][1] + u*(n[7][1] - n[6][1])));
			hn[2] = (1 - w)*((1 - v)*(n[0][2] + u*(n[1][2] - n[0][2])) + v*(n[2][2] + u*(n[3][2] - n[2][2]))) + w*((1 - v)*(n[4][2] + u*(n[5][2] - n[4][2])) + v*(n[6][2] + u*(n[7][2] - n[6][2])));
			// normalize normal
			const double factor = std::sqrt(hn[0] * hn[0] + hn[1] * hn[1] + hn[2] * hn[2]);
			hn[0] = hn[0] / factor;
			hn[1] = hn[1] / factor;
			hn[2] = hn[2] / factor;
			tg_idx[i] = (uint)m_points.size();
			m_points.push_back(hp);
			m_pnorms.push_back(hn);
			p_ih[i] = hp;
		}

		// triangulate contours with inner hexagon
		std::vector<Triangle> lt_; // remember triangles for debugging
		unsigned char tcon_[12];
		for (int i = 0; i < (int)cnt_; i++) {
			if (_not_tunnel != i) { // contour belongs to tunnel
				const int cnt_sz = (int)get_cnt_size(i, c_);
				for (int r = 0; r < cnt_sz; r++) {
					uint index = -1;
					double dist = 1000.;
					uint ci = get_c(i, r, c_);
					const double u_edge = e_vert(ci, 0);
					const double v_edge = e_vert(ci, 1);
					const double w_edge = e_vert(ci, 2);
					for (int s = 0; s < 6; s++) {
						const double uval = u_edge - hvt[s][0];
						const double vval = v_edge - hvt[s][1];
						const double wval = w_edge - hvt[s][2];
						double val = uval*uval + vval*vval + wval*wval;
						if (dist > val) {
							index = s;
							dist = val;
						}
					}
					tcon_[ci] = (unsigned char)index;
				}

				// correspondence between vertices found
				// create triangles
				// needs some functions
				auto distanceRingIntsModulo = [](const int d1, const int d2) {
					const int r = (d1 - d2) < 0 ? d2 - d1 : d1 - d2;
					return (r > 2 ? 6 - r : r);
				};
				auto midpointRingIntModulo = [](const int d1, const int d2) {
					const int dmax = (d1 > d2) ? d1 : d2;
					const int dmin = (d1 < d2) ? d1 : d2;
					return ((dmax + 2) % 6 == dmin) ? (dmax + 1) % 6 : (dmax + dmin) / 2;
				};

				for (int r = 0; r < cnt_sz; r++) {
					const uint tid1 = get_c(i, r, c_);
					const uint tid2 = get_c(i, ((r + 1) % cnt_sz), c_);
					const uint cid1 = tcon_[tid1];
					const uint cid2 = tcon_[tid2];
					// compute index distance
					const int dst = distanceRingIntsModulo(cid1, cid2);
					switch (dst)
					{
					case 0:
					{
						Triangle tr;
						tr.v[0] = vertices[tid1].g_idx;
						tr.v[1] = vertices[tid2].g_idx;
						tr.v[2] = tg_idx[cid1];
						m_triangles.push_back(tr);
						lt_.push_back(tr);
					}
					break;
					case 1:
					{
						// measure diagonals
						// triangulate along shortest diagonal
						double u_edge = e_vert(tid1, 0);
						double v_edge = e_vert(tid1, 1);
						double w_edge = e_vert(tid1, 2);
						const double l1 = (u_edge - hvt[cid2][0])*(u_edge - hvt[cid2][0]) + (v_edge - hvt[cid2][1])*(v_edge - hvt[cid2][1]) + (w_edge - hvt[cid2][2])*(w_edge - hvt[cid2][2]);
						u_edge = e_vert(tid2, 0);
						v_edge = e_vert(tid2, 1);
						w_edge = e_vert(tid2, 2);
						const double l2 = (u_edge - hvt[cid1][0])*(u_edge - hvt[cid1][0]) + (v_edge - hvt[cid1][1])*(v_edge - hvt[cid1][1]) + (w_edge - hvt[cid1][2])*(w_edge - hvt[cid1][2]);
						Triangle t1;
						Triangle t2;
						if (l1 < l2) {
							t1.v[0] = vertices[tid1].g_idx;
							t1.v[1] = vertices[tid2].g_idx;
							t1.v[2] = tg_idx[cid2];
							t2.v[0] = vertices[tid1].g_idx;
							t2.v[1] = tg_idx[cid2];
							t2.v[2] = tg_idx[cid1];
						}
						else {
							t1.v[0] = vertices[tid1].g_idx;
							t1.v[1] = vertices[tid2].g_idx;
							t1.v[2] = tg_idx[cid1];
							t2.v[0] = vertices[tid2].g_idx;
							t2.v[1] = tg_idx[cid2];
							t2.v[2] = tg_idx[cid1];
						}
						m_triangles.push_back(t1);
						m_triangles.push_back(t2);
						lt_.push_back(t1);
						lt_.push_back(t2);
					}
					break;
					case 2:
					{
						const int cidm = midpointRingIntModulo(cid1, cid2);
						Triangle t1;
						Triangle t2;
						Triangle t3;
						t1.v[0] = vertices[tid1].g_idx;
						t1.v[1] = vertices[tid2].g_idx;
						t1.v[2] = tg_idx[cidm];

						t2.v[0] = vertices[tid1].g_idx;
						t2.v[1] = tg_idx[cidm];
						t2.v[2] = tg_idx[cid1];

						t3.v[0] = vertices[tid2].g_idx;
						t3.v[1] = tg_idx[cid2];
						t3.v[2] = tg_idx[cidm];

						m_triangles.push_back(t1);
						m_triangles.push_back(t2);
						m_triangles.push_back(t3);
						lt_.push_back(t1);
						lt_.push_back(t2);
						lt_.push_back(t3);
					}
					break;
					} // switch
				} // for loop over the vertices of the contour
			} // if (_not_tunnel)
		} // for loop over contours
		if (cnt_ == 1) {
			m_ccases_tunnel = m_ccases_tunnel - 1;
			m_ccases_12cont++;
			// there is a single contour
			// triangulate and close inner hexagon
			// triangle must have the correct orientation
			// use asymptotic_decider() to see if positive vertices
			// are separated, in thic case orientation must be changed
			const bool s_ = (asymptotic_decider(F[0], F[1], F[2], F[3]) <= i0);
			const bool of_ = (wi[1] < wi[0]) ? s_ : !s_;
			Triangle t1;
			Triangle t2;
			Triangle t3;
			Triangle t4;
			if (!of_) {
				t1.v[0] = tg_idx[0]; t1.v[1] = tg_idx[2]; t1.v[2] = tg_idx[1];
				t2.v[0] = tg_idx[2]; t2.v[1] = tg_idx[4]; t2.v[2] = tg_idx[3];
				t3.v[0] = tg_idx[0]; t3.v[1] = tg_idx[5]; t3.v[2] = tg_idx[4];
				t4.v[0] = tg_idx[0]; t4.v[1] = tg_idx[4]; t4.v[2] = tg_idx[2];
			}
			else {
				t1.v[0] = tg_idx[0]; t1.v[1] = tg_idx[1]; t1.v[2] = tg_idx[2];
				t2.v[0] = tg_idx[2]; t2.v[1] = tg_idx[3]; t2.v[2] = tg_idx[4];
				t3.v[0] = tg_idx[0]; t3.v[1] = tg_idx[4]; t3.v[2] = tg_idx[5];
				t4.v[0] = tg_idx[0]; t4.v[1] = tg_idx[2]; t4.v[2] = tg_idx[4];
			}
			m_triangles.push_back(t1);
			m_triangles.push_back(t2);
			m_triangles.push_back(t3);
			m_triangles.push_back(t4);
		}
	}
	else {
		// there is no tunnel
		// handle case with no saddle point as simple polygons with 3, 4, 5 or six vertices
		const unsigned char nr_u{ (unsigned char)nrQSolFace(0, q_sol) };
		const unsigned char nr_v{ (unsigned char)nrQSolFace(1, q_sol) };
		const unsigned char nr_w{ (unsigned char)nrQSolFace(2, q_sol) };
		const unsigned char nr_t{ (unsigned char)(nr_u + nr_v + nr_w) };
		if (nr_t == nr_u || nr_t == nr_v || nr_t == nr_w) {
			// loop over all contours
			for (int i = 0; i < (int)cnt_; i++) {
				switch (get_cnt_size(i, c_)) {
				case 3:
				{
					Triangle t1;
					t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
					t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
					m_triangles.push_back(t1);
					m_ccases_3++;
				}
				break;
				case 4:
				{
					Triangle t1;
					Triangle t2;
					t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
					t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
					t2.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t2.v[1] = vertices[get_c(i, 2, c_)].g_idx;
					t2.v[2] = vertices[get_c(i, 3, c_)].g_idx;
					m_triangles.push_back(t1);
					m_triangles.push_back(t2);
					m_ccases_4++;
				}
				break;
				case 5:
				{
					Triangle t1;
					Triangle t2;
					Triangle t3;
					t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
					t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
					t2.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t2.v[1] = vertices[get_c(i, 2, c_)].g_idx;
					t2.v[2] = vertices[get_c(i, 3, c_)].g_idx;
					t3.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t3.v[1] = vertices[get_c(i, 3, c_)].g_idx;
					t3.v[2] = vertices[get_c(i, 4, c_)].g_idx;
					m_triangles.push_back(t1);
					m_triangles.push_back(t2);
					m_triangles.push_back(t3);
					m_ccases_5++;
				}
				break;
				case 6:
				{
					Triangle t1;
					Triangle t2;
					Triangle t3;
					Triangle t4;
					t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
					t1.v[2] = vertices[get_c(i, 3, c_)].g_idx;
					t2.v[0] = vertices[get_c(i, 1, c_)].g_idx;
					t2.v[1] = vertices[get_c(i, 2, c_)].g_idx;
					t2.v[2] = vertices[get_c(i, 3, c_)].g_idx;
					t3.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t3.v[1] = vertices[get_c(i, 3, c_)].g_idx;
					t3.v[2] = vertices[get_c(i, 4, c_)].g_idx;
					t4.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t4.v[1] = vertices[get_c(i, 4, c_)].g_idx;
					t4.v[2] = vertices[get_c(i, 5, c_)].g_idx;
					m_triangles.push_back(t1);
					m_triangles.push_back(t2);
					m_triangles.push_back(t3);
					m_triangles.push_back(t4);
					m_ccases_6++;
				}
				break;
				} // switch over size of contour
			} // loop over contorus
		} // thre are no saddle points
		else {
			// there are saddle points
			//fc1 = fs(1, 1)*fs(2, 1) + fs(1, 2)*fs(2, 2);
			//fc2 = fs(1, 1)*fs(3, 1) + fs(1, 2)*fs(3, 2);
			//fc3 = fs(2, 1)*fs(3, 2) + fs(2, 2)*fs(3, 1);
			unsigned char fs[3][2]{{(uchar)(q_sol & 1), (uchar)((q_sol >> 1) & 1)}, { (uchar)((q_sol >> 2) & 1), (uchar)((q_sol >> 3) & 1) }, { (uchar)((q_sol >> 4) & 1), (uchar)((q_sol >> 5) & 1) }};

			const unsigned char fc1 = fs[0][0] * fs[1][0] + fs[0][1] * fs[1][1];
			const unsigned char fc2 = fs[0][0] * fs[2][0] + fs[0][1] * fs[2][1];
			const unsigned char fc3 = fs[1][0] * fs[2][1] + fs[1][1] * fs[2][0];
			const unsigned char c_faces = fc1 + fc2 + fc3;
			double ucoord{};
			double vcoord{};
			double wcoord{};
			switch (c_faces) {
			case 2:
			{
				if (fc1 == 0) {
					ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
					vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
					wcoord = fs[1][0] * wi[1] + fs[1][1] * wi[0];
				}
				else if (fc2 == 0) {
					ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
					vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
					wcoord = fs[0][0] * wi[1] + fs[0][1] * wi[0];
				}
				else if (fc3 == 0) {
					ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
					vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
					wcoord = fs[1][0] * wi[0] + fs[1][1] * wi[1];
				}
			}
			break;
			case 3:
			{
				ucoord = (fs[0][0] * ui[0] + fs[0][1] * ui[1]) / (fs[0][0] + fs[0][1]);
				vcoord = (fs[1][0] * vi[0] + fs[1][1] * vi[1]) / (fs[1][0] + fs[1][1]);
				wcoord = (fs[2][0] * wi[0] + fs[2][1] * wi[1]) / (fs[2][0] + fs[2][1]);
			}
			break;
			case 4:
			{
				const unsigned char nr_u = fs[0][0] + fs[0][1];
				const unsigned char nr_v = fs[1][0] + fs[1][1];
				const unsigned char nr_w = fs[2][0] + fs[2][1];
				if (nr_w == 1) {
					ucoord = fs[2][0] * ui[0] + fs[2][1] * ui[1];
					vcoord = fs[2][1] * vi[0] + fs[2][0] * vi[1];
					wcoord = fs[2][0] * wi[0] + fs[2][1] * wi[1];
				}
				else if (nr_v == 1) {
					ucoord = fs[1][0] * ui[0] + fs[1][1] * ui[1];
					vcoord = fs[1][0] * vi[0] + fs[1][1] * vi[1];
					wcoord = fs[1][1] * wi[0] + fs[1][0] * wi[1];
				}
				else if (nr_u == 1) {
					ucoord = fs[0][0] * ui[0] + fs[0][1] * ui[1];
					vcoord = fs[0][0] * vi[0] + fs[0][1] * vi[1];
					wcoord = fs[0][0] * wi[0] + fs[0][1] * wi[1];
				}
			}
			break;
			} // switch(c_faces)

			// create inner vertex
			Point  ip;
			Normal in;
			ip[0] = (1 - wcoord)*((1 - vcoord)*(p[0][0] + ucoord*(p[1][0] - p[0][0])) + vcoord*(p[2][0] + ucoord*(p[3][0] - p[2][0]))) + wcoord*((1 - vcoord)*(p[4][0] + ucoord*(p[5][0] - p[4][0])) + vcoord*(p[6][0] + ucoord*(p[7][0] - p[6][0])));
			ip[1] = (1 - wcoord)*((1 - vcoord)*(p[0][1] + ucoord*(p[1][1] - p[0][1])) + vcoord*(p[2][1] + ucoord*(p[3][1] - p[2][1]))) + wcoord*((1 - vcoord)*(p[4][1] + ucoord*(p[5][1] - p[4][1])) + vcoord*(p[6][1] + ucoord*(p[7][1] - p[6][1])));
			ip[2] = (1 - wcoord)*((1 - vcoord)*(p[0][2] + ucoord*(p[1][2] - p[0][2])) + vcoord*(p[2][2] + ucoord*(p[3][2] - p[2][2]))) + wcoord*((1 - vcoord)*(p[4][2] + ucoord*(p[5][2] - p[4][2])) + vcoord*(p[6][2] + ucoord*(p[7][2] - p[6][2])));
			in[0] = (1 - wcoord)*((1 - vcoord)*(n[0][0] + ucoord*(n[1][0] - n[0][0])) + vcoord*(n[2][0] + ucoord*(n[3][0] - n[2][0]))) + wcoord*((1 - vcoord)*(n[4][0] + ucoord*(n[5][0] - n[4][0])) + vcoord*(n[6][0] + ucoord*(n[7][0] - n[6][0])));
			in[1] = (1 - wcoord)*((1 - vcoord)*(n[0][1] + ucoord*(n[1][1] - n[0][1])) + vcoord*(n[2][1] + ucoord*(n[3][1] - n[2][1]))) + wcoord*((1 - vcoord)*(n[4][1] + ucoord*(n[5][1] - n[4][1])) + vcoord*(n[6][1] + ucoord*(n[7][1] - n[6][1])));
			in[2] = (1 - wcoord)*((1 - vcoord)*(n[0][2] + ucoord*(n[1][2] - n[0][2])) + vcoord*(n[2][2] + ucoord*(n[3][2] - n[2][2]))) + wcoord*((1 - vcoord)*(n[4][2] + ucoord*(n[5][2] - n[4][2])) + vcoord*(n[6][2] + ucoord*(n[7][2] - n[6][2])));
			// normalize normal
			const double factor = std::sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
			in[0] = in[0] / factor;
			in[1] = in[1] / factor;
			in[2] = in[2] / factor;
			const uint g_index = (uint)m_points.size();
			// loop over the contorus
			bool pt_used = false;
			for (int i = 0; i < (int)cnt_; i++) {
				const unsigned char cnt_sz = (unsigned char)get_cnt_size(i, c_);
				if (cnt_sz == 3) {
					Triangle t1;
					t1.v[0] = vertices[get_c(i, 0, c_)].g_idx;
					t1.v[1] = vertices[get_c(i, 1, c_)].g_idx;
					t1.v[2] = vertices[get_c(i, 2, c_)].g_idx;
					m_triangles.push_back(t1);
				}
				else {
					pt_used = true;
					for (int t = 0; t < cnt_sz; t++) {
						Triangle ts;
						ts.v[0] = vertices[get_c(i, t, c_)].g_idx;
						ts.v[1] = vertices[get_c(i, (t + 1) % cnt_sz, c_)].g_idx;
						ts.v[2] = g_index;
						m_triangles.push_back(ts);
					}
				}
				
				switch (cnt_sz) {
				case 3:
					m_ccases_3++;
					break;
				case 4:
					m_ccases_4++;
					break;
				case 5:
					m_ccases_5++;
					break;
				case 6:
					m_ccases_6a++;
					break;
				case 7:
					m_ccases_7++;
					break;
				case 8:
					m_ccases_8++;
					break;
				case 9:
					m_ccases_9++;
					break;
				default:
					break;
				}
			}
			if (pt_used) {
				m_points.push_back(ip);
				m_pnorms.push_back(in);
			}
		} // else - there are saddle points
	}

} // void p_slice()


//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
void
tmc::MarchingCubes::UGrid::init(const int nx,const int ny,const int nz)
{
    // set grid size
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    // compute spacing
    m_dx = 1. / (m_nx - 1);
    m_dy = 1. / (m_ny - 1);
    m_dz = 1. / (m_nz - 1);
    // total number of grid points
    size_t tot_size = m_nx*m_ny*m_nz;
    // initialize scalar fields
    m_scalars.resize(tot_size,0);
    m_normals.resize(tot_size);
    for (auto& p : m_normals) {
        p[0] = 0;
        p[1] = 0;
        p[2] = 0;
    }
    // create bounding box
    // fastest index is x, then y and slowest index is z
    for(int i = 0; i <= 1; i++) {
        for(int j = 0; j <= 1; j++) {
            for(int k = 0; k <= 1; k++) {
                int index = i*4 + j*2 + k;
				m_bbox[index] = Point{ { double(k), double(j), double(i) } };
            }
        }
    }
}

//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
void
tmc::MarchingCubes::UGrid::init(const int nx,const int ny,const int nz,std::array<Point, 8>& bb)
{
    // set grid size
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    // total number of grid points
    size_t tot_size = m_nx*m_ny*m_nz;
	m_scalars.clear();
	m_scalars.resize(tot_size,0);
	m_normals.clear();
	m_normals.resize(tot_size);
	for (auto& n : m_normals) {
		n[0] = 0;
		n[1] = 0;
		n[2] = 0;
	}

    // create bounding box
    // fastest index is x, then y and slowest index is z
    for(int i = 0; i <= 1; i++) {
        for(int j = 0; j <= 1; j++) {
            for(int k = 0; k <= 1; k++) {
                int index = i*4 + j*2 + k;
                m_bbox[index] = bb[index];
            }
        }
    }
    // compute spacing
    double x_space = m_bbox[7][0] - m_bbox[0][0];
    double y_space = m_bbox[7][1] - m_bbox[0][1];
    double z_space = m_bbox[7][2] - m_bbox[0][2];
    m_dx = x_space / (m_nx - 1);
    m_dy = y_space / (m_ny - 1);
    m_dz = z_space / (m_nz - 1);
}


void
tmc::MarchingCubes::UGrid::set(const int nx,const int ny,const int nz,float* tmp_a)
{
    // set grid size
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    // compute spacing
    m_dx = 1. / (m_nx - 1);
    m_dy = 1. / (m_ny - 1);
    m_dz = 1. / (m_nz - 1);
    // total number of grid points
    size_t tot_size = m_nx*m_ny*m_nz;
    // initialize scalar fields
    m_scalars.resize(tot_size,0);
    m_normals.resize(tot_size);

    // create bounding box
    // fastest index is x, then y and slowest index is z
    for(int i = 0; i <= 1; i++) {
        for(int j = 0; j <= 1; j++) {
            for(int k = 0; k <= 1; k++) {
                int index = i*4 + j*2 + k;
                m_bbox[index] = Point{double(k),double(j),double(i)};
            }
        }
    }

    for (int i = 0; i < tot_size; i++) {
        m_scalars[i] = (double) tmp_a[i];
    }
    gradient();
}


//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
// compute the divergence of the vector field defined on the grid
void
tmc::MarchingCubes::UniformGrid::divergence() {
    // compute inner nodes
    for (int k = 1; k < m_nz-1; k++) {
        for (int j = 1; j < m_ny-1; j++) {
            for (int i = 1; i < m_nx-1; i++) {
                // dfx/dx
                int idf = global_index(i+1, j, k);
                int idb = global_index(i-1, j, k);
                double fx = (m_normals[idf][0] - m_normals[idb][0]) / (2*m_dx);
                // dfy/dy
                idf = global_index(i, j+1, k);
                idb = global_index(i, j-1, k);
                double fy = (m_normals[idf][1] - m_normals[idb][1]) / (2*m_dy);
                // dfz/dz
                idf = global_index(i, j, k+1);
                idb = global_index(i, j, k-1);
                double fz = (m_normals[idf][2] - m_normals[idb][2]) / (2*m_dz);
                const int index = global_index(i,j,k);
                m_scalars[index]= fx+fy+fz;
            }
        }
    }
}


//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
// compute the gradient of the scalar field at all nodes of the ugrid
// do not consider boundary nodes, they are not needed for post processing
void
tmc::MarchingCubes::UniformGrid::gradient()
{
    // compute inner nodes
    double fx{0};
    double fy{0};
    double fz{0};
    for (int k = 0; k < m_nz; k++) {
        for (int j = 0; j < m_ny; j++) {
            for (int i = 0; i < m_nx; i++) {
                // dfx/dx
                if (i == 0) {
                    int idf = global_index(i+1, j, k);
                    int idb = global_index(i, j, k);
                    fx = (m_scalars[idf] - m_scalars[idb]) / m_dx;
                } else if (i == m_nx-1) {
                    int idf = global_index(i, j, k);
                    int idb = global_index(i-1, j, k);
                    fx = (m_scalars[idf] - m_scalars[idb]) / m_dx;
                } else {
                    int idf = global_index(i+1, j, k);
                    int idb = global_index(i-1, j, k);
                    fx = (m_scalars[idf] - m_scalars[idb]) / (2*m_dx);
                }
                // dfy/dy
                if (j == 0) {
                    int idf = global_index(i, j+1, k);
                    int idb = global_index(i, j, k);
                    fy = (m_scalars[idf] - m_scalars[idb]) / m_dy;
                } else if (j == m_ny-1) {
                    int idf = global_index(i, j, k);
                    int idb = global_index(i, j-1, k);
                    fy = (m_scalars[idf] - m_scalars[idb]) / m_dy;
                } else {
                    int idf = global_index(i, j+1, k);
                    int idb = global_index(i, j-1, k);
                    fy = (m_scalars[idf] - m_scalars[idb]) / (2*m_dy);
                }
                // dfz/dz
                if (k == 0) {
                    int idf = global_index(i, j, k+1);
                    int idb = global_index(i, j, k);
                    fz = (m_scalars[idf] - m_scalars[idb]) / m_dz;
                } else if (k == m_nz-1) {
                    int idf = global_index(i, j, k);
                    int idb = global_index(i, j, k-1);
                    fz = (m_scalars[idf] - m_scalars[idb]) / m_dz;
                } else {
                    int idf = global_index(i, j, k+1);
                    int idb = global_index(i, j, k-1);
                    fz = (m_scalars[idf] - m_scalars[idb]) / (2*m_dz);
                }
                const int index = global_index(i,j,k);
                m_normals[index][0] = fx;
                m_normals[index][1] = fy;
                m_normals[index][2] = fz;
            }
        }
    }
}


//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
double
tmc::MarchingCubes::UniformGrid::interpolate_scalar(const Point& p)
{
    double val{0};
    if (in_bbox(p)) {
        Index idx = cell_index(p);
        const int i = idx[0];
        const int j = idx[1];
        const int k = idx[2];
        Point vp = point(i,j,k);
        const double u = (p[0]-vp[0]) / m_dx;
        const double v = (p[1]-vp[1]) / m_dy;
        const double w = (p[2]-vp[2]) / m_dz;
        int gindex = global_index(i,j,k);
        const double s000 = m_scalars[gindex];
        gindex = global_index(i+1,j,k);
        const double s001 = m_scalars[gindex];
        gindex = global_index(i,j+1,k);
        const double s010 = m_scalars[gindex];
        gindex = global_index(i+1,j+1,k);
        const double s011 = m_scalars[gindex];
        gindex = global_index(i,j,k+1);
        const double s100 = m_scalars[gindex];
        gindex = global_index(i+1,j,k+1);
        const double s101 = m_scalars[gindex];
        gindex = global_index(i,j+1,k+1);
        const double s110 = m_scalars[gindex];
        gindex = global_index(i+1,j+1,k+1);
        const double s111 = m_scalars[gindex];

        val = (1-w)*((1-v)*((1-u)*s000 + u*s001) + v*((1-u)*s010 + u*s011))
        + w*((1-v)*((1-u)*s100 + u*s101) + v*((1-u)*s110 + u*s111));
    }

    return val;
}



//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
// interpolate the vector data at input position
bool
tmc::MarchingCubes::UniformGrid::interpolate_normal(const Point& p, Normal& n)
{
    // check if point is in ugrid
    if (in_bbox(p)) {
        Index id = cell_index(p);
        // collect indices
        const int i = id[0];
        const int j = id[1];
        const int k = id[2];
        Point vp = point(i,j,k);
        const double u = (p[0]-vp[0]) / m_dx;
        const double v = (p[1]-vp[1]) / m_dy;
        const double w = (p[2]-vp[2]) / m_dz;
        int gindex = global_index(i,j,k);
        Normal n000 = m_normals[gindex];
        gindex = global_index(mod(i+1,m_nx),j,k);
        Normal n001 = m_normals[gindex];
        gindex = global_index(i,mod(j+1,m_ny),k);
        Normal n010 = m_normals[gindex];
        gindex = global_index(mod(i+1,m_nx),mod(j+1,m_ny),k);
        Normal n011 = m_normals[gindex];
        gindex = global_index(i,j,mod(k+1,m_nz));
        Normal n100 = m_normals[gindex];
        gindex = global_index(mod(i+1,m_nx),j,mod(k+1,m_nz));
        Normal n101 = m_normals[gindex];
        gindex = global_index(i,mod(j+1,m_ny),mod(k+1,m_nz));
        Normal n110 = m_normals[gindex];
        gindex = global_index(mod(i+1,m_nx),mod(j+1,m_ny),mod(k+1,m_nz));
        Normal n111 = m_normals[gindex];

        n[0] = (1-w)*((1-v)*((1-u)*n000[0]+ u*n001[0]) + v*((1-u)*n010[0] + u*n011[0]))
             + w*((1-v)*((1-u)*n100[0] + u*n101[0]) + v*((1-u)*n110[0] + u*n111[0]));
        n[1] = (1-w)*((1-v)*((1-u)*n000[1]+ u*n001[1]) + v*((1-u)*n010[1] + u*n011[1]))
             + w*((1-v)*((1-u)*n100[1] + u*n101[1]) + v*((1-u)*n110[1] + u*n111[1]));
        n[2] = (1-w)*((1-v)*((1-u)*n000[2]+ u*n001[2]) + v*((1-u)*n010[2] + u*n011[2]))
             + w*((1-v)*((1-u)*n100[2] + u*n101[2]) + v*((1-u)*n110[2] + u*n111[2]));
        return true;
    } else {
		n = Normal{ { 0, 0, 0 } };
        return false;
    }
}


//*******************************************************************************************************************************************
//  IMPLEMENTATION t_mc
//*******************************************************************************************************************************************
void
tmc::MarchingCubes::s_mc(const double i0, UGrid& ugrid, std::vector<Point>& v_list, std::vector<Normal>& n_list, std::vector<Triangle>& t_list)
{
	// edges are uniquely characterized by the two end vertices, which have a unique vertex id
	// the end vertices of the edge are computed in the cell by giving the indices (i,j,k).
	// These indices are obtained from the cell index by adding 0 or 1 to i, j or k respectively
	// Example: edge 0: (i,j,k) - (i+1,j,k)
	//          edge 1: (i+1,j,k) - (i+1,j+1,k)
	// The first 3 indices are for the first vertex and the second 3 for the second vertex.
	// there are 12 edges, assign to each vertex three edges, the global edge numbering
	// consist of 3*global_vertex_id + edge_offset.
	const int global_edge_id[][4] = { { 0, 0, 0, 0 }, { 1, 0, 0, 1 }, { 0, 1, 0, 0 }, { 0, 0, 0, 1 },
	{ 0, 0, 1, 0 }, { 1, 0, 1, 1 }, { 0, 1, 1, 0 }, { 0, 0, 1, 1 },
	{ 0, 0, 0, 2 }, { 1, 0, 0, 2 }, { 1, 1, 0, 2 }, { 0, 1, 0, 2 } };
	// the end vertices of an edge
	int l_edges[12][2] = { { 0, 1 }, { 1, 3 }, { 2, 3 }, { 0, 2 },
	{ 4, 5 }, { 5, 7 }, { 6, 7 }, { 4, 6 },
	{ 0, 4 }, { 1, 5 }, { 3, 7 }, { 2, 6 } };
	// compute sizes
	const int nx = ugrid.x_size();
	const int ny = ugrid.y_size();
	const int nz = ugrid.z_size();

	// we need to compute up to 3 vertices at the interior of a cell, therefore
	// the cell shift factor is set to 3+3 = 6, i.e. 3 edges assigned to a cell for global numberig
	// and 3 vertices in the interior of the cell
	m_cell_shift_factor = 3;

	// there can be at most 12 intersections
	std::vector<Vertex> vertices(12);
	std::vector<Point>  ip(12);
	std::vector<Normal> in(12);
	// compute a unique global index for vertices
	// use as key the unique edge number
	std::map<int, int> v_map;

	timer.start();
	// marching cubes
	for (int k = 0; k < (nz - 1); k++) {
		m_kindex = k;
		for (int j = 0; j < (ny - 1); j++) {
			m_jindex = j;
			for (int i = 0; i < (nx - 1); i++) {
				m_iindex = i;
				// slice hex
				// collect function values and build index
				double u[8];
				Point  p[8];
				Normal n[8];
				int vi{ 0 };
				std::bitset<8> index = 0;
				for (int kl = 0; kl <= 1; kl++) {
					for (int jl = 0; jl <= 1; jl++) {
						for (int il = 0; il <= 1; il++) {
							// collect scalar values and computex index
							p[vi] = ugrid.point(i + il, j + jl, k + kl);
							u[vi] = ugrid.scalar(i + il, j + jl, k + kl);
							if (u[vi] >= i0) {
								//index.set(VertexMapping[vi]);
								index.set(vi);

							}
							// probably better get normals here
							n[vi] = ugrid.normal(i + il, j + jl, k + kl);
							// next cell vertex
							vi++;
						}
					}
				}

				// collect edges from table and
				// interpolate triangle vertex positon
				int i_case = int(index.to_ullong());
				// compute for this case the vertices
				ushort flag = 1;
				for (int eg = 0; eg < 12; eg++) {
					if (flag & e_pattern[i_case]) {
						// the edge global index is given by the vertex global index + the edge offset
						const int ix = i + global_edge_id[eg][0];
						const int iy = j + global_edge_id[eg][1];
						const int iz = k + global_edge_id[eg][2];
						vertices[eg].g_edg = uint(m_cell_shift_factor*ugrid.global_index(ix, iy, iz) + global_edge_id[eg][3]);
						// generate vertex here, do not care at this point if vertex already exist
						int* vert = l_edges[eg];
						// interpolation weight
						const int v0 = vert[0];
						const int v1 = vert[1];
						double l = (i0 - u[v0]) / (u[v1] - u[v0]);
						// interpolate vertex
						ip[eg][0] = (1 - l)*p[v0][0] + l*p[v1][0];
						ip[eg][1] = (1 - l)*p[v0][1] + l*p[v1][1];
						ip[eg][2] = (1 - l)*p[v0][2] + l*p[v1][2];

						// interpolate normal
						in[eg][0] = (1 - l)*n[v0][0] + l*n[v1][0];
						in[eg][1] = (1 - l)*n[v0][1] + l*n[v1][1];
						in[eg][2] = (1 - l)*n[v0][2] + l*n[v1][2];
						const double nlength = std::sqrt(in[eg][0] * in[eg][0] + in[eg][1] * in[eg][1] + in[eg][2] * in[eg][2]);
						in[eg][0] = in[eg][0] / nlength;
						in[eg][1] = in[eg][1] / nlength;
						in[eg][2] = in[eg][2] / nlength;

						// set vertex index
						auto s_index = v_map.find(vertices[eg].g_edg);
						if (s_index == v_map.end()) {
							// index not found! Add index to hash map
							const int g_idx = (int)v_list.size();
							v_map[vertices[eg].g_edg] = g_idx;
							vertices[eg].g_idx = g_idx;
							v_list.push_back(ip[eg]);
							n_list.push_back(in[eg]);
						}
						else {
							vertices[eg].g_idx = s_index->second; // this is vertex global index g_idx
						}
					}
					flag <<= 1;
				}

				// construct triangles
				for (int t = 0; t < 16; t += 3) {
					const int t_index = i_case * 16 + t;
					//if (e_tris_list[t_index] == 0x7f)
					if (t_pattern[t_index] == -1)
						break;
					Triangle tri;
					const int eg0 = t_pattern[t_index];
					const int eg1 = t_pattern[t_index + 1];
					const int eg2 = t_pattern[t_index + 2];
					tri.v[0] = (int)vertices[eg0].g_idx;
					tri.v[1] = (int)vertices[eg1].g_idx;
					tri.v[2] = (int)vertices[eg2].g_idx;

					// insert new triangle in list
					t_list.push_back(tri);
				}
			}
		}
	}

	timer.stop();
	std::cout << "Standard Marching Cubes: \n";
	timer.print();

	std::cout << "tot. nr. of vertices:        " << v_list.size() << std::endl;
	std::cout << "tot. nr. of triangles:       " << t_list.size() << std::endl;
}




//*******************************************************************************************************************************************
//  IMPLEMENTATION UniformGrid
//*******************************************************************************************************************************************
double
tmc::MarchingCubes::UGrid::trilinear(const double u, const double v, const double w, const std::array<double, 8>& f)
{
    double val = (1-w)*((1-v)*(f[0]*(1-u)+f[1]*u)+v*(f[2]*(1-u)+f[3]*u)) + w*((1-v)*(f[4]*(1-u)+f[5]*u)+v*(f[6]*(1-u)+f[7]*u));
    return val;
}



void tmc::MarchingCubes::connectivity()
{
	const int nr_v = (int)m_points.size();
	// check if a triangle is degenerated, i.e. it has to vertices with the same index
	std::cout << " ... checking triangle indices\n";
	auto c_index = [](Triangle& t) {
		if (t.v[0] == t.v[1])
			return false;
		if (t.v[0] == t.v[2])
			return false;
		if (t.v[1] == t.v[2])
			return false;
		return true;
	};

	// assign triangle to pos in triangle array
	for (int t = 0; t < (int)m_triangles.size(); t++) {
		m_triangles[t].id = t;
		if (!c_index(m_triangles[t])) {
			std::cout << "triangle " << t << ", has two vertices equal: " << m_triangles[t].v[0] << ", " << m_triangles[t].v[1] << ", " << m_triangles[t].v[2] << std::endl;
		}

	}
	// check vertices ids
	m_vboundary.resize(nr_v);
	for (int i = 0; i < nr_v; i++) {
		m_vboundary[i] = false;
	}

	// check triangle orientation
	// collect for each vertex its neighbor triangles
	std::vector<std::vector<int>> v_tris(nr_v);
	for (auto t : m_triangles) {
		v_tris[t.v[0]].push_back(t.id);
		v_tris[t.v[1]].push_back(t.id);
		v_tris[t.v[2]].push_back(t.id);
	}

	// check if there are triangles with the same three vertices
	// this might happens when triangulating cells in the case
	// singularities, e.g. asymptotic decider does not work 
	// properly because isovalue equals decider value
	auto compare_vertexIds = [ & ](const std::vector<int>& tris) {
		for (int i = 0; i < (int)tris.size(); i++) {
			const int t = tris[i];
			const int t0 = m_triangles[i].v[0];
			const int t1 = m_triangles[i].v[1];
			const int t2 = m_triangles[i].v[2];
			for (int j = (i + 1); j < (int) tris.size(); j++) {
				const int n = tris[j];
				const int n0 = m_triangles[n].v[0];
				const int n1 = m_triangles[n].v[1];
				const int n2 = m_triangles[n].v[2];
				const bool f0_ = (t0 == n0) || (t0 == n1) || (t0 == n2);
				const bool f1_ = (t1 == n0) || (t1 == n1) || (t1 == n2);
				const bool f2_ = (t2 == n0) || (t2 == n1) || (t2 == n2);
				if (f0_ && f1_ && f2_) {
					return true;
				}
			}
		}
		return false;
	};
	std::cout << " ... checking for triangles with the same vertex indices\n";
	for (int i = 0; i < nr_v; i++) {
		if (compare_vertexIds(v_tris[i])) {
			std::cout << "ERROR: there are at least two triangles with the same vertex indices\n";
			std::cout << "       check the triangles: \n";
			for (auto t : v_tris[i]) {
				std::cout << "       triangle id: " << t << std::endl;
			}
		}
	}



	// check vertex fan
	// collect adjacent vertices
	std::cout << " ... checking for triangle's neighborhood\n";
	for (int v = 0; v < nr_v; v++) {
		std::vector<int> v_fan;
		for (auto t : v_tris[v]) {
			int t0 = m_triangles[t].v[0];
			int t1 = m_triangles[t].v[1];
			int t2 = m_triangles[t].v[2];
			if (t0 != v) v_fan.push_back(t0);
			if (t1 != v) v_fan.push_back(t1);
			if (t2 != v) v_fan.push_back(t2);
		}
		// is this an isolated vertex?
		if (v_fan.size() == 0) {
			std::cerr << "ERROR: vertex " << v << " is isolated\n";
			continue;
		}
		// check triangle fan
		if (v_fan.size() % 2 != 0) {
			std::cout << "ERROR: wrong number of vertices in triangle fan\n";
		}
		std::sort(v_fan.begin(), v_fan.end());
		std::vector<int> v_neigh(v_fan.size());
		if (v_fan[0] == v_fan[1])
			v_neigh[0] = 2;
		else
			v_neigh[0] = 1;
		for (int n = 1; n < (v_fan.size() - 1); n++) {
			if (v_fan[n] == v_fan[n - 1] || v_fan[n] == v_fan[n + 1])
				v_neigh[n] = 2;
			else
				v_neigh[n] = 1;
		}
		if (v_fan[v_fan.size() - 1] == v_fan[v_fan.size() - 2])
			v_neigh[v_fan.size() - 1] = 2;
		else
			v_neigh[v_fan.size() - 1] = 1;

		// check nr. of vertices
		int v_one{ 0 };
		for (auto n : v_neigh) {
			if (n == 0) {
				std::cout << "ERROR: vertex v0: " << v << ", no neighbors ZERO\n";
			}
			else if (n == 1) {
				m_vboundary[v] = true;
				v_one = v_one + 1;
			}
		}
		if (v_one == 1 || v_one > 2) {
			std::cout << "ERROR: wrong number of neighbors\n";
		}
	}

	// check if edge (v1,v2) is in triangle (t0,t1,t2)
	auto v_compare = [](const int v1, const int v2, const int t0, const int t1, const int t2) {
		bool face0 = (v1 == t1 && v2 == t2) || (v1 == t2 && v2 == t1);
		bool face1 = (v1 == t0 && v2 == t2) || (v1 == t2 && v2 == t0);
		bool face2 = (v1 == t1 && v2 == t0) || (v1 == t0 && v2 == t1);

		if (face0 || face1 || face2) {
			return true;
		}
		else {
			return false;
		}
	};

	// create face neighbors
	// loop over all triangle and set neighors
	std::vector<std::array<int, 3>> t_neigh(m_triangles.size());
	for (auto t : m_triangles) {
        const int v0 = t.v[0];
        const int v1 = t.v[1];
        const int v2 = t.v[2];
        t_neigh[t.id][0] = -1;
        t_neigh[t.id][1] = -1;
        t_neigh[t.id][2] = -1;


        // face 0
        int t_topol{0};
        bool f_neigh{false};
        for (auto tt : v_tris[t.v[1]]) {
            if (tt == t.id)
                continue;
            int t0 = m_triangles[tt].v[0];
            int t1 = m_triangles[tt].v[1];
            int t2 = m_triangles[tt].v[2];
            if (v_compare(v1,v2,t0,t1,t2)) {
                t_neigh[t.id][0] = tt;
                f_neigh = true;
                break;
            }
        }
        if (!f_neigh) {
            // this face has no neighbor, vertices must be boundary vertices
            if (!m_vboundary[v1]) {
                std::cout << "ERROR: triangle " << t.id << ", has no neighbor at face 0 but vertex v1 is not boundary!\n";
            }
            if (!m_vboundary[v2]) {
                std::cout << "ERROR: triangle " << t.id << ", has not neighbor at face 0 but vertex v2 is not boundary!\n";
            }
		}
		else {
			t_topol += 1;
		}

        // face 1
		f_neigh = false;
        for (auto tt : v_tris[t.v[2]]) {
            if (tt == t.id)
                continue;
            int t0 = m_triangles[tt].v[0];
            int t1 = m_triangles[tt].v[1];
            int t2 = m_triangles[tt].v[2];
            if (v_compare(v2,v0,t0,t1,t2)) {
                t_neigh[t.id][1] = tt;
                f_neigh = true;
                break;
            }
        }
        if (!f_neigh) {
            // this face has no neighbor, vertices must be boundary vertices
            if (!m_vboundary[v0]) {
                std::cout << "ERROR: triangle " << t.id << ", has no neighbor at face 1 but vertex v0 is not boundary!\n";
            }
            if (!m_vboundary[v2]) {
                std::cout << "ERROR: triangle " << t.id << ", has not neighbor at face 1 but vertex v2 is not boundary!\n";
            }
        }
		else {
			t_topol += 1;
		}

        // face 2
		f_neigh = false;
        for (auto tt : v_tris[t.v[0]]) {
            if (tt == t.id)
                continue;
            int t0 = m_triangles[tt].v[0];
            int t1 = m_triangles[tt].v[1];
            int t2 = m_triangles[tt].v[2];
            if (v_compare(v0,v1,t0,t1,t2)) {
                t_neigh[t.id][2] = tt;
                f_neigh = true;
                break;
            }
        }
        if (!f_neigh) {
            // this face has no neighbor, vertices must be boundary vertices
            if (!m_vboundary[v0]) {
                std::cout << "ERROR: triangle " << t.id << ", has not neighbor at face 2 but vertex v0 is not boundary!\n";
            }
            if (!m_vboundary[v1]) {
                std::cout << "ERROR: triangle " << t.id << ", has no neighbor at face 2 but vertex v1 is not boundary!\n";
            }
        }
		else {
			t_topol += 1;
		}

        if (t_topol == 0) {
            std::cout << "ERROR: triangle " << t.id << " has no neighbors at all, it is an isolated triangle\n";
        }
    }

    // check neighborhood
	// at this points all neighbors were set
	// that means, the neighbors should also have the
	// triangle as neighbor
    for (auto t : m_triangles) {
        int v0 = t.v[0];
        int v1 = t.v[1];
        int v2 = t.v[2];
        int n0 = t_neigh[t.id][0];
        int n1 = t_neigh[t.id][1];
        int n2 = t_neigh[t.id][2];
        if (n0 != -1) {
            int t0 = m_triangles[n0].v[0];
            int t1 = m_triangles[n0].v[1];
            int t2 = m_triangles[n0].v[2];
            if (!v_compare(v1,v2,t0,t1,t2)) {
                std::cout << "ERROR: wrong neighbor for triangle " << t.id << ", neighbor 0 = " << n0 << std::endl;
            }
        }
        if (n1 != -1) {
            int t0 = m_triangles[n1].v[0];
            int t1 = m_triangles[n1].v[1];
            int t2 = m_triangles[n1].v[2];
            if (!v_compare(v2,v0,t0,t1,t2)) {
                std::cout << "ERROR: wrong neighbor for triangle " << t.id << ", neighbor 1 = " << n1 << std::endl;
            }
        }
        if (n2 != -1) {
            int t0 = m_triangles[n2].v[0];
            int t1 = m_triangles[n2].v[1];
            int t2 = m_triangles[n2].v[2];
            if (!v_compare(v0,v1,t0,t1,t2)) {
                std::cout << "ERROR: wrong neighbor for triangle " << t.id << ", neighbor 2 = " << n2 << std::endl;
            }
        }
    }


    // Check and compute triangle orientation
	std::cout << " ... checking triangle orientation\n";
    std::vector<bool> t_oriented(m_triangles.size(),false);
    std::stack<int> t_queue;
    int nr_ocorrections{0};
    auto t_compare = [&nr_ocorrections,&t_neigh,t_oriented] (const int v1, const int v2, Triangle& tri) {
        const int t0 = tri.v[0];
        const int t1 = tri.v[1];
        const int t2 = tri.v[2];
        bool face0 = (v1 == t2 && v2 == t1);
        bool face1 = (v1 == t0 && v2 == t2);
        bool face2 = (v1 == t1 && v2 == t0);

        if (face0 || face1 || face2) {
            return true;
        }
        else {
            // needs to change orientation
            if (t_oriented[tri.id]) {
                std::cout << "ERROR: should not reorient triangle\n";
            }
            nr_ocorrections = nr_ocorrections + 1;
            tri.v[1] = t2;
            tri.v[2] = t1;
            int trin1 = t_neigh[tri.id][1];
            int trin2 = t_neigh[tri.id][2];
            t_neigh[tri.id][1] = trin2;
            t_neigh[tri.id][2] = trin1;
            return true;
        }
        // something wrong
        std::cout << "ERROR: can't change or test triangle orientation\n";
        return false;
    };

	// loop over all triangles, set element as oriented, this is the
	// default orientation of the connected component in the triangle mesh.
	// if triangle was not processed, i.e. it is not oriented
	// collect neighbors into  stack, neighbors are correctly oriented
	// such that elements in the stack have the same orientation 
	// as the first element. Process the elements in the stack the same
	// way, collect neighbors, correctly orient neighbors if necessary
	// and push them into the stack, then mark element as oriented
	// such that it will not processed again.
	// When the stack is emtpy, a connected component of the
	// triangle mesh was processed, then go to the next component.
    int t_nrcomps{0};
    for (auto t : m_triangles) {
        if (t_oriented[t.id] == true) {
            continue;
        }
        t_nrcomps = t_nrcomps + 1;
        // triangle still not processed
        t_oriented[t.id] = true;
        int v0 = t.v[0];
        int v1 = t.v[1];
        int v2 = t.v[2];
		// collect neighbors sharing a face
        int n0 = t_neigh[t.id][0];
        int n1 = t_neigh[t.id][1];
        int n2 = t_neigh[t.id][2];
        // face 0
        if (n0 != -1) {
			if (t_oriented[n0]) {
				std::cout << "for triangle " << t.id << ", neighbor " << n0 << " is oriented" << std::endl;
			}
			else {
				t_compare(v1, v2, m_triangles[n0]);
				t_queue.push(n0);
			}
        }
        // face 1
        if (n1 != -1) {
			if (t_oriented[n1]) {
				std::cout << "for triangle " << t.id << ", neighbor " << n1 << " is oriented" << std::endl;
			} 
			else {
				t_compare(v2, v0, m_triangles[n1]);
				t_queue.push(n1);
			}
        }

        // face 2
        if (n2 != -1) {
			if (t_oriented[n2]) {
				std::cout << "for triangle " << t.id << ", neighbor " << n2 << " is oriented" << std::endl;
			}
			else {
				t_compare(v0, v1, m_triangles[n2]);
				t_queue.push(n2);
			}
        }
        // set stack
        while (!t_queue.empty()) {
            int c_tri = t_queue.top();
            t_queue.pop();
            if (!t_oriented[c_tri]) {
                // triangle still not processed
                t_oriented[c_tri] = true;
                Triangle tr_ = m_triangles[c_tri];
                int t0 = tr_.v[0];
                int t1 = tr_.v[1];
                int t2 = tr_.v[2];
                int tn0 = t_neigh[tr_.id][0];
                int tn1 = t_neigh[tr_.id][1];
                int tn2 = t_neigh[tr_.id][2];
                // face 0
                if (tn0 != -1 && !t_oriented[tn0]) {
                    t_compare(t1,t2,m_triangles[tn0]);
                    t_queue.push(tn0);
                }
                // face 1
                if (tn1 != -1 && !t_oriented[tn1]) {
                    t_compare(t2,t0,m_triangles[tn1]);
                    t_queue.push(tn1);
                }
                // face 2
                if (tn2 != -1 && !t_oriented[tn2]) {
                    t_compare(t0,t1,m_triangles[tn2]);
                    t_queue.push(tn2);
                }
            } // if () // triangle still not processed
        } // while
    }

    // check orientation
    auto t_checko = [](const int v1,const int v2, const int t0,const int t1, const int t2) {
        bool face0 = (v1 == t2 && v2 == t1);
        bool face1 = (v1 == t0 && v2 == t2);
        bool face2 = (v1 == t1 && v2 == t0);

        if (face0 || face1 || face2) {
            return true;
        } else {
            return false;
        }
    };
    for (auto t : m_triangles) {
        const int v0 = t.v[0];
        const int v1 = t.v[1];
        const int v2 = t.v[2];
        const int n0 = t_neigh[t.id][0];
        const int n1 = t_neigh[t.id][1];
        const int n2 = t_neigh[t.id][2];
        if (n0 != -1) {
            const int t0 = m_triangles[n0].v[0];
            const int t1 = m_triangles[n0].v[1];
            const int t2 = m_triangles[n0].v[2];
            if (!t_checko(v1,v2,t0,t1,t2)) {
                std::cout << "ERROR: neighbor n0 not oriented for t " << t.id << std::endl;
            }
        }

        if (n1 != -1) {
            const int t0 = m_triangles[n1].v[0];
            const int t1 = m_triangles[n1].v[1];
            const int t2 = m_triangles[n1].v[2];
            if (!t_checko(v2,v0,t0,t1,t2)) {
                std::cout << "ERROR: neighbor n1 not oriented for t " << t.id << std::endl;
            }
        }

        if (n2 != -1) {
            const int t0 = m_triangles[n2].v[0];
            const int t1 = m_triangles[n2].v[1];
            const int t2 = m_triangles[n2].v[2];
            if (!t_checko(v0,v1,t0,t1,t2)) {
                std::cout << "ERROR: neighbor n2 not oriented for t " << t.id << std::endl;
            }
        }
    }

    std::cout << " connectivity has changed " << nr_ocorrections << " triangle orientations\n";
    std::cout << " connectivity found " << t_nrcomps << " components" << std::endl;

}
