//  MarchingCubes.cpp
//  Quitte
//
//  Created by Roberto Grosso on 05.06.16.
//  Copyright © 2016 Roberto Grosso. All rights reserved.
//

#include "MarchingCubes.h"


void
tmc::MarchingCubes::operator() (const std::string& i_file,const std::string& o_objF,const std::string& o_offF)
{
    std::cout << " ... reading data \n";

	/*std::FILE* f{ nullptr };
	errno_t status = fopen_s(&f, i_file.c_str(), "rb"); 
	if (status != 0) {
		std::cerr << "ERROR: can't open file " << i_file.c_str() << std::endl;
		exit(1);
	}*/
	std::FILE* f = fopen(i_file.c_str(), "rb");
    short x_size;
    short y_size;
    short z_size;
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
    //m_ugrid.gradient();

    // invert gradient
    //m_ugrid.invert_normals();
    m_ugrid.invert_normals();

    // compute isosurface
    std::cout << " ... computing isosurface\n";
	const double i0 = 700.01;
    t_mc(i0);


    const int nr_v = (int)m_points.size();
    const int nr_t = (int)m_triangles.size();

    std::cout << "tot. nr. of triangles: " << nr_t << std::endl;
    std::cout << "tot. nr. of vertices:  " << nr_v << std::endl;


    // Check mesh topology and correct triangle orientation if necessary
    std::cout << " ... check and correct triangle orientation\n";
    connectivity();


    // write obj
    std::cout << " ... write obj file\n";
    std::ofstream objF;
    objF.open(o_objF.c_str());
    objF << "# Topologically correct and manifold isosurface\n";
    for (int i = 0; i < nr_v; i++) {
        objF << "v " << std::setprecision(7) << std::fixed << m_points[i][0] << " " << m_points[i][1] << " " << m_points[i][2] << std::endl;
    }
    for (int n = 0; n < nr_v; n++) {
        objF << "vn " << std::setprecision(7) << std::fixed << m_pnorms[n][0] << " " << m_pnorms[n][1] << " " << m_pnorms[n][2] << std::endl;
    }
    for (auto t : m_triangles) {
        objF << "f " << (t.v[0]+1) << " " << (t.v[1]+1) << " " << (t.v[2]+1) << std::endl;
    }
    objF.close();

    // write obj output file
    std::cout << " ... write OFF file\n";
    std::ofstream offF;
    offF.open(o_offF.c_str());
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
                        const int v0 = (int)vertices[eg0].g_idx;
                        const int v1 = (int)vertices[eg1].g_idx;
                        const int v2 = (int)vertices[eg2].g_idx;
                        tri.v[0] = (int)vertices[eg0].g_idx;
                        tri.v[1] = (int)vertices[eg1].g_idx;
                        tri.v[2] = (int)vertices[eg2].g_idx;

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
//  IMPLEMENTATION
//*******************************************************************************************************************************************
void
tmc::MarchingCubes::t_slice(const int i,const int j,const int k,const double i0,double* F, Point* p,Normal* n,const int i_case)
{
    const int face_e[6][4] = {{0,1,2,3},{4,5,6,7},{0,9,4,8},{2,10,6,11},{3,11,7,8},{1,10,5,9}};
    //unsigned short face_e_[6] = {12816, 30292, 33936, 46754, 34739, 38305};
    const int face_v[6][4] = {{0,1,2,3},{4,5,6,7},{0,1,4,5},{2,3,6,7},{0,2,4,6},{1,3,5,7}};
    //unsigned short face_v_[6] = {12816, 30292, 21520, 30258, 25632, 30001};

    // there are 12 edges, assign to each vertex three edges, the global edge numbering
    // consist of 3*global_vertex_id + edge_offset.
    const int global_edge_id[][4] = {{0,0,0,0},{1,0,0,1},{0,1,0,0},{0,0,0,1},
        {0,0,1,0},{1,0,1,1},{0,1,1,0},{0,0,1,1},
        {0,0,0,2},{1,0,0,2},{1,1,0,2},{0,1,0,2}};
    // the end vertices of an edge
    int l_edges[12][2] = {{0,1}, {1,3}, {2,3}, {0,2},
        {4,5}, {5,7}, {6,7}, {4,6},
        {0,4}, {1,5}, {3,7}, {2,6}};
    // remember if an edge is intersected by the level set
    bool e_set[12];
    // there might be up to 15 vertices and normals including three vertices or normals at the interior
    std::vector<Vertex> vertices(12);
    std::vector<Point>  ip(12);
    std::vector<Normal> in(12);
    std::vector<double> ecoord(12); // there are 12 coordinates along the edges, these are the intersections



    // debugging
    //    const int nr_vr_0 = m_tvertices.size();
    //    const int nr_tr_0 = m_ttriangles.size();

    // collect vertices
    ushort   flag{1};
    for (int eg = 0; eg < 12; eg++) {
        if (flag & e_pattern[i_case]) {
            e_set[eg] = true;
            // the edge global index is given by the vertex global index + the edge offset
            const int ix = i + global_edge_id[eg][0];
            const int iy = j + global_edge_id[eg][1];
            const int iz = k + global_edge_id[eg][2];
            vertices[eg].g_edg = int(m_cell_shift_factor*m_ugrid.global_index(ix, iy, iz) + global_edge_id[eg][3]);
            // generate vertex here, do not care at this point if vertex already exist
            int* vert = l_edges[eg];
            // interpolation weight
            int v0 = vert[0];
            int v1 = vert[1];
            double l = (i0 - F[v0]) / (F[v1] - F[v0]);
            ecoord[eg] = l;
            // interpolate vertex
            // interpolate vertex
            ip[eg][0] = (1 - l)*p[v0][0] + l*p[v1][0];
            ip[eg][1] = (1 - l)*p[v0][1] + l*p[v1][1];
            ip[eg][2] = (1 - l)*p[v0][2] + l*p[v1][2];

            // interpolate normal
            in[eg][0] = (1 - l)*n[v0][0] + l*n[v1][0];
            in[eg][1] = (1 - l)*n[v0][1] + l*n[v1][1];
            in[eg][2] = (1 - l)*n[v0][2] + l*n[v1][2];
            const double nlength = std::sqrt(in[eg][0]*in[eg][0] + in[eg][1]*in[eg][1] + in[eg][2]*in[eg][2]);
            in[eg][0] = in[eg][0] / nlength;
            in[eg][1] = in[eg][1] / nlength;
            in[eg][2] = in[eg][2] / nlength;

            // set vertex index
            auto s_index = m_vertices.find(vertices[eg].g_edg);
            if (s_index == m_vertices.end()) {
                const int g_idx = (int)m_points.size();
                vertices[eg].g_idx = g_idx;
                m_vertices[vertices[eg].g_edg] = g_idx;
                m_points.push_back(ip[eg]);
                m_pnorms.push_back(in[eg]);
            } else {
                vertices[eg].g_idx = s_index->second;
            }
        } else {
            e_set[eg] = false;
        }
        //next edge
        flag <<=1;
    }

    // build up segments
    int segments[12][2] = {{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1},{-1,-1}};
    int spos[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    for (int f = 0; f < 6; f++) {
        int pos = 0;
        int s[4] = {0,0,0,0};
        for (int eg = 0; eg < 4; eg++) {
            int eid = face_e[f][eg];
            if (e_set[eid])
                s[pos++] = eid;
        }
        if (pos == 2) {
            segments[s[0]][spos[s[0]]] = s[1];
            segments[s[1]][spos[s[1]]] = s[0];
            spos[s[0]] += 1;
            spos[s[1]] += 1;
        } else if (pos == 4) {
            // build up contour segments for this face
            // need asymptotes
            const double f0 = F[face_v[f][0]];
            const double f1 = F[face_v[f][1]];
            const double f2 = F[face_v[f][2]];
            const double f3 = F[face_v[f][3]];
            const double eta = f0+f3-f1-f2;
            const double u0 = (f0-f2)/eta;
            //const double v0 = (f0-f1)/eta;
            //
            if (ecoord[face_e[f][0]] < u0) {
                segments[face_e[f][0]][spos[face_e[f][0]]]   = face_e[f][3];
                segments[face_e[f][3]][spos[face_e[f][3]]]   = face_e[f][0];
                segments[face_e[f][1]][spos[face_e[f][1]]]   = face_e[f][2];
                segments[face_e[f][2]][spos[face_e[f][2]]]   = face_e[f][1];

            } else {
                segments[face_e[f][0]][spos[face_e[f][0]]]   = face_e[f][1];
                segments[face_e[f][1]][spos[face_e[f][1]]]   = face_e[f][0];
                segments[face_e[f][2]][spos[face_e[f][2]]]   = face_e[f][3];
                segments[face_e[f][3]][spos[face_e[f][3]]]   = face_e[f][2];

            }
            spos[face_e[f][0]] += 1;
            spos[face_e[f][1]] += 1;
            spos[face_e[f][2]] += 1;
            spos[face_e[f][3]] += 1;
        }
    } //

    // there might give at most 4 counturs with a maximum of 9 vertices
    int contours[4][14] = {{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
        {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}};
    // loop over all edges
    int cnt = 0; // cnt is number of contours
    for (int e = 0; e < 12; e++) {
        if (segments[e][0] == -1)
            continue;

        contours[cnt][0] = e;
        int se = e;
        contours[cnt][1] = segments[e][0];
        int ce = segments[e][0];
        int pe = se;
        int cpos = 2;
        segments[e][0] = -1;
        while (true) {
            int s1 = segments[ce][0];
            int s2 = segments[ce][1];
            segments[ce][0] = -1;
            if (s1 == pe) {
                contours[cnt][cpos] = s2;
                pe = ce;
                ce = s2;
            } else if (s2 == pe) {
                contours[cnt][cpos] = s1;
                pe = ce;
                ce = s1;
            }
            cpos++;
            if (ce == se) {
                cnt++;
                break;
            }
        }

    }

    // compute intersection of opposite faces
    double ui[3][2];
    double vi[3][2];
    int    ic[3];
    for (int f = 0; f < 6; f += 2) {
        ic[f/2] = 0;
        double f1 = F[face_v[f][0]];
        double f2 = F[face_v[f][1]];
        double f3 = F[face_v[f][2]];
        double f4 = F[face_v[f][3]];
        double h1 = F[face_v[f+1][0]];
        double h2 = F[face_v[f+1][1]];
        double h3 = F[face_v[f+1][2]];
        double h4 = F[face_v[f+1][3]];
        const double a = (f1-f2)*(-h3+h4+h1-h2)-(h1-h2)*(-f3+f4+f1-f2);
        const double b = (i0-f1)*(-h3+h4+h1-h2)+(f1-f2)*(h3-h1)-(i0-h1)*(-f3+f4+f1-f2)-(h1-h2)*(f3-f1);
        const double c = (i0-f1)*(h3-h1)-(i0-h1)*(f3-f1);
        double d = b*b - 4*a*c;
        ui[f/2][0] = 0;
        ui[f/2][1] = 0;
        vi[f/2][0] = 0;
        vi[f/2][1] = 0;
        if (d > 0) {
            d = std::sqrt(d);
            double u1 = (-b-d) / (2*a);
            double u2 = (-b+d) / (2*a);
            if (u1 > u2) {
                double t = u1;
                u1 = u2;
                u2 = t;
            }
            double g1 = f1*(1-u1) + f2*u1;
            double g2 = f3*(1-u1) + f4*u1;
            double v1 = (i0 - g1)/(g2-g1);
            g1 = f1*(1-u2) + f2*u2;
            g2 = f3*(1-u2) + f4*u2;
            double v2 = (i0 - g1)/(g2-g1);
            if ((0 < u1 && u1 < 1) && (0 < v1 && v1 < 1)) {
                ui[f/2][ic[f/2]] = u1;
                vi[f/2][ic[f/2]] = v1;
                ic[f/2] += 1;
            }
            if ((0 < u2 && u2 < 1) && (0 < v2 && v2 < 1)) {
                ui[f/2][ic[f/2]] = u2;
                vi[f/2][ic[f/2]] = v2;
                ic[f/2] += 1;
            }
        }
    }


    // compute size of contours
    int csz[4] = {0,0,0,0};
    for (int ii = 0; ii < cnt; ii++) {
        int pos = 0;
        while (contours[ii][pos] != -1)
            pos++;
        csz[ii] = pos-1;
    }

    // triangulate contours
    // cnt is the number of contours, the last vertex is equal the first indicating the closed
    // contour, then there is a -1
    int tc = ic[0] + ic[1] + ic[2];

    // debugging
    auto c_index = [](Triangle& t) {
        if (t.v[0] == t.v[1])
            return false;
        if (t.v[0] == t.v[2])
            return false;
        if (t.v[1] == t.v[2])
            return false;
        return true;
    };

    // check if there is a tunnel or a contour with 12 vertices
    bool btunnel[] = {false,false,false,false};
    if (tc == 6) {
        m_ccases_tunnel += 1;
        // there is a tunnel, triangulate and mark used contours
        // if there is a tunnel, there are at most three contours
        // if there are only two contours, both build up the contour
        // if there are three contours, exclude the contour of length 3 which does not
        // belong to the tunnel
        btunnel[0] = true;
        btunnel[1] = true;
        if (cnt == 3) {
            btunnel[2] = true;
            // need coordinates
            double vc[12];
            vc[0]  = ecoord[0];
            vc[1]  = 1;
            vc[2]  = ecoord[2];
            vc[3]  = 0;
            vc[4]  = ecoord[4];
            vc[5]  = 1;
            vc[6]  = ecoord[6];
            vc[7]  = 0;
            vc[8]  = 0;
            vc[9]  = 1;
            vc[10] = 1;
            vc[11] = 0;
            for (int t = 0; t < cnt; t++) {
                if (csz[t] == 3) {
                    // check if countour does not belong to tunnel
                    // only 3, unroll loop
                    double umin = 2;
                    double umax = -2;
                    umin = (vc[contours[t][0]] < umin) ? vc[contours[t][0]] : umin;
                    umin = (vc[contours[t][1]] < umin) ? vc[contours[t][1]] : umin;
                    umin = (vc[contours[t][2]] < umin) ? vc[contours[t][2]] : umin;
                    umax = (vc[contours[t][0]] > umax) ? vc[contours[t][0]] : umax;
                    umax = (vc[contours[t][1]] > umax) ? vc[contours[t][1]] : umax;
                    umax = (vc[contours[t][2]] > umax) ? vc[contours[t][2]] : umax;
                    if (ui[0][0] > umax || ui[0][1] < umin) {
                        // outside the range, exclude
                        btunnel[t] = false;
                    }
                }
            }
        }
        // compute space hexagone
        // create 6 vertices
        double cvt[6][3];
        // face 1,2 and face 3,4 common coordinate is u
        double u = ui[0][0];
        double v = vi[0][0];
        double w = vi[1][0];
        uint p1 = 0;
        uint p2 = 0;
        uint p3 = 0;
        cvt[0][0] = u; cvt[0][1] = v; cvt[0][2] = w;
        // face 3,4 and face 5,6 common coord is w
        p3 = (std::fabs(w - vi[2][1]) < 0.000001) ? 1 : 0;

        // connect alternating in p1, p2 and p3
        // get new v coordinate from face 4,5
        v = ui[2][p3];
        cvt[1][0] = u; cvt[1][1] = v; cvt[1][2] = w;
        // get new u coordinate from face 0,1
        p1 = (p1+1)%2;
        u = ui[0][p1];
        cvt[2][0] = u; cvt[2][1] = v; cvt[2][2] = w; //cvt(3,:) = [u,v,w];
        // get new w coordinate from face 2,3
        p2 = (p2+1)%2;
        w = vi[1][p2];
        cvt[3][0] = u; cvt[3][1] = v; cvt[3][2] = w; //cvt(4,:) = [u,v,w];
        // get new v coordinate from face 4,5
        p3 = (p3+1)%2;
        v = ui[2][p3];
        cvt[4][0] = u; cvt[4][1] = v; cvt[4][2] = w; //cvt(5,:) = [u,v,w];
        // get nuew u coordinate from face 0,1
        p1 = (p1+1)%2;
        u = ui[0][p1];
        cvt[5][0] = u; cvt[5][1] = v; cvt[5][2] = w; //cvt(6,:) = [u,v,w];

        // compute triangulation
        // needs the space hexagon
        Point hex_p[6];
        for (int t = 0; t < 6; t++) {
            u = cvt[t][0]; v = cvt[t][1]; w = cvt[t][2];
            hex_p[t][0] = (1-w)*(1-v)*(1-u)*p[0][0] + (1-w)*(1-v)*u*p[1][0] + (1-w)*v*(1-u)*p[2][0] + (1-w)*v*u*p[3][0]
                        + w*(1-v)*(1-u)*p[4][0] + w*(1-v)*u*p[5][0] + w*v*(1-u)*p[6][0] + w*v*u*p[7][0];
            hex_p[t][1] = (1-w)*(1-v)*(1-u)*p[0][1] + (1-w)*(1-v)*u*p[1][1] + (1-w)*v*(1-u)*p[2][1] + (1-w)*v*u*p[3][1]
                        + w*(1-v)*(1-u)*p[4][1] + w*(1-v)*u*p[5][1] + w*v*(1-u)*p[6][1] + w*v*u*p[7][1];
            hex_p[t][2] = (1-w)*(1-v)*(1-u)*p[0][2] + (1-w)*(1-v)*u*p[1][2] + (1-w)*v*(1-u)*p[2][2] + (1-w)*v*u*p[3][2]
                        + w*(1-v)*(1-u)*p[4][2] + w*(1-v)*u*p[5][2] + w*v*(1-u)*p[6][2] + w*v*u*p[7][2];
        }
        // compute the three vertices of the triangulation and assign a global number
        Point  tv[3];
        Normal tn[3];
        int tg_edg[3];
        int tg_idx[3];
        for (int t = 0; t < 3; t++) {
            // compute vertices
            tv[t][0] = 0.5*(hex_p[2*t][0]+hex_p[2*t+1][0]);
            tv[t][1] = 0.5*(hex_p[2*t][1]+hex_p[2*t+1][1]);
            tv[t][2] = 0.5*(hex_p[2*t][2]+hex_p[2*t+1][2]);
            // compute normals
            u = (cvt[2*t][0]+cvt[2*t+1][0])/2; v = (cvt[2*t][1]+cvt[2*t+1][1])/2; w = (cvt[2*t][2]+cvt[2*t+1][2])/2;
            tn[t][0] = (1-w)*(1-v)*(1-u)*n[0][0] + (1-w)*(1-v)*u*n[1][0] + (1-w)*v*(1-u)*n[2][0] + (1-w)*v*u*n[3][0]
                     + w*(1-v)*(1-u)*n[4][0] + w*(1-v)*u*n[5][0] + w*v*(1-u)*n[6][0] + w*v*u*n[7][0];
            tn[t][1] = (1-w)*(1-v)*(1-u)*n[0][1] + (1-w)*(1-v)*u*n[1][1] + (1-w)*v*(1-u)*n[2][1] + (1-w)*v*u*n[3][1]
                     + w*(1-v)*(1-u)*n[4][1] + w*(1-v)*u*n[5][1] + w*v*(1-u)*n[6][1] + w*v*u*n[7][1];
            tn[t][2] = (1-w)*(1-v)*(1-u)*n[0][2] + (1-w)*(1-v)*u*n[1][2] + (1-w)*v*(1-u)*n[2][2] + (1-w)*v*u*n[3][2]
                     + w*(1-v)*(1-u)*n[4][2] + w*(1-v)*u*n[5][2] + w*v*(1-u)*n[6][2] + w*v*u*n[7][2];

            tg_edg[t] = int(m_cell_shift_factor*m_ugrid.global_index(i, j, k) + (3+t)); // this is the first interior vertex
            tg_idx[t] = (int)m_points.size();
            m_points.push_back(tv[t]);
            m_pnorms.push_back(tn[t]);
        }

        // triangulate:
        // 1. collect vertex pairs by using distance
        // 2. triangulate
        int dedgin[12]{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
        for (int t = 0; t < 4; t++) {
            if (!btunnel[t])
                continue;
            for (int r = 0; r < csz[t]; r++) {
                uint index;
                double dist = 4;
                for (int s = 0; s < 6; s++) {
                    double uval = ip[contours[t][r]][0] - hex_p[s][0];
                    double vval = ip[contours[t][r]][1] - hex_p[s][1];
                    double wval = ip[contours[t][r]][2] - hex_p[s][2];
                    double val = uval*uval + vval*vval + wval*wval;
                    if (dist > val) {
                        index = s;
                        dist = val;
                    }
                }
                dedgin[contours[t][r]] = index;
            }
        }

        // triangulate
        for (int t = 0; t < cnt; t++) {
            if (!btunnel[t])
                continue;
            for (int r = 0; r < csz[t]; r++) {
                // collect indices
                uint tid1 = contours[t][r];
                uint tid2 = contours[t][(r+1)%csz[t]];
                uint cid1 = dedgin[tid1];
                uint cid2 = dedgin[tid2];
                if (cid1 == cid2 || (cid1/2 == cid2/2)) {
                    Triangle tr;
                    tr.v[0] = vertices[tid1].g_idx;
                    tr.v[1] = vertices[tid2].g_idx;
                    tr.v[2] = tg_idx[cid1/2];
                    if (!c_index(tr))
                        std::cout << "wrong index\n";
                    m_triangles.push_back(tr);
                } else {
                    // measure distance and compute
                    // better triangulation: user shorter diag to divide
                    double uval = ip[tid1][0] - tv[cid2/2][0];
                    double vval = ip[tid1][1] - tv[cid2/2][1];
                    double wval = ip[tid1][2] - tv[cid2/2][2];
                    double val1 = uval*uval + vval*vval + wval*wval;
                    // second diag
                    uval = ip[tid2][0] - tv[cid1/2][0];
                    vval = ip[tid2][1] - tv[cid1/2][1];
                    wval = ip[tid2][2] - tv[cid1/2][2];
                    double val2 = uval*uval + vval*vval + wval*wval;
                    if (val1 < val2) {
                        Triangle tr;
                        tr.v[0] = vertices[tid1].g_idx;
                        tr.v[1] = tg_idx[cid2/2];
                        tr.v[2] = tg_idx[cid1/2];
                        if (!c_index(tr))
                            std::cout << "wrong index\n";
                        m_triangles.push_back(tr);
                        tr.v[0] = vertices[tid1].g_idx;
                        tr.v[1] = vertices[tid2].g_idx;
                        tr.v[2] = tg_idx[cid2/2];
                        if (!c_index(tr))
                            std::cout << "wrong index\n";
                        m_triangles.push_back(tr);
                    } else {
                        Triangle tr;
                        tr.v[0] = vertices[tid1].g_idx;
                        tr.v[1] = vertices[tid2].g_idx;
                        tr.v[2] = tg_idx[cid1/2];
                        m_triangles.push_back(tr);
                        tr.v[0] = vertices[tid2].g_idx;
                        tr.v[1] = tg_idx[cid2/2];
                        tr.v[2] = tg_idx[cid1/2];
                        m_triangles.push_back(tr);
                    }
                }
            }
        }
        // if there is a unique contour, then add triangle at the midpoint
        if (cnt == 1) {
            m_ccases_tunnel -= 1;
            m_ccases_12cont += 1;
            Triangle tr;
            tr.v[0] = tg_idx[0];
            tr.v[1] = tg_idx[1];
            tr.v[2] = tg_idx[2];
            if (!c_index(tr))
                std::cout << "wrong index\n";
            m_triangles.push_back(tr);
        }
    } else {
        // there is no tunnel
        if ((tc == 2 && ic[0] == 2) || (tc == 2 && ic[1] == 2)  ||  (tc == 2 && ic[2] == 2) || tc < 2) {
            for (int t = 0; t < cnt; t++) {
                // triangulation
                switch (csz[t]) {
                    case 3:
                    {
                        m_ccases_3 += 1;
                        Triangle tri;
                        tri.v[0] = vertices[contours[t][0]].g_idx; tri.v[1] = vertices[contours[t][1]].g_idx; tri.v[2] = vertices[contours[t][2]].g_idx;
                        m_triangles.push_back(tri);
                    }
                        break;
                    case 4:
                    {
                        m_ccases_4 += 1;
                        Triangle tri41;
                        Triangle tri42;
                        double d1 = distance(ip[contours[t][0]], ip[contours[t][2]]);
                        double d2 = distance(ip[contours[t][1]], ip[contours[t][3]]);
                        if (d1 < d2) {
                            tri41.v[0] = vertices[contours[t][0]].g_idx;
                            tri41.v[1] = vertices[contours[t][1]].g_idx;
                            tri41.v[2] = vertices[contours[t][2]].g_idx;
                            tri42.v[0] = vertices[contours[t][0]].g_idx;
                            tri42.v[1] = vertices[contours[t][2]].g_idx;
                            tri42.v[2] = vertices[contours[t][3]].g_idx;
                        } else {
                            tri41.v[0] = vertices[contours[t][0]].g_idx;
                            tri41.v[1] = vertices[contours[t][1]].g_idx;
                            tri41.v[2] = vertices[contours[t][3]].g_idx;
                            tri42.v[0] = vertices[contours[t][1]].g_idx;
                            tri42.v[1] = vertices[contours[t][2]].g_idx;
                            tri42.v[2] = vertices[contours[t][3]].g_idx;
                        }
                        m_triangles.push_back(tri41);
                        m_triangles.push_back(tri42);
                    }
                        break;
                    case 5:
                    {
                        m_ccases_5 += 1;
                        Triangle tri51;
                        Triangle tri52;
                        Triangle tri53;
                        tri51.v[0] = vertices[contours[t][0]].g_idx; tri51.v[1] = vertices[contours[t][1]].g_idx; tri51.v[2] = vertices[contours[t][2]].g_idx;
                        tri52.v[0] = vertices[contours[t][0]].g_idx; tri52.v[1] = vertices[contours[t][2]].g_idx; tri52.v[2] = vertices[contours[t][3]].g_idx;
                        tri53.v[0] = vertices[contours[t][0]].g_idx; tri53.v[1] = vertices[contours[t][3]].g_idx; tri53.v[2] = vertices[contours[t][4]].g_idx;
                        m_triangles.push_back(tri51);
                        m_triangles.push_back(tri52);
                        m_triangles.push_back(tri53);
                    }
                        break;
                    case 6:
                        m_ccases_6 += 1;
                        Triangle tri61;
                        Triangle tri62;
                        Triangle tri63;
                        Triangle tri64;
                        tri61.v[0] = vertices[contours[t][0]].g_idx; tri61.v[1] = vertices[contours[t][4]].g_idx; tri61.v[2] = vertices[contours[t][5]].g_idx;
                        tri62.v[0] = vertices[contours[t][0]].g_idx; tri62.v[1] = vertices[contours[t][2]].g_idx; tri62.v[2] = vertices[contours[t][4]].g_idx;
                        tri63.v[0] = vertices[contours[t][0]].g_idx; tri63.v[1] = vertices[contours[t][1]].g_idx; tri63.v[2] = vertices[contours[t][2]].g_idx;
                        tri64.v[0] = vertices[contours[t][2]].g_idx; tri64.v[1] = vertices[contours[t][3]].g_idx; tri64.v[2] = vertices[contours[t][4]].g_idx;
                        m_triangles.push_back(tri61);
                        m_triangles.push_back(tri62);
                        m_triangles.push_back(tri63);
                        m_triangles.push_back(tri64);

                } // switch (csz(t))
            }
        } else {
            // ambiguous face
            // compute inner vertex
            double u{0},v{0},w{0};
            if (ic[0] == 2) {
                // face 0,1 has not asymptotes, common coordinate is w
                u = ui[1][0];
                v = ui[2][0];
                w = vi[1][0];
            } else if (ic[1] == 2) {
                // face 2,3 have no asymptotes, common coordinate is v
                u = ui[0][0];
                v = vi[0][0];
                w = vi[2][0];
            } else if (ic[2] == 2){
                // face 4,5 have no asymptotes, common coordinate is u
                u = ui[0][0];
                v = vi[0][0];
                w = vi[1][0];
            } else {
                u = ui[0][0] + ui[0][1] + ui[1][0] + ui[1][1];
                v = vi[0][0] + vi[0][1] + ui[2][0] + ui[2][1];
                w = vi[1][0] + vi[1][1] + vi[2][0] + vi[2][1];
                u = u/(ic[0] + ic[1]);
                v = v/(ic[0] + ic[2]);
                w = w/(ic[1] + ic[2]);

            }
            // compute vertex: trilinear interpolate
			Point cv;
			Normal cn;
            cv[0] = (1-w)*((1-v)*((1-u)*p[0][0]+u*p[1][0]) + v*((1-u)*p[2][0]+u*p[3][0]))
                  + w*((1-v)*((1-u)*p[4][0]+u*p[5][0]) + v*((1-u)*p[6][0]+u*p[7][0]));
            cv[1] = (1-w)*((1-v)*((1-u)*p[0][1]+u*p[1][1]) + v*((1-u)*p[2][1]+u*p[3][1]))
                  + w*((1-v)*((1-u)*p[4][1]+u*p[5][1]) + v*((1-u)*p[6][1]+u*p[7][1]));
            cv[2] = (1-w)*((1-v)*((1-u)*p[0][2]+u*p[1][2]) + v*((1-u)*p[2][2]+u*p[3][2]))
                  + w*((1-v)*((1-u)*p[4][2]+u*p[5][2]) + v*((1-u)*p[6][2]+u*p[7][2]));

            cn[0] = (1-w)*((1-v)*((1-u)*n[0][0]+u*n[1][0]) + v*((1-u)*n[2][0]+u*n[3][0]))
                  + w*((1-v)*((1-u)*n[4][0]+u*n[5][0]) + v*((1-u)*n[6][0]+u*n[7][0]));
            cn[1] = (1-w)*((1-v)*((1-u)*n[0][1]+u*n[1][1]) + v*((1-u)*n[2][1]+u*n[3][1]))
                  + w*((1-v)*((1-u)*n[4][1]+u*n[5][1]) + v*((1-u)*n[6][1]+u*n[7][1]));
            cn[2] = (1-w)*((1-v)*((1-u)*n[0][2]+u*n[1][2]) + v*((1-u)*n[2][2]+u*n[3][2]))
                  + w*((1-v)*((1-u)*n[4][2]+u*n[5][2]) + v*((1-u)*n[6][2]+u*n[7][2]));

            bool flag = false;
            for (int s = 0; s < cnt; s++) { // for each contour
                if (csz[s] == 3) {
                    m_ccases_3 += 1;
                    Triangle tri;
                    tri.v[0] = vertices[contours[s][0]].g_idx;
                    tri.v[1] = vertices[contours[s][1]].g_idx;
                    tri.v[2] = vertices[contours[s][2]].g_idx;
                    m_triangles.push_back(tri);
                } else if (csz[s] == 5) {
                    m_ccases_5 += 1;
                    Triangle tri51;
                    Triangle tri52;
                    Triangle tri53;
                    tri51.v[0] = vertices[contours[s][0]].g_idx; tri51.v[1] = vertices[contours[s][1]].g_idx; tri51.v[2] = vertices[contours[s][2]].g_idx;
                    tri52.v[0] = vertices[contours[s][0]].g_idx; tri52.v[1] = vertices[contours[s][2]].g_idx; tri52.v[2] = vertices[contours[s][3]].g_idx;
                    tri53.v[0] = vertices[contours[s][0]].g_idx; tri53.v[1] = vertices[contours[s][3]].g_idx; tri53.v[2] = vertices[contours[s][4]].g_idx;
                    m_triangles.push_back(tri51);
                    m_triangles.push_back(tri52);
                    m_triangles.push_back(tri53);
                }
                else {
                    const int c_gidx = m_points.size();
                    m_points.push_back(cv);
                    m_pnorms.push_back(cn);
                    for (int l = 0; l < csz[s]; l++) {
                        int v0 = l;
                        int v1 = (l+1)%csz[s];
                        Triangle tri;
                        tri.v[0] = vertices[contours[s][v0]].g_idx;
                        tri.v[1] = vertices[contours[s][v1]].g_idx;
                        tri.v[2] = c_gidx;
                        m_triangles.push_back(tri);
                    }
                    // statistics
                    switch (csz[s]) {
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
                            std::cout << "ERROR: wrong contour length: " << csz[s] << std::endl;
                            break;
                    }

                }

            }
        }
    }
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
			if (val >= i0) {
				set_segm(e3, 0, e0, segm_);
				set_segm(e0, 1, e3, segm_);
				set_segm(e1, 0, e2, segm_);
				set_segm(e2, 1, e1, segm_);
			}
			else {
				set_segm(e1, 0, e0, segm_);
				set_segm(e0, 1, e1, segm_);
				set_segm(e3, 0, e2, segm_);
				set_segm(e2, 1, e3, segm_);
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
			if (val >= i0){
				set_segm(e0, 0, e1, segm_);
				set_segm(e1, 1, e0, segm_);
				set_segm(e2, 0, e3, segm_);
				set_segm(e3, 1, e2, segm_);
			}
			else {
				set_segm(e0, 0, e3, segm_);
				set_segm(e3, 1, e0, segm_);
				set_segm(e2, 0, e1, segm_);
				set_segm(e1, 1, e2, segm_);
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
		g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
		g2 = F[2] * (1 - ui[1]) + F[3] * ui[1];
		vi[1] = (i0 - g1) / (g2 - g1);
		// compute w-coordinates of solutions
		g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
		g2 = F[4] * (1 - ui[0]) + F[5] * ui[0];
		wi[0] = (i0 - g1) / (g2 - g1);
		g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
		g2 = F[4] * (1 - ui[1]) + F[5] * ui[1];
		wi[1] = (i0 - g1) / (g2 - g1);

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
		}


		// triangulate contours with inner hexagon
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
			Triangle t1;
			Triangle t2;
			Triangle t3;
			Triangle t4;
			t1.v[0] = tg_idx[0]; t1.v[1] = tg_idx[2]; t1.v[2] = tg_idx[1];
			t2.v[0] = tg_idx[2]; t2.v[1] = tg_idx[4]; t2.v[2] = tg_idx[3];
			t3.v[0] = tg_idx[0]; t3.v[1] = tg_idx[5]; t3.v[2] = tg_idx[4];
			t4.v[0] = tg_idx[0]; t4.v[1] = tg_idx[4]; t4.v[2] = tg_idx[2];
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
			m_points.push_back(ip);
			m_pnorms.push_back(in);
			// loop over the contorus
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



void tmc::MarchingCubes::connectivity( )
{
    const int nr_v = (int)m_points.size();
    // debugging
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

    // reorient triangles
    std::vector<std::vector<int>> v_tris(nr_v);
    for (auto t : m_triangles) {
        v_tris[t.v[0]].push_back(t.id);
        v_tris[t.v[1]].push_back(t.id);
        v_tris[t.v[2]].push_back(t.id);
    }

    // check vertex fan
    for (int v = 0; v < nr_v; v++) {
        std::vector<int> v_fan;
        for (auto t : v_tris[v]) {
            int t0 = m_triangles[t].v[0]; //m_ttriangles[t].v[0];
            int t1 = m_triangles[t].v[1];
            int t2 = m_triangles[t].v[2];
            if (t0 != v) v_fan.push_back(t0);
            if (t1 != v) v_fan.push_back(t1);
            if (t2 != v) v_fan.push_back(t2);
        }
        // check triangle fan
        if (v_fan.size()%2 != 0) {
            std::cout << "ERROR: wrong number of vertices in triangle fan\n";
        }
        std::sort (v_fan.begin(), v_fan.end());
        std::vector<int> v_neigh(v_fan.size());
        if (v_fan[0] == v_fan[1])
            v_neigh[0] = 2;
        else
            v_neigh[0] = 1;
        for (int n = 1; n < (v_fan.size()-1); n++) {
            if (v_fan[n] == v_fan[n-1] || v_fan[n] == v_fan[n+1])
                v_neigh[n] = 2;
            else
                v_neigh[n] = 1;
        }
        if (v_fan[v_fan.size()-1] == v_fan[v_fan.size()-2])
            v_neigh[v_fan.size()-1] = 2;
        else
            v_neigh[v_fan.size()-1] = 1;

        // check nr. of vertices
        int v_one{0};
        for (auto n : v_neigh) {
            if (n == 0) {
                std::cout << "ERROR: vertex v0: " << v << ", no neighbors ZERO\n";
            } else if (n == 1) {
                m_vboundary[v] = true;
                v_one = v_one + 1;
            }
        }
        if (v_one == 1 || v_one > 2) {
            std::cout << "ERROR: wrong number of neighbors\n";
        }
    }


    auto v_compare = [] (const int v1,const int v2, const int t0,const int t1,const int t2) {
        bool face0 = (v1 == t1 && v2 == t2) || (v1 == t2 && v2 == t1);
        bool face1 = (v1 == t0 && v2 == t2) || (v1 == t2 && v2 == t0);
        bool face2 = (v1 == t1 && v2 == t0) || (v1 == t0 && v2 == t1);

        if (face0 || face1 || face2) {
            return true;
        } else {
            return false;
        }
    };

    // create face neighbors
    // loop over all triangle and set neighors
    std::vector<std::array<int,3>> t_neigh(m_triangles.size());
    for (auto t : m_triangles) {
        int v0 = t.v[0];
        int v1 = t.v[1];
        int v2 = t.v[2];
        t_neigh[t.id][0] = -1;
        t_neigh[t.id][1] = -1;
        t_neigh[t.id][2] = -1;


        // face 0
        bool t_topol{false};
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
                std::cout << "ERROR: triangle " << t.id << ", face 0: vertex v1 is not boundary!\n";
            }
            if (!m_vboundary[v2]) {
                std::cout << "ERROR: triangle " << t.id << ", face 0: vertex v2 is not boundary!\n";
            }
        }
        t_topol = f_neigh;

        // face 1
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
                std::cout << "ERROR: triangle " << t.id << ", face 1: vertex v0 is not boundary!\n";
            }
            if (!m_vboundary[v2]) {
                std::cout << "ERROR: triangle " << t.id << ", face 1: vertex v2 is not boundary!\n";
            }
        }
        t_topol = f_neigh;

        // face 2
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
                std::cout << "ERROR: triangle " << t.id << ", face 2: vertex v0 is not boundary!\n";
            }
            if (!m_vboundary[v1]) {
                std::cout << "ERROR: triangle " << t.id << ", face 2: vertex v1 is not boundary!\n";
            }
        }
        t_topol = f_neigh;

        if (!t_topol) {
            std::cout << "ERROR: triangle " << t.id << " has no neighbor at face\n";
        }
    }

    // check neighborhood
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

    int t_nrcomps{0};
    for (auto t : m_triangles) {
        if (t_oriented[t.id] == true) {
            continue;
        }
        if (t.id == 20) {
            std::cout << "stop here\n";
        }
        t_nrcomps = t_nrcomps + 1;
        // triangle still not processed
        t_oriented[t.id] = true;
        int v0 = t.v[0];
        int v1 = t.v[1];
        int v2 = t.v[2];
        int n0 = t_neigh[t.id][0];
        int n1 = t_neigh[t.id][1];
        int n2 = t_neigh[t.id][2];
        // face 0
        if (n0 != -1) {
            if (t_oriented[n0])
                std::cout << "for triangle " << t.id << ", neighbor " << n0 << " is oriented" << std::endl;
            t_compare(v1,v2,m_triangles[n0]);
            t_queue.push(n0);
        }
        // face 1
        if (n1 != -1) {
            if (t_oriented[n1])
                std::cout << "for triangle " << t.id << ", neighbor " << n1 << " is oriented" << std::endl;
            t_compare(v2,v0,m_triangles[n1]);
            t_queue.push(n1);
        }

        // face 2
        if (n2 != -1) {
            if (t_oriented[n2])
                std::cout << "for triangle " << t.id << ", neighbor " << n2 << " is oriented" << std::endl;
            t_compare(v0,v1,m_triangles[n2]);
            t_queue.push(n2);
        }
        // set stack
        while (!t_queue.empty()) {
            int c_tri = t_queue.top();
            t_queue.pop();
            if (!t_oriented[c_tri]) {
                // triangle still not processed
                t_oriented[c_tri] = true;
                Triangle t = m_triangles[c_tri];
                int t0 = t.v[0];
                int t1 = t.v[1];
                int t2 = t.v[2];
                int tn0 = t_neigh[t.id][0];
                int tn1 = t_neigh[t.id][1];
                int tn2 = t_neigh[t.id][2];
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

