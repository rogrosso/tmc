//
//  Created by Roberto Grosso on 05.06.16.
//  Copyright © 2016 Roberto Grosso. All rights reserved.
//

#include "MarchingCubes.h"

//
//  MarchingCubes.cpp
//
//  Created by Roberto Grosso on 07/03/16.
//  Copyright © 2016 Roberto Grosso. All rights reserved.
//

#include "MarchingCubes.h"


void
tmc::MarchingCubes::operator() (const std::string& i_file,const std::string& o_objF,const std::string& o_offF)
{
    std::cout << " ... reading data \n";
    
    std::FILE* f = std::fopen(i_file.c_str(), "r");
    short x_size;
    short y_size;
    short z_size;
    std::fread(&x_size, sizeof(short), 1, f);
    std::fread(&y_size, sizeof(short), 1, f);
    std::fread(&z_size, sizeof(short), 1, f);
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
    bb[0] = Point{0,0,0};
    bb[1] = Point{xmax,0,0};
    bb[2] = Point{0,ymax,0};
    bb[3] = Point{xmax,ymax,0};
    bb[4] = Point{0,0,zmax};
    bb[5] = Point{xmax,0,zmax};
    bb[6] = Point{0,ymax,zmax};
    bb[7] = Point{xmax,ymax,zmax};
    
    m_ugrid.init(m_nx, m_ny, m_nz,bb);
    int m_size = m_nx * m_ny * m_nz;
    unsigned short* t_buff = new unsigned short[x_size*y_size*z_size];
    std::fread(&t_buff[0],sizeof(unsigned short),x_size*y_size*z_size,f);
    for (int i = 0; i < m_size; i++) {
        m_ugrid.scalar(i,double(t_buff[i]));
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
    t_mc(1000);

    
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
//  IMPLEMENTATION
//*******************************************************************************************************************************************
//void
//tmc::MarchingCubes::create_test()
//{
//    const int nx = m_ugrid.x_size();
//    const int ny = m_ugrid.y_size();
//    const int nz = m_ugrid.z_size();
//    Point center{0.5,0.5,0.5};
//    for (int k = 0; k < nz; k++) {
//        for (int j = 0; j < ny; j++) {
//            for (int i = 0; i < nx; i++) {
//                Point p = m_ugrid.point(i,j,k);
//                double r2 = std::sqrt(dot(p-center,p-center));
//                m_ugrid.scalar(i,j,k,r2);
//            }
//        }
//    }
//}


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
                    t_slice(i,j,k,i0,u,p,n,i_case);
                    //p_slice(i,j,k,i0,u,p,n,i_case);
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
    unsigned short face_e_[6] = {12816, 30292, 33936, 46754, 34739, 38305};
    const int face_v[6][4] = {{0,1,2,3},{4,5,6,7},{0,1,4,5},{2,3,6,7},{0,2,4,6},{1,3,5,7}};
    unsigned short face_v_[6] = {12816, 30292, 21520, 30258, 25632, 30001};
    
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
    for (int i = 0; i < cnt; i++) {
        int pos = 0;
        while (contours[i][pos] != -1)
            pos++;
        csz[i] = pos-1;
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
            Point cv{0,0,0};
            Normal cn{0,0,0};
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
//void
//tmc::MarchingCubes::p_slice(const int i,const int j,const int k,const double i0,double* F,Point* p,Normal* n,const int i_case)
//{
//    unsigned short face_e_[6] = {12816, 30292, 33936, 46754, 34739, 38305};
//    unsigned short face_v_[6] = {12816, 30292, 21520, 30258, 25632, 30001};
//    
//    // reading edge from face
//    auto get_face_e = [face_e_](const int f,const int e) { return ((face_e_[f]>>(4*e))&0xF); };
//    auto get_face_v = [face_v_](const int f,const int e) { return ((face_v_[f]>>(4*e))&0xF); };
//    
//    // there are 12 edges, assign to each vertex three edges, the global edge numbering
//    // consist of 3*global_vertex_id + edge_offset.
//    const unsigned long long gei_pattern_ = 670526590282893600ull;
//    const uint axis_mask = 1;
//    const uint offs_mask = 3;
//    //    const int global_edge_id[][4] = {{0,0,0,0},{1,0,0,1},{0,1,0,0},{0,0,0,1},
//    //        {0,0,1,0},{1,0,1,1},{0,1,1,0},{0,0,1,1},
//    //        {0,0,0,2},{1,0,0,2},{1,1,0,2},{0,1,0,2}};
//    // the end vertices of an edge
//    //    int l_edges[12][2] = {{0,1}, {1,3}, {2,3}, {0,2},
//    //        {4,5}, {5,7}, {6,7}, {4,6},
//    //        {0,4}, {1,5}, {3,7}, {2,6}};
//    const unsigned char l_edges_[12] = {16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98};
//    // the two faces sharing an edge
//    struct vertex {
//        uint g_edg; // global edge id
//        uint g_idx; // final index in vertex list
//    };
//    // there might be up to 15 vertices and normals including three vertices or normals at the interior
//    std::vector<vertex> vertices(12);
//    std::vector<Point>  ip(12);
//    std::vector<Normal> in(12);
//    std::vector<double> ecoord{0,0,0,0,0,0,0,0,0,0,0,0}; // there are 12 coordinates along the edges, these are the intersections
//    
//    
//    // edge intersection pattern
//    const ushort e_pattern_ = edge_pattern[i_case];
//    
//    
//    // collect vertices
//    ushort   flag{1};
//    for (int eg = 0; eg < 12; eg++) {
//        vertices[eg].g_edg = -1;
//        vertices[eg].g_idx = -1;
//        if (flag & e_pattern[i_case]) {
//            // the edge global index is given by the vertex global index + the edge offset
//            uint shift = 5*eg;
//            const int ix = i + (int)((gei_pattern_>>shift)&1); // global_edge_id[eg][0];
//            const int iy = j + (int)((gei_pattern_>>(shift+1))&1); // global_edge_id[eg][1];
//            const int iz = k + (int)((gei_pattern_>>(shift+2))&1); // global_edge_id[eg][2];
//            const int off_val = (int)((gei_pattern_>>(shift+3))&3);
//            
//            vertices[eg].g_edg = int(m_cell_shift_factor*m_ugrid.global_index(ix, iy, iz) + off_val);
//            
//            // generate vertex here, do not care at this point if vertex already exist
//            //int* vert = l_edges[eg];
//            // interpolation weight
//            uint v0 = (l_edges_[eg]&0xF);
//            uint v1 = (l_edges_[eg]>>4)&0xF;
//            //            int v0 = vert[0];
//            //            int v1 = vert[1];
//            double l = (i0 - F[v0]) / (F[v1] - F[v0]);
//            ecoord[eg] = l;
//            // interpolate vertex
//            ip[eg] = (1 - l)*p[v0] + l*p[v1];
//            // interpolate normal
//            in[eg] = (1 - l)*n[v0] + l*n[v1];
//            in[eg] /= std::sqrt(dot(in[eg], in[eg]));
//            // set vertex index
//            int g_idx = set_tvertex(vertices[eg].g_edg);
//            if (g_idx == -1) {
//                TVertex tvertex;
//                tvertex.p = ip[eg];
//                tvertex.n = in[eg];
//                g_idx = (int)m_tvertices.size();
//                tvertex.g_idx = g_idx;
//                tvertex.g_edg = vertices[eg].g_edg;
//                m_tvertices.push_back(tvertex);
//                EdgeVertex ev{tvertex.g_edg,tvertex.g_idx};
//                m_ehash[hashbucket_index(m_nbuckets, tvertex.g_edg)].push_back(ev);
//            }
//            vertices[eg].g_idx = g_idx;
//        }
//        flag <<=1;
//    }
//    
//    // build up segments
//    unsigned char segm_[12]={0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF};
//    auto set_segm = [] (const int& e, const int& val,unsigned char segm_[12]) {
//        if ((segm_[e]&0xF) == 0xF)
//            segm_[e] = (unsigned char)val&0xF;
//        else
//            segm_[e] |= val<<4;
//    };
//    auto get_segm = [] (const int& e, const int pos,unsigned char segm_[12]) {
//        if (pos == 1)
//            return (int)((segm_[e]>>4)&0xF);
//        else
//            return (int)(segm_[e]&0xF);
//    };
//    
//    for (int f = 0; f < 6; f++) {
//        int pos = 0;
//        int p_ = 0;
//        int s[4] = {0,0,0,0};
//        unsigned short s_ = 0xFFFF;
//        auto set_e = [] (const int e, unsigned short& s_) {
//            if ((s_&0xF) == 0xF)
//                s_ &= 0xFFF0^(e);
//            else if (((s_>>4)&0xF) == 0xF)
//                s_ &= 0xFF0F^(e<<4);
//            else if (((s_>>8)&0xF) == 0xF)
//                s_ &= 0xF0FF^(e<<8);
//            else
//                s_ &= 0x0FFF^(e<<12);
//        };
//        auto get_e = [] (const int p,unsigned short s_) { return (int)((s_>>(4*p))&0xF); };
//        
//        for (int eg = 0; eg < 4; eg++) {
//            ushort eid_ = get_face_e(f,eg);
//            if ((e_pattern_>>eid_)&1) {
//                p_ += 1;
//                set_e(eid_,s_);
//            }
//        }
//        
//        // test new revolutionary method
//        if (p_ == 2) {
//            set_segm(get_e(0,s_),get_e(1,s_),segm_);
//            set_segm(get_e(1,s_),get_e(0,s_),segm_);
//        } else if (p_ == 4){
//            const double f0 = F[get_face_v(f,0)];
//            const double f1 = F[get_face_v(f,1)];
//            const double f2 = F[get_face_v(f,2)];
//            const double f3 = F[get_face_v(f,3)];
//            const double eta = f0+f3-f1-f2;
//            const double u0 = (f0-f2)/eta;
//            if (ecoord[get_face_e(f,0)] < u0) {
//                set_segm(get_face_e(f,0),get_face_e(f,3),segm_);
//                set_segm(get_face_e(f,3),get_face_e(f,0),segm_);
//                set_segm(get_face_e(f,1),get_face_e(f,2),segm_);
//                set_segm(get_face_e(f,2),get_face_e(f,1),segm_);
//            } else {
//                set_segm(get_face_e(f,0),get_face_e(f,1),segm_);
//                set_segm(get_face_e(f,1),get_face_e(f,0),segm_);
//                set_segm(get_face_e(f,2),get_face_e(f,3),segm_);
//                set_segm(get_face_e(f,3),get_face_e(f,2),segm_);
//            }
//        }
//    } //
//    
//    unsigned long long c_ = 0xFFFFFFFFFFFFFFFF;
//    auto set_c = [] (const int e,const int val,unsigned long long &c_) {
//        c_ &= ~(0xFull<<(e*4));
//        c_ |= (((unsigned long long)val)<<(e*4));
//    };
//    auto get_c = [] (const int e, unsigned long long c_) {
//        return (int)((c_>>(4*e))&0xF);
//    };
//    unsigned char c_sz[4]={0,0,0,0};
//    int cnt_ = 0; // cnt is number of contours
//    int pos_ = 0;
//    for (int e = 0; e < 12; e++) {
//        if (segm_[e] == 0xFF) {
//            continue;
//        }
//        
//        // contour found
//        unsigned char le = get_segm(e,0,segm_);
//        unsigned char ce = e;
//        unsigned char ne = get_segm(e,1,segm_);
//        
//        set_c(pos_+c_sz[cnt_],le,c_);
//        set_c(pos_+c_sz[cnt_]+1,ce,c_);
//        set_c(pos_+c_sz[cnt_]+2,ne,c_);
//        c_sz[cnt_] += 3;
//        segm_[e] = 0xFF;
//        while (true) {
//            int s1_ = get_segm(ne,0,segm_);
//            int s2_ = get_segm(ne,1,segm_);
//            segm_[ne] = 0xFF;
//            if (s1_ == ce) {
//                if (s2_ == le) {
//                    segm_[le] = 0xFF;
//                    pos_ += c_sz[cnt_];
//                    cnt_ += 1;
//                    break;
//                }
//                set_c(pos_+c_sz[cnt_],s2_,c_);
//                ce = ne;
//                ne = s2_;
//                c_sz[cnt_] += 1;
//            } else if (s2_ == ce) {
//                if (s1_ == le) {
//                    segm_[le] = 0xFF;
//                    pos_ += c_sz[cnt_];
//                    cnt_ += 1;
//                    break;
//                }
//                set_c(pos_+c_sz[cnt_],s1_,c_);
//                ce = ne;
//                ne = s1_;
//                c_sz[cnt_] += 1;
//            }
//        }
//        
//    }
//    
//    // compute intersection of opposite faces
//    double ui[3][2];
//    double vi[3][2];
//    int    ic[3];
//    for (int f = 0; f < 6; f += 2) {
//        ic[f/2] = 0;
//        double f1 = F[get_face_v(f,0)];
//        double f2 = F[get_face_v(f,1)];
//        double f3 = F[get_face_v(f,2)];
//        double f4 = F[get_face_v(f,3)];
//        double h1 = F[get_face_v(f+1,0)];
//        double h2 = F[get_face_v(f+1,1)];
//        double h3 = F[get_face_v(f+1,2)];
//        double h4 = F[get_face_v(f+1,3)];
//        const double a = (f1-f2)*(-h3+h4+h1-h2)-(h1-h2)*(-f3+f4+f1-f2);
//        const double b = (i0-f1)*(-h3+h4+h1-h2)+(f1-f2)*(h3-h1)-(i0-h1)*(-f3+f4+f1-f2)-(h1-h2)*(f3-f1);
//        const double c = (i0-f1)*(h3-h1)-(i0-h1)*(f3-f1);
//        double d = b*b - 4*a*c;
//        ui[f/2][0] = 0;
//        ui[f/2][1] = 0;
//        vi[f/2][0] = 0;
//        vi[f/2][1] = 0;
//        if (d > 0) {
//            d = std::sqrt(d);
//            double u1 = (-b-d) / (2*a);
//            double u2 = (-b+d) / (2*a);
//            if (u1 > u2) {
//                double t = u1;
//                u1 = u2;
//                u2 = t;
//            }
//            double g1 = f1*(1-u1) + f2*u1;
//            double g2 = f3*(1-u1) + f4*u1;
//            double v1 = (i0 - g1)/(g2-g1);
//            g1 = f1*(1-u2) + f2*u2;
//            g2 = f3*(1-u2) + f4*u2;
//            double v2 = (i0 - g1)/(g2-g1);
//            if ((0 < u1 && u1 < 1) && (0 < v1 && v1 < 1)) {
//                ui[f/2][ic[f/2]] = u1;
//                vi[f/2][ic[f/2]] = v1;
//                ic[f/2] += 1;
//            }
//            if ((0 < u2 && u2 < 1) && (0 < v2 && v2 < 1)) {
//                ui[f/2][ic[f/2]] = u2;
//                vi[f/2][ic[f/2]] = v2;
//                ic[f/2] += 1;
//            }
//        }
//    }
//    
//    
//    // triangulate contours
//    // cnt is the number of contours, the last vertex is equal the first indicating the closed
//    // contour, then there is a -1
//    int tc = ic[0] + ic[1] + ic[2];
//    
//    // check if there is a tunnel or a contour with 12 vertices
//    bool tun_[4] = {false,false,false,false};
//    if (tc == 6) {
//        m_ccases_tunnel += 1;
//        // there is a tunnel, triangulate and mark used contours
//        // if there is a tunnel, there are at most three contours
//        // if there are only two contours, both build up the contour
//        // if there are three contours, exclude the contour of length 3 which does not
//        // belong to the tunnel
//        // check which contour does not belong to the tunnel
//        tun_[0] = true;
//        if (cnt_ == 2) {
//            tun_[1] = true;
//        } else if (cnt_ == 3) {
//            tun_[2] = true;
//            double vc[12];
//            vc[0]  = ecoord[0];
//            vc[1]  = 1;
//            vc[2]  = ecoord[2];
//            vc[3]  = 0;
//            vc[4]  = ecoord[4];
//            vc[5]  = 1;
//            vc[6]  = ecoord[6];
//            vc[7]  = 0;
//            vc[8]  = 0;
//            vc[9]  = 1;
//            vc[10] = 1;
//            vc[11] = 0;
//            pos_ = 0;
//            for (int t = 0; t < cnt_; t++) {
//                if (c_sz[t] == 3) {
//                    // check if countour does not belong to tunnel
//                    // only 3, unroll loop
//                    double umin = 2;
//                    double umax = -2;
//                    uint e0 = get_c(pos_,c_);
//                    uint e1 = get_c(pos_+1,c_);
//                    uint e2 = get_c(pos_+2,c_);
//                    umin = (vc[e0] < umin) ? vc[e0] : umin;
//                    umin = (vc[e1] < umin) ? vc[e1] : umin;
//                    umin = (vc[e2] < umin) ? vc[e2] : umin;
//                    umax = (vc[e0] > umax) ? vc[e0] : umax;
//                    umax = (vc[e1] > umax) ? vc[e1] : umax;
//                    umax = (vc[e2] > umax) ? vc[e2] : umax;
//                    if (ui[0][0] > umax || ui[0][1] < umin) {
//                        // outside the range, exclude
//                        tun_[t] = false;
//                        break;
//                    }
//                }
//                pos_ += c_sz[t];
//            }
//        }
//        // compute space hexagone
//        // create 6 vertices
//        double cvt[6][3];
//        // face 1,2 and face 3,4 common coordinate is u
//        double u = ui[0][0];
//        double v = vi[0][0];
//        double w = vi[1][0];
//        uint p1 = 0;
//        uint p2 = 0;
//        uint p3 = 0;
//        cvt[0][0] = u; cvt[0][1] = v; cvt[0][2] = w;
//        // face 3,4 and face 5,6 common coord is w
//        p3 = (std::fabs(w - vi[2][1]) < 0.00005) ? 1 : 0;
//        
//        // connect alternating in p1, p2 and p3
//        // get new v coordinate from face 4,5
//        v = ui[2][p3];
//        cvt[1][0] = u; cvt[1][1] = v; cvt[1][2] = w;
//        // get new u coordinate from face 0,1
//        p1 = (p1+1)%2;
//        u = ui[0][p1];
//        cvt[2][0] = u; cvt[2][1] = v; cvt[2][2] = w; //cvt(3,:) = [u,v,w];
//        // get new w coordinate from face 2,3
//        p2 = (p2+1)%2;
//        w = vi[1][p2];
//        cvt[3][0] = u; cvt[3][1] = v; cvt[3][2] = w; //cvt(4,:) = [u,v,w];
//        // get new v coordinate from face 4,5
//        p3 = (p3+1)%2;
//        v = ui[2][p3];
//        cvt[4][0] = u; cvt[4][1] = v; cvt[4][2] = w; //cvt(5,:) = [u,v,w];
//        // get nuew u coordinate from face 0,1
//        p1 = (p1+1)%2;
//        u = ui[0][p1];
//        cvt[5][0] = u; cvt[5][1] = v; cvt[5][2] = w; //cvt(6,:) = [u,v,w];
//        
//        // compute triangulation
//        // needs the space hexagon
//        Point hex_p[6];
//        for (int t = 0; t < 6; t++) {
//            u = cvt[t][0]; v = cvt[t][1]; w = cvt[t][2];
//            hex_p[t] = (1-w)*(1-v)*(1-u)*p[0] + (1-w)*(1-v)*u*p[1] + (1-w)*v*(1-u)*p[2] + (1-w)*v*u*p[3] + w*(1-v)*(1-u)*p[4] + w*(1-v)*u*p[5] + w*v*(1-u)*p[6] + w*v*u*p[7];
//        }
//        // compute the three vertices of the triangulation and assign a global number
//        Point  tv[3];
//        Normal tn[3];
//        int tg_edg[3];
//        int tg_idx[3];
//        for (int t = 0; t < 3; t++) {
//            // compute vertices
//            tv[t] = 0.5*(hex_p[2*t]+hex_p[2*t+1]);
//            // compute normals
//            u = (cvt[2*t][0]+cvt[2*t+1][0])/2; v = (cvt[2*t][1]+cvt[2*t+1][1])/2; w = (cvt[2*t][2]+cvt[2*t+1][2])/2;
//            tn[t] = (1-w)*(1-v)*(1-u)*n[0] + (1-w)*(1-v)*u*n[1] + (1-w)*v*(1-u)*n[2] + (1-w)*v*u*n[3] + w*(1-v)*(1-u)*n[4] + w*(1-v)*u*n[5] + w*v*(1-u)*n[6] + w*v*u*n[7];
//            
//            tg_edg[t] = int(m_cell_shift_factor*m_ugrid.global_index(i, j, k) + (3+t)); // this is the first interior vertex
//            tg_idx[t] = set_tvertex(tg_edg[t]);
//            if (tg_idx[t] == -1) {
//                TVertex tvertex;
//                tvertex.p = tv[t];
//                tvertex.n = tn[t];
//                tg_idx[t] = (int)m_tvertices.size();
//                tvertex.g_idx = tg_idx[t];
//                tvertex.g_edg = tg_edg[t];
//                m_tvertices.push_back(tvertex);
//                EdgeVertex ev{tg_edg[t],tg_idx[t]};
//                m_ehash[hashbucket_index(m_nbuckets, tg_edg[t])].push_back(ev);
//            }
//        }
//        
//        // triangulate:
//        unsigned char tcon_[12];
//        pos_ = 0;
//        for (int t = 0; t < cnt_; t++) {
//            if (tun_[t] != false) { // contour belongs to tunnel
//                for (int r = 0; r < c_sz[t]; r++) {
//                    uint index = -1;
//                    double dist = 4;
//                    uint ci = get_c(pos_+r,c_);
//                    for (int s = 0; s < 6; s++) {
//                        double uval = ip[ci][0] - hex_p[s][0];
//                        double vval = ip[ci][1] - hex_p[s][1];
//                        double wval = ip[ci][2] - hex_p[s][2];
//                        double val = uval*uval + vval*vval + wval*wval;
//                        if (dist > val) {
//                            index = s;
//                            dist = val;
//                        }
//                    }
//                    tcon_[ci] = (unsigned char)index;
//                }
//            }
//            pos_ += c_sz[t];
//        }
//        
//        // triangulate
//        pos_ = 0;
//        for (int t = 0; t < cnt_; t++) {
//            if (tun_[t]) {
//                for (int r = 0; r < c_sz[t]; r++) {
//                    // collect indices
//                    uint tid1 = get_c(pos_+r,c_);
//                    uint tid2 = get_c(pos_+((r+1)%c_sz[t]),c_);
//                    uint cid1 = tcon_[tid1];
//                    uint cid2 = tcon_[tid2];
//                    
//                    if (cid1 == cid2 ) {
//                        Triangle tr;
//                        tr.v[0] = vertices[tid1].g_idx;
//                        tr.v[1] = vertices[tid2].g_idx;
//                        tr.v[2] = tg_idx[cid1/2];
//                        //                        if (tr.v[1] > 1575043)
//                        //                            std::cout << "ERROR: wrong vertex index\n";
//                        m_ttriangles.push_back(tr);
//                    } else {
//                        // measure distance and compute
//                        // better triangulation: user shorter diag to divide
//                        double uval = ip[tid1][0] - tv[cid2/2][0];
//                        double vval = ip[tid1][1] - tv[cid2/2][1];
//                        double wval = ip[tid1][2] - tv[cid2/2][2];
//                        double val1 = uval*uval + vval*vval + wval*wval;
//                        // second diag
//                        uval = ip[tid2][0] - tv[cid1/2][0];
//                        vval = ip[tid2][1] - tv[cid1/2][1];
//                        wval = ip[tid2][2] - tv[cid1/2][2];
//                        double val2 = uval*uval + vval*vval + wval*wval;
//                        if (val1 < val2) {
//                            Triangle tr;
//                            tr.v[0] = vertices[tid1].g_idx;
//                            tr.v[1] = tg_idx[cid2/2];
//                            tr.v[2] = tg_idx[cid1/2];
//                            //                            if (tr.v[1] > 1575043)
//                            //                                std::cout << "ERROR: wrong vertex index\n";
//                            m_ttriangles.push_back(tr);
//                            tr.v[0] = vertices[tid1].g_idx;
//                            tr.v[1] = vertices[tid2].g_idx;
//                            tr.v[2] = tg_idx[cid2/2];
//                            //                            if (tr.v[1] > 1575043)
//                            //                                std::cout << "ERROR: wrong vertex index\n";
//                            m_ttriangles.push_back(tr);
//                        } else {
//                            Triangle tr;
//                            tr.v[0] = vertices[tid1].g_idx;
//                            tr.v[1] = vertices[tid2].g_idx;
//                            tr.v[2] = tg_idx[cid1/2];
//                            //                            if (tr.v[1] > 1575043)
//                            //                                std::cout << "ERROR: wrong vertex index\n";
//                            m_ttriangles.push_back(tr);
//                            tr.v[0] = vertices[tid2].g_idx;
//                            tr.v[1] = tg_idx[cid2/2];
//                            tr.v[2] = tg_idx[cid1/2];
//                            //                            if (tr.v[1] > 1575043)
//                            //                                std::cout << "ERROR: wrong vertex index\n";
//                            m_ttriangles.push_back(tr);
//                        }
//                    }
//                }
//            }
//            pos_ += c_sz[t];
//        }
//        // if there is a unique contour, then add triangle at the midpoint
//        if (cnt_ == 1) {
//            m_ccases_tunnel -= 1;
//            m_ccases_12cont += 1;
//            //std::cout << "12 contour at " << m_iindex << ", " << m_jindex << ", " << m_kindex << std::endl;
//            Triangle tr;
//            tr.v[0] = tg_idx[0];
//            tr.v[1] = tg_idx[1];
//            tr.v[2] = tg_idx[2];
//            //            if (tr.v[1] > 1575043)
//            //                std::cout << "ERROR: wrong vertex index\n";
//            m_ttriangles.push_back(tr);
//        }
//    }
//    pos_ = 0;
//    for (int t = 0; t < cnt_; t++) {
//        if (tun_[t]==false) {
//            // triangulation
//            uint csz = (uint) c_sz[t];
//            switch (csz) {
//                case 3:
//                {
//                    m_ccases_3 += 1;
//                    Triangle tri;
//                    tri.v[0] = vertices[get_c(pos_,c_)].g_idx; tri.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri.v[2] = vertices[get_c(pos_+2,c_)].g_idx;
//                    m_ttriangles.push_back(tri);
//                }
//                    break;
//                case 4:
//                {
//                    m_ccases_4 += 1;
//                    Triangle tri41;
//                    Triangle tri42;
//                    double d1 = distance(ip[get_c(pos_,c_)], ip[get_c(pos_+2,c_)]);
//                    double d2 = distance(ip[get_c(pos_+1,c_)], ip[get_c(pos_+3,c_)]);
//                    if (d1 < d2) {
//                        tri41.v[0] = vertices[get_c(pos_,c_)].g_idx;
//                        tri41.v[1] = vertices[get_c(pos_+1,c_)].g_idx;
//                        tri41.v[2] = vertices[get_c(pos_+2,c_)].g_idx;
//                        tri42.v[0] = vertices[get_c(pos_,c_)].g_idx;
//                        tri42.v[1] = vertices[get_c(pos_+2,c_)].g_idx;
//                        tri42.v[2] = vertices[get_c(pos_+3,c_)].g_idx;
//                    } else {
//                        tri41.v[0] = vertices[get_c(pos_,c_)].g_idx;
//                        tri41.v[1] = vertices[get_c(pos_+1,c_)].g_idx;
//                        tri41.v[2] = vertices[get_c(pos_+3,c_)].g_idx;
//                        tri42.v[0] = vertices[get_c(pos_+1,c_)].g_idx;
//                        tri42.v[1] = vertices[get_c(pos_+2,c_)].g_idx;
//                        tri42.v[2] = vertices[get_c(pos_+3,c_)].g_idx;
//                    }
//                    m_ttriangles.push_back(tri41);
//                    m_ttriangles.push_back(tri42);
//                }
//                    break;
//                case 5:
//                {
//                    m_ccases_5 += 1;
//                    Triangle tri51;
//                    Triangle tri52;
//                    Triangle tri53;
//                    tri51.v[0] = vertices[get_c(pos_,c_)].g_idx; tri51.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri51.v[2] = vertices[get_c(pos_+2,c_)].g_idx;
//                    tri52.v[0] = vertices[get_c(pos_,c_)].g_idx; tri52.v[1] = vertices[get_c(pos_+2,c_)].g_idx; tri52.v[2] = vertices[get_c(pos_+3,c_)].g_idx;
//                    tri53.v[0] = vertices[get_c(pos_,c_)].g_idx; tri53.v[1] = vertices[get_c(pos_+3,c_)].g_idx; tri53.v[2] = vertices[get_c(pos_+4,c_)].g_idx;
//                    m_ttriangles.push_back(tri51);
//                    m_ttriangles.push_back(tri52);
//                    m_ttriangles.push_back(tri53);
//                }
//                    break;
//                case 6:
//                {
//                    // check if there are asymptotes
//                    if (tc == 2 && ic[0] != 2 && ic[1] != 2 && ic[2] != 2) {
//                        m_ccases_6a += 1;
//                        Point  cv;
//                        Normal cn;
//                        double u{0},v{0},w{0};
//                        if (ic[0] == 0) {
//                            // face 0,1 has not asymptotes, common coordinate is w
//                            u = ui[1][0];
//                            v = ui[2][0];
//                            w = vi[1][0];
//                        } else if (ic[1] == 0) {
//                            // face 2,3 have no asymptotes, common coordinate is v
//                            u = ui[0][0];
//                            v = vi[0][0];
//                            w = vi[2][0];
//                        } else {
//                            // face 4,5 have no asymptotes, common coordinate is u
//                            u = ui[0][0];
//                            v = vi[0][0];
//                            w = vi[1][0];
//                        }
//                        // compute vertex: trilinear interpolate
//                        cv = (1-w)*((1-v)*((1-u)*p[0]+u*p[1]) + v*((1-u)*p[2]+u*p[3])) + w*((1-v)*((1-u)*p[4]+u*p[5]) + v*((1-u)*p[6]+u*p[7]));
//                        cn = (1-w)*((1-v)*((1-u)*n[0]+u*n[1]) + v*((1-u)*n[2]+u*n[3])) + w*((1-v)*((1-u)*n[4]+u*n[5]) + v*((1-u)*n[6]+u*n[7]));
//                        // store vertex in list
//                        // compute unique vertex id
//                        int g_edg = int(m_cell_shift_factor*m_ugrid.global_index(i, j, k) + 3); // this is the first interior vertex
//                        int g_idx = set_tvertex(g_edg);
//                        if (g_idx == -1) {
//                            TVertex tvertex;
//                            tvertex.p = cv;
//                            tvertex.n = cn;
//                            g_idx = (int)m_tvertices.size();
//                            tvertex.g_idx = g_idx;
//                            tvertex.g_edg = g_edg;
//                            m_tvertices.push_back(tvertex);
//                            EdgeVertex ev{g_edg,g_idx};
//                            m_ehash[hashbucket_index(m_nbuckets, g_edg)].push_back(ev);
//                        }
//                        // triangulate
//                        Triangle tri61;
//                        Triangle tri62;
//                        Triangle tri63;
//                        Triangle tri64;
//                        Triangle tri65;
//                        Triangle tri66;
//                        tri61.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri61.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri61.v[2] = g_idx;
//                        tri62.v[0] = vertices[get_c(pos_+1,c_)].g_idx; tri62.v[1] = vertices[get_c(pos_+2,c_)].g_idx; tri62.v[2] = g_idx;
//                        tri63.v[0] = vertices[get_c(pos_+2,c_)].g_idx; tri63.v[1] = vertices[get_c(pos_+3,c_)].g_idx; tri63.v[2] = g_idx;
//                        tri64.v[0] = vertices[get_c(pos_+3,c_)].g_idx; tri64.v[1] = vertices[get_c(pos_+4,c_)].g_idx; tri64.v[2] = g_idx;
//                        tri65.v[0] = vertices[get_c(pos_+4,c_)].g_idx; tri65.v[1] = vertices[get_c(pos_+5,c_)].g_idx; tri65.v[2] = g_idx;
//                        tri66.v[0] = vertices[get_c(pos_+5,c_)].g_idx; tri66.v[1] = vertices[get_c(pos_,c_)].g_idx;   tri66.v[2] = g_idx;
//                        m_ttriangles.push_back(tri61); m_ttriangles.push_back(tri62); m_ttriangles.push_back(tri63);
//                        m_ttriangles.push_back(tri64); m_ttriangles.push_back(tri65); m_ttriangles.push_back(tri66);
//                    } else {
//                        m_ccases_6 += 1;
//                        Triangle tri61;
//                        Triangle tri62;
//                        Triangle tri63;
//                        Triangle tri64;
//                        tri61.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri61.v[1] = vertices[get_c(pos_+4,c_)].g_idx; tri61.v[2] = vertices[get_c(pos_+5,c_)].g_idx;
//                        tri62.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri62.v[1] = vertices[get_c(pos_+2,c_)].g_idx; tri62.v[2] = vertices[get_c(pos_+4,c_)].g_idx;
//                        tri63.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri63.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri63.v[2] = vertices[get_c(pos_+2,c_)].g_idx;
//                        tri64.v[0] = vertices[get_c(pos_+2,c_)].g_idx; tri64.v[1] = vertices[get_c(pos_+3,c_)].g_idx; tri64.v[2] = vertices[get_c(pos_+4,c_)].g_idx;
//                        m_ttriangles.push_back(tri61);
//                        m_ttriangles.push_back(tri62);
//                        m_ttriangles.push_back(tri63);
//                        m_ttriangles.push_back(tri64);
//                    }
//                }
//                    break;
//                case 7:
//                {
//                    m_ccases_7 += 1;
//                    Point  cv;
//                    Normal cn;
//                    double u{0},v{0},w{0};
//                    if (ic[0] == 0) {
//                        // face 0,1 has not asymptotes, common coordinate is w
//                        u = ui[1][0];
//                        v = ui[2][0];
//                        w = vi[1][0];
//                    } else if (ic[1] == 0) {
//                        // face 2,3 have no asymptotes, common coordinate is v
//                        u = ui[0][0];
//                        v = vi[0][0];
//                        w = vi[2][0];
//                    } else {
//                        // face 4,5 have no asymptotes, common coordinate is u
//                        u = ui[0][0];
//                        v = vi[0][0];
//                        w = vi[1][0];
//                    }
//                    // compute vertex: trilinear interpolate
//                    cv = (1-w)*((1-v)*((1-u)*p[0]+u*p[1]) + v*((1-u)*p[2]+u*p[3])) + w*((1-v)*((1-u)*p[4]+u*p[5]) + v*((1-u)*p[6]+u*p[7]));
//                    cn = (1-w)*((1-v)*((1-u)*n[0]+u*n[1]) + v*((1-u)*n[2]+u*n[3])) + w*((1-v)*((1-u)*n[4]+u*n[5]) + v*((1-u)*n[6]+u*n[7]));
//                    // store vertex in list
//                    // compute unique vertex id
//                    int g_edg = int(m_cell_shift_factor*m_ugrid.global_index(i, j, k) + 3); // this is the first interior vertex
//                    int g_idx = set_tvertex(g_edg);
//                    if (g_idx == -1) {
//                        TVertex tvertex;
//                        tvertex.p = cv;
//                        tvertex.n = cn;
//                        g_idx = (int)m_tvertices.size();
//                        tvertex.g_idx = g_idx;
//                        tvertex.g_edg = g_edg;
//                        m_tvertices.push_back(tvertex);
//                        EdgeVertex ev{g_edg,g_idx};
//                        m_ehash[hashbucket_index(m_nbuckets, g_edg)].push_back(ev);
//                    }
//                    // triangulate
//                    Triangle tri71;
//                    Triangle tri72;
//                    Triangle tri73;
//                    Triangle tri74;
//                    Triangle tri75;
//                    Triangle tri76;
//                    Triangle tri77;
//                    tri71.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri71.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri71.v[2] = g_idx;
//                    tri72.v[0] = vertices[get_c(pos_+1,c_)].g_idx; tri72.v[1] = vertices[get_c(pos_+2,c_)].g_idx; tri72.v[2] = g_idx;
//                    tri73.v[0] = vertices[get_c(pos_+2,c_)].g_idx; tri73.v[1] = vertices[get_c(pos_+3,c_)].g_idx; tri73.v[2] = g_idx;
//                    tri74.v[0] = vertices[get_c(pos_+3,c_)].g_idx; tri74.v[1] = vertices[get_c(pos_+4,c_)].g_idx; tri74.v[2] = g_idx;
//                    tri75.v[0] = vertices[get_c(pos_+4,c_)].g_idx; tri75.v[1] = vertices[get_c(pos_+5,c_)].g_idx; tri75.v[2] = g_idx;
//                    tri76.v[0] = vertices[get_c(pos_+5,c_)].g_idx; tri76.v[1] = vertices[get_c(pos_+6,c_)].g_idx; tri76.v[2] = g_idx;
//                    tri77.v[0] = vertices[get_c(pos_+6,c_)].g_idx; tri77.v[1] = vertices[get_c(pos_,c_)].g_idx;   tri77.v[2] = g_idx;
//                    m_ttriangles.push_back(tri71); m_ttriangles.push_back(tri72); m_ttriangles.push_back(tri73);
//                    m_ttriangles.push_back(tri74); m_ttriangles.push_back(tri75); m_ttriangles.push_back(tri76);
//                    m_ttriangles.push_back(tri77);
//                }
//                    break;
//                case 8:
//                {
//                    m_ccases_8 += 1;
//                    // collect u,v,w, there are one which has two values, take always mean
//                    double u{0};
//                    double v{0};
//                    double w{0};
//                    if (tc == 4) {
//                        if (ic[0] == 2) {
//                            // face 0,1 has not asymptotes, common coordinate is w
//                            u = ui[1][0];
//                            v = ui[2][0];
//                            w = vi[1][0];
//                        } else if (ic[1] == 2) {
//                            // face 2,3 have no asymptotes, common coordinate is v
//                            u = ui[0][0];
//                            v = vi[0][0];
//                            w = vi[2][0];
//                        } else {
//                            // face 4,5 have no asymptotes, common coordinate is u
//                            u = ui[0][0];
//                            v = vi[0][0];
//                            w = vi[1][0];
//                        }
//                    } else {
//                        u = (ui[0][0] + ui[1][0]) / 2.;
//                        v = (vi[0][0] + ui[2][0]) / 2.;
//                        w = (vi[1][0] + vi[2][0]) / 2.;
//                    }
//                    // compute vertex: trilinear interpolate
//                    Point  cv = (1-w)*((1-v)*((1-u)*p[0]+u*p[1]) + v*((1-u)*p[2]+u*p[3])) + w*((1-v)*((1-u)*p[4]+u*p[5]) + v*((1-u)*p[6]+u*p[7]));
//                    Normal cn = (1-w)*((1-v)*((1-u)*n[0]+u*n[1]) + v*((1-u)*n[2]+u*n[3])) + w*((1-v)*((1-u)*n[4]+u*n[5]) + v*((1-u)*n[6]+u*n[7]));
//                    // store vertex in list
//                    // compute unique vertex id
//                    int g_edg = int(m_cell_shift_factor*m_ugrid.global_index(i, j, k) + 3); // this is the first interior vertex
//                    int g_idx = set_tvertex(g_edg);
//                    if (g_idx == -1) {
//                        TVertex tvertex;
//                        tvertex.p = cv;
//                        tvertex.n = cn;
//                        g_idx = (int)m_tvertices.size();
//                        tvertex.g_idx = g_idx;
//                        tvertex.g_edg = g_edg;
//                        m_tvertices.push_back(tvertex);
//                        EdgeVertex ev{g_edg,g_idx};
//                        m_ehash[hashbucket_index(m_nbuckets, g_edg)].push_back(ev);
//                    }
//                    Triangle tri81;
//                    Triangle tri82;
//                    Triangle tri83;
//                    Triangle tri84;
//                    Triangle tri85;
//                    Triangle tri86;
//                    Triangle tri87;
//                    Triangle tri88;
//                    tri81.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri81.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri81.v[2] = g_idx;
//                    tri82.v[0] = vertices[get_c(pos_+1,c_)].g_idx; tri82.v[1] = vertices[get_c(pos_+2,c_)].g_idx; tri82.v[2] = g_idx;
//                    tri83.v[0] = vertices[get_c(pos_+2,c_)].g_idx; tri83.v[1] = vertices[get_c(pos_+3,c_)].g_idx; tri83.v[2] = g_idx;
//                    tri84.v[0] = vertices[get_c(pos_+3,c_)].g_idx; tri84.v[1] = vertices[get_c(pos_+4,c_)].g_idx; tri84.v[2] = g_idx;
//                    tri85.v[0] = vertices[get_c(pos_+4,c_)].g_idx; tri85.v[1] = vertices[get_c(pos_+5,c_)].g_idx; tri85.v[2] = g_idx;
//                    tri86.v[0] = vertices[get_c(pos_+5,c_)].g_idx; tri86.v[1] = vertices[get_c(pos_+6,c_)].g_idx; tri86.v[2] = g_idx;
//                    tri87.v[0] = vertices[get_c(pos_+6,c_)].g_idx; tri87.v[1] = vertices[get_c(pos_+7,c_)].g_idx; tri87.v[2] = g_idx;
//                    tri88.v[0] = vertices[get_c(pos_+7,c_)].g_idx; tri88.v[1] = vertices[get_c(pos_,c_)].g_idx; tri88.v[2] = g_idx;
//                    m_ttriangles.push_back(tri81);
//                    m_ttriangles.push_back(tri82);
//                    m_ttriangles.push_back(tri83);
//                    m_ttriangles.push_back(tri84);
//                    m_ttriangles.push_back(tri85);
//                    m_ttriangles.push_back(tri86);
//                    m_ttriangles.push_back(tri87);
//                    m_ttriangles.push_back(tri88);
//                    //                    std::cout << "case 8: ";
//                    //                    for (int c8 = 0; c8 < 8; c8++) {
//                    //                        std::cout << (int)get_c(pos_+c8,c_) << ", ";
//                    //                    }
//                    //                    std::cout << std::endl;
//                }
//                    break;
//                case 9:
//                {
//                    m_ccases_9 += 1;
//                    // collect u,v,w,
//                    double u = (ui[0][0] + ui[0][1] + ui[1][0] + ui[1][1]) / (4-ic[2]);
//                    double v = (vi[0][0] + vi[0][1] + ui[2][0] + ui[2][1]) / (4-ic[1]);
//                    double w = (vi[1][0] + vi[1][1] + vi[2][0] + vi[2][1]) / (4-ic[0]);
//                    // compute vertex: trilinear interpolate
//                    Point  cv = (1-w)*((1-v)*((1-u)*p[0]+u*p[1]) + v*((1-u)*p[2]+u*p[3])) + w*((1-v)*((1-u)*p[4]+u*p[5]) + v*((1-u)*p[6]+u*p[7]));
//                    Normal cn = (1-w)*((1-v)*((1-u)*n[0]+u*n[1]) + v*((1-u)*n[2]+u*n[3])) + w*((1-v)*((1-u)*n[4]+u*n[5]) + v*((1-u)*n[6]+u*n[7]));
//                    // store vertex in list
//                    // compute unique vertex id
//                    int g_edg = int(m_cell_shift_factor*m_ugrid.global_index(i, j, k) + 3); // this is the first interior vertex
//                    int g_idx = set_tvertex(g_edg);
//                    if (g_idx == -1) {
//                        TVertex tvertex;
//                        tvertex.p = cv;
//                        tvertex.n = cn;
//                        g_idx = (uint)m_tvertices.size();
//                        tvertex.g_idx = g_idx;
//                        tvertex.g_edg = g_edg;
//                        m_tvertices.push_back(tvertex);
//                        EdgeVertex ev{g_edg,g_idx};
//                        m_ehash[hashbucket_index(m_nbuckets, g_edg)].push_back(ev);
//                    }
//                    Triangle tri91;
//                    Triangle tri92;
//                    Triangle tri93;
//                    Triangle tri94;
//                    Triangle tri95;
//                    Triangle tri96;
//                    Triangle tri97;
//                    Triangle tri98;
//                    Triangle tri99;
//                    tri91.v[0] = vertices[get_c(pos_,c_)].g_idx;   tri91.v[1] = vertices[get_c(pos_+1,c_)].g_idx; tri91.v[2] = g_idx;
//                    tri92.v[0] = vertices[get_c(pos_+1,c_)].g_idx; tri92.v[1] = vertices[get_c(pos_+2,c_)].g_idx; tri92.v[2] = g_idx;
//                    tri93.v[0] = vertices[get_c(pos_+2,c_)].g_idx; tri93.v[1] = vertices[get_c(pos_+3,c_)].g_idx; tri93.v[2] = g_idx;
//                    tri94.v[0] = vertices[get_c(pos_+3,c_)].g_idx; tri94.v[1] = vertices[get_c(pos_+4,c_)].g_idx; tri94.v[2] = g_idx;
//                    tri95.v[0] = vertices[get_c(pos_+4,c_)].g_idx; tri95.v[1] = vertices[get_c(pos_+5,c_)].g_idx; tri95.v[2] = g_idx;
//                    tri96.v[0] = vertices[get_c(pos_+5,c_)].g_idx; tri96.v[1] = vertices[get_c(pos_+6,c_)].g_idx; tri96.v[2] = g_idx;
//                    tri97.v[0] = vertices[get_c(pos_+6,c_)].g_idx; tri97.v[1] = vertices[get_c(pos_+7,c_)].g_idx; tri97.v[2] = g_idx;
//                    tri98.v[0] = vertices[get_c(pos_+7,c_)].g_idx; tri98.v[1] = vertices[get_c(pos_+8,c_)].g_idx; tri98.v[2] = g_idx;
//                    tri99.v[0] = vertices[get_c(pos_+8,c_)].g_idx; tri99.v[1] = vertices[get_c(pos_,c_)].g_idx; tri99.v[2] = g_idx;
//                    m_ttriangles.push_back(tri91);
//                    m_ttriangles.push_back(tri92);
//                    m_ttriangles.push_back(tri93);
//                    m_ttriangles.push_back(tri94);
//                    m_ttriangles.push_back(tri95);
//                    m_ttriangles.push_back(tri96);
//                    m_ttriangles.push_back(tri97);
//                    m_ttriangles.push_back(tri98);
//                    m_ttriangles.push_back(tri99);
//                    //                    std::cout << "case 9: ";
//                    //                    for (int c9 = 0; c9 < 9; c9++) {
//                    //                        std::cout << (int)get_c(pos_+c9,c_) << ", ";
//                    //                    }
//                    //                    std::cout << std::endl;
//                }
//                    break;
//                default:
//                    std::cout << "WARNING: " << c_sz[t] << " vertices in contour " << t << " for i_case: " << i_case << std::endl;
//                    break;
//            }
//        }
//        pos_ += c_sz[t];
//    }
//} // void p_slice()


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
                m_bbox[index] = Point{double(k),double(j),double(i)};
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
    // initialize scalar fields
    m_scalars.clear();
    m_normals.clear();
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
        n = Normal{0,0,0};
        return false;
    }
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





////*******************************************************************************************************************************************
////  IMPLEMENTATION function for test purposes
////*******************************************************************************************************************************************
//// interpolate the vector data at input position
//void
//tmc::MarchingCubes::processOneCell(std::vector<float>& bbox_v,std::vector<float>& bbox_n,std::vector<float>& bbox_c,std::vector<uint>& bbox_e,
//                                                      std::vector<std::vector<float>>& v,  std::vector<std::vector<float>>& n,
//                                                      std::vector<std::vector<float>>& c,  std::vector<std::vector<uint>>&  e,
//                                                      std::vector<std::vector<float>>& tv, std::vector<std::vector<float>>& tn,
//                                                      std::vector<std::vector<float>>& tc, std::vector<std::vector<uint>>&  te,
//                                                      std::vector<float>& points,std::vector<float>& pnorms)
//{
//    
//    double xmax = 1;
//    double ymax = 1;
//    double zmax = 0.3 / 0.4;
//    // set bounding box
//    bbox_v.resize(3*8);
//    int pos = 0;
//    for (int k = 0; k < 2; k++) {
//        for (int j = 0; j < 2; j++) {
//            for (int i = 0; i < 2; i++) {
//                bbox_v[3*pos + 0] = (float)i*xmax;
//                bbox_v[3*pos + 1] = (float)j*ymax;
//                bbox_v[3*pos + 2] = (float)k*zmax;
//                pos++;
//            }
//        }
//    }
//    // compute bbox normals at 8 vertices
//    bbox_n.resize(3*8);
//    pos = 0;
//    for (int k = 0; k < 2; k++) {
//        for (int j = 0; j < 2; j++) {
//            for (int i = 0; i < 2; i++) {
//                float vn[3];
//                vn[0] = (i == 0) ? -xmax : xmax;
//                vn[1] = (j == 0) ? -ymax : ymax;
//                vn[2] = (k == 0) ? -ymax : zmax;
//                //float sz = std::sqrt(3);
//                float sz = std::sqrt(vn[0]*vn[0] + vn[1]*vn[1] + vn[2]*vn[2]);
//                bbox_n[3*pos + 0] = vn[0]/sz;
//                bbox_n[3*pos + 1] = vn[1]/sz;
//                bbox_n[3*pos + 2] = vn[2]/sz;
//                pos++;
//            }
//        }
//    }
//    // compute colors
//    bbox_c.resize(3*8);
//    pos = 0;
//    float bbc[8][3] = {{0.1,0.1,0.1},{1,0,0},{0,1,0},{0.8,0.8,0.8},{0,0,1},{0.8,0.8,0.8},{0.8,0.8,0.8},{0.9,0.9,0.9}};
//    for (int k = 0; k < 2; k++) {
//        for (int j = 0; j < 2; j++) {
//            for (int i = 0; i < 2; i++) {
//                bbox_c[3*pos + 0] = bbc[pos][0];
//                bbox_c[3*pos + 1] = bbc[pos][1];
//                bbox_c[3*pos + 2] = bbc[pos][2];
//                pos++;
//            }
//        }
//    }
//    // compute lines segemests
//    bbox_e.resize(2*12); // there are 12 edges
//    pos = 0;
//    bbox_e[pos++]  = 0; bbox_e[pos++] = 1;
//    bbox_e[pos++]  = 1; bbox_e[pos++] = 3;
//    bbox_e[pos++]  = 3; bbox_e[pos++] = 2;
//    bbox_e[pos++]  = 2; bbox_e[pos++] = 0;
//    
//    bbox_e[pos++]  = 4; bbox_e[pos++] = 5;
//    bbox_e[pos++]  = 5; bbox_e[pos++] = 7;
//    bbox_e[pos++]  = 7; bbox_e[pos++] = 6;
//    bbox_e[pos++]  = 6; bbox_e[pos++] = 4;
//    
//    bbox_e[pos++]  = 0; bbox_e[pos++] = 4;
//    bbox_e[pos++]  = 1; bbox_e[pos++] = 5;
//    bbox_e[pos++]  = 3; bbox_e[pos++] = 7;
//    bbox_e[pos++]  = 2; bbox_e[pos++] = 6;
//    
//    // create scalar field
//    int mc = 105;
//    std::array<double,8> F{};
//    if (mc == 30) {
//        // tunnel, the first has even two tunnels
//        // the first example has 5 tunnels
//        // F = std::array<double, 8>{-3.37811990337124,0.473258332744286,2.54344310345736,7.87658724379480,4.38700713005133,-1.49950251870885,-4.21025867362045,-1.00233824192217};
//        // f = {-16.2669614745437,4.54362949424187,3.73698472065357,1.30355753541773,3.44918920037055,-2.37395197855922,-2.22401498175641,-3.43725135840084};
//        //f = {-0.484294567734286	0.676624087742047	1.59797165904310	1.97597515037871	0.319095109444789	-0.237879783226055	-0.703236632667890	-0.376471662697547};
//        //f ={-5.19491074997163,12.1411370147787,4.51687508029188,23.7520169338784,18.1237823075420,-11.0794063241597,-12.0776863106043,-1.00933986344768};
//        // case 30 fancy but without tunnel
//        // this case contains contours with 8 vertices (8 edges)
//        F = std::array<double, 8>{-11.2776683522162,4.34264464125961,2.20398713449477,4.94242873504933,5.56284709122434,-2.11044784460422,-10.3112450831516,-3.95693011723305};
//        //f = {-28.0660802670280,11.5348985096636,0.661165156427608,9.24807607673638,3.37369161205236,-12.4710878273712,-16.1585752223895,-27.8037975287738};
//    } else if (mc == 22) {
//        // case 22, one surface with 9 edges
//        F = std::array<double,8>{-15.6504952739285,2.90290077342601,24.5454566157887,-24.5274127623786,21.6741877710053,-4.49696327433901,-19.7891575872492,-15.5588482753161};
//        //f = {-0.166648729499781,1.20496388280327,0.526942569080289,-0.655079098476782,1.37942900628002,-0.749151592823710,-0.451541598502498,-0.0848213779969326};
//        //f = {-8.07458279195668,12.6860684502642,16.4371270364453,-28.2831095283080,12.5333231294999,-29.4925739940957,-9.04464846136196,-21.0339626770278};
//        // case 22, other case with one surface  with 9 edges
//        //f = {-0.107652770180584,1.92479616171011,0.0102684482681349,-0.775910464711502,1.63560644130687,-0.869694705363510,-0.0854358455109103,-0.400782649098897};
//        // case 22 with tunnel
//        //F = std::array<double,8>{-0.870292207640089,1.16040917473114,1.10072040367266,-0.145954798223727,1.70706223544379,-0.623055131485066,-0.351952380892271,-0.514249539867053};
//        //f = {-8.50512808246324,2.25743809995132,1.15136722493571,-5.45696139315505,5.36214351039946,-5.72025480562082,-2.96855742891816,-7.83964259861431};
//        // case 22 with 3 tunnels
//        //F = std::array<double,8>{-3.42744283804455,0.621278122151001,4.48110777981235,-1.95551129669134,2.30448107596369,-1.04182240925489,-3.51087814405650,-6.44976786808517};
//    } else if (mc == 60) {
//        // case 60 with tunnel
//        //F = std::array<double,8>{-0.100000000000000,-6.11000000000000,2,10.2000000000000,10.8000000000000,1.80000000000000,-8.20000000000000,-0.180000000000000};
//        F = std::array<double,8>{9.9985934885536665,9.9998695572230147,9.9999045831713928,9.999316745478131,9.9986117521866866,9.9998754368055813,9.9999031760062458,9.9992041920402936};
//    } else if (mc == 105) {
//        // case 105 with tunnel
//        //F = std::array<double,8>{2.74742516087490,-3.39187542578189,-12.5297639669456,0.431517989649243,-6.92460546400188,2.52228314017858,14.6950568276448,-10.0732624062474};
//        //F = std::array<double,8>{9.62943583411392,-9.57067290101394,-15.7184912369771,1.63719021729391,-9.00821395927326,12.9838964358508,27.1294747489964,-11.3443003256802};
//        // example of a countour with intersects twice but itself, and
//        // therefore there is no tunnel, contour with 9 vertices/edges
//        //F = std::array<double,8> {4.69314856679690,-25.6666841753773,-19.3439361061026,11.2891663083649,-5.72871085708909,12.8485897893816,14.4616618309557,-3.61934839891487};
//        F = std::array<double,8> {893,1135,1115,971,1111,778,919,1046};
//        //f = {2.56647391270132,-7.87546704094998,-24.0314386830922,0.877608326864389,-27.8666241843413,21.9109258856636,14.6592692141074,-17.3567518307032};
//    } else if (mc == 24) {
//        F = std::array<double,8>{-7.70146936482581,-3.21868369245987,-5.44023748418735,15.6051950593180,12.7611835388515,-4.46952393442309,-11.7240576326183,-9.23038948829007};
//    } else if (mc == 25) {
//        F = std::array<double,8>{9.9998593195995547,9.9993381282115549,9.9979160205452544,9.9986053863704142,9.9999374908631235,9.999424800002032,9.9983922749132219,9.999579324965488};
//    }
//    
//    
//    // apply rotation
//    //    double rot1[] = {4,5,0,1,6,7,2,3};
//    //    double rot2[] = {1,3,0,2,5,7,4,6};
//    //    std::array<double, 8> rF(F);
//    //    for (int i = 0; i < 8; i++){
//    //        rF[i] = F[rot1[i]];
//    //    }
//    //    F = rF;
//    //    // rotate again
//    //    for (int i = 0; i < 8; i++) {
//    //        rF[i] = F[rot2[i]];
//    //    }
//    //    F = rF;
//    
//    // discretize cell and compute isosurfaces
//    // init uniform grid
//    std::array<Point,8> bb{Point{0,0,0},Point{1,0,0},Point{0,1,0},Point{1,1,0},Point{0,0,zmax},Point{1,0,zmax},Point{0,1,zmax},Point{1,1,zmax}};
//    const int grsz = 50;
//    m_ugrid.init(grsz,grsz,grsz,bb);
//    double du = 1. / (grsz-1.);
//    double dv = du;
//    double dw = du;
//    
//    double u = 0;
//    for (int i = 0; i < grsz; i++) {
//        double v = 0;
//        for (int j = 0; j < grsz; j++) {
//            double w = 0;
//            for (int k = 0; k < grsz; k++) {
//                const double val = m_ugrid.trilinear(u,v,w,F);
//                m_ugrid.scalar(i, j, k, val);
//                w += dw;
//            }
//            v += dv;
//        }
//        u += du;
//    }
//    
//    m_ugrid.gradient();
//    
//    // set grid for tcm
//    m_ugrid.init(2,2,2,bb);
//    std::bitset<3> t_index;
//    for (int i = 0; i < 2; i++) {
//        t_index.set(0,i);
//        for (int j = 0; j < 2; j++) {
//            t_index.set(1,j);
//            for (int k = 0; k < 2; k++) {
//                t_index.set(2,k);
//                m_ugrid.scalar(i, j, k, F[t_index.to_ullong()]);
//            }
//        }
//    }
//    m_ugrid.gradient();
//    
//    
//    //double i0 = 9.99946; // case 60
//    //double i0 = 9.9994608191478135; // case 145, i.e 25
//    double i0 = 1000.;
//    std::array<double, 8> sF = F;
//    std::sort(sF.begin(),sF.end());
//    double maxv = -std::numeric_limits<double>::max();
//    double minv = -maxv;
//    for (auto f : sF) {
//        if (f < i0 && f > maxv) maxv = f;
//        if (f > i0 && f < minv) minv = f;
//    }
//    
//    // compute isosurfaces
//    int nIso = 1;
//    // generate nIso colors
//    std::vector<std::vector<float>> colors(10);
//    double minC = 0;
//    double maxC = 213;
//    double dCol = (maxC - minC) / (nIso - 1.);
//    double colV = minC;
//    for (int col = 0; col < nIso; col++) {
//        quitte::mesh::Vector<float,4> val{(float)colV,1,1,1};
//        quitte::mesh::Vector<float,4> rgbC;
//        rgbC = quitte::utils::hsv2rgb(val);
//        colors[col].resize(3);
//        colors[col][0] = rgbC[0];
//        colors[col][1] = rgbC[1];
//        colors[col][2] = rgbC[2];
//        colV += dCol;
//    }
//    // set fields
//    v.resize(nIso);
//    n.resize(nIso);
//    c.resize(nIso);
//    e.resize(nIso);
//    tv.resize(nIso);
//    tn.resize(nIso);
//    tc.resize(nIso);
//    te.resize(nIso);
//    
//    double di0 = (minv - maxv) / (nIso+1.);
//    double isov = maxv + di0;
//    for (int s = 0 ; s < nIso; s++) {
//        if (s == 0)
//            isov = 1000;
//        iso_surface(isov);
//        const int nr_v = (int)m_vertices.size();
//        const int nr_t = (int)m_triangles.size();
//        v[s].resize(3*nr_v);
//        n[s].resize(3*nr_v);
//        c[s].resize(3*nr_v);
//        e[s].resize(3*nr_t);
//        std::vector<int> t_nrv(nr_v,0);
//        for (auto vt : m_vertices) {
//            auto id = vt->i;
//            if (id >= nr_v)
//                std::cout << "ERROR: id: " << id << ", nr_v: " << nr_v << std::endl;
//            v[s][3*id]     = (float)vt->p[0];
//            v[s][3*id + 1] = (float)vt->p[1];
//            v[s][3*id + 2] = (float)vt->p[2];
//            
//            n[s][3*id]     = -(float)vt->n[0];
//            n[s][3*id + 1] = -(float)vt->n[1];
//            n[s][3*id + 2] = -(float)vt->n[2];
//            
//            c[s][3*id]   = colors[s][0];
//            c[s][3*id+1] = colors[s][1];
//            c[s][3*id+2] = colors[s][2];
//            t_nrv[id] += 1;
//        }
//        for (auto p : t_nrv) {
//            if (p != 1)
//                std::cout << "ERROR: wrong nr. of vertices\n";
//        }
//        int count = 0;
//        for (auto t : m_triangles) {
//            e[s][3*count]     = t.v[0];
//            e[s][3*count + 1] = t.v[1];
//            e[s][3*count + 2] = t.v[2];
//            count++;
//        }
//        
//        // compute cell intersection
//        t_mc(isov);
//        const size_t nr_tv = m_tvertices.size();
//        const size_t nr_tt = m_ttriangles.size();
//        if (nr_tv > 0) {
//            tv[s].resize(3*nr_tv);
//            tn[s].resize(3*nr_tv);
//            tc[s].resize(3*nr_tv);
//            te[s].resize(3*nr_tt);
//            for (auto vt : m_tvertices) {
//                auto id = vt.g_idx;
//                tv[s][3*id]     = (float)vt.p[0];
//                tv[s][3*id + 1] = (float)vt.p[1];
//                tv[s][3*id + 2] = (float)vt.p[2];
//                
//                tn[s][3*id]     = -(float)vt.n[0];
//                tn[s][3*id + 1] = -(float)vt.n[1];
//                tn[s][3*id + 2] = -(float)vt.n[2];
//                
//                tc[s][3*id]   = colors[s][0];
//                tc[s][3*id+1] = colors[s][1];
//                tc[s][3*id+2] = colors[s][2];
//                
//            }
//            count = 0;
//            for (auto t : m_ttriangles) {
//                te[s][3*count]     = t.v[0];
//                te[s][3*count + 1] = t.v[1];
//                te[s][3*count + 2] = t.v[2];
//                count++;
//            }
//        }
//        
//        // check if there are points
//        size_t nr_pts = m_points.size();
//        if (nr_pts > 0) {
//            points.resize(3*nr_pts);
//            pnorms.resize(3*nr_pts);
//            int idx = 0;
//            for (int i = 0; i < (int) m_points.size(); i++) {
//                points[3*i]   = m_points[i][0];
//                points[3*i+1] = m_points[i][1];
//                points[3*i+2] = m_points[i][2];
//                pnorms[3*i]   = m_pnorms[i][0];
//                pnorms[3*i+1] = m_pnorms[i][1];
//                pnorms[3*i+2] = m_pnorms[i][2];
//                idx++;
//            }
//        }
//        // update isovalue
//        isov += di0;
//        
//        // clear all values
//        clear_all();
//    }
//    
//}


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



//void
//tmc::MarchingCubes::testData(int& nr_v,float** vertices,float** normals,int& nr_t,int** triangles, Mesh& mesh)
//{
//    std::cout << " ... create test data \n";
//    // set grid sizes
//    m_nx = 128;
//    m_ny = 128;
//    m_nz = 128;
//    
//    // set domain
//    std::array<Point, 8> bb;
//    bb[0] = Point{0,0,0};
//    bb[1] = Point{1,0,0};
//    bb[2] = Point{0,1,0};
//    bb[3] = Point{1,1,0};
//    bb[4] = Point{0,0,1};
//    bb[5] = Point{1,0,1};
//    bb[6] = Point{0,1,1};
//    bb[7] = Point{1,1,1};
//    
//    m_ugrid.init(m_nx, m_ny, m_nz,bb);
//    m_dx = m_ugrid.dx();
//    m_dy = m_ugrid.dy();
//    m_dz = m_ugrid.dz();
//    // the test data set
//    auto ft = [] (double x, double y, double z) {
//        double alpha = 1.6;
//        x = alpha * (2*x - 1);
//        y = alpha * (2*y - 1);
//        z = alpha * (2*z - 1);
//        //        double c1[]{-1,0,0};
//        //        double c2[]{ 1,0,0};
//        //        double c3[]{1.8,1.8,1.8};
//        //        double c4[]{-1.8,-1.8,-1.8};
//        //        double val1 = (x-c1[0])*(x-c1[0]) + (y-c1[1])*(y-c1[1]) + (z-c1[2])*(z-c1[2]);
//        //        double val2 = (x-c2[0])*(x-c2[0]) + (y-c2[1])*(y-c2[1]) + (z-c2[2])*(z-c2[2]);
//        //        double val3 = (x-c3[0])*(x-c3[0]) + (y-c3[1])*(y-c3[1]) + (z-c3[2])*(z-c3[2]);
//        //        double val4 = (x-c4[0])*(x-c4[0]) + (y-c4[1])*(y-c4[1]) + (z-c4[2])*(z-c4[2]);
//        //        double beta = 5;
//        //        return std::exp(-beta*val1) + std::exp(-beta*val2) + std::exp(-beta*val3) + std::exp(-beta*val4);
//        return 2*y*(y*y - 3*x*x)*(1-z*z) + (x*x+y*y)*(x*x+y*y) - (9*z*z -1) * (1-z*z);
//    };
//    
//    int m_size = m_nx * m_ny * m_nz;
//    int x_size = m_nx;
//    int y_size = m_ny;
//    int z_size = m_nz;
//    auto g_index = [x_size,y_size,z_size]  (const int i, const int j, const int k) { return (k*y_size*x_size+j*x_size+i); };
//    double z = bb[0][0];
//    for (int k = 0; k < z_size; k++) {
//        double y = bb[0][1];
//        for (int j = 0; j < y_size; j++) {
//            double x = bb[0][0];
//            for (int i = 0; i < x_size; i++) {
//                //const double val = ft(x,y,z);
//                m_ugrid.scalar(g_index(i,j,k), ft(x,y,z));
//                x = x + m_dx;
//            }
//            y = y + m_dy;
//        }
//        z = z + m_dz;
//    }
//    
//    
//    // compute gradient for normals
//    m_ugrid.gradient();
//    // invert gradient
//    //m_ugrid.invert_normals();
//    
//    // compute isosurface
//    std::cout << " ... computing isosurface\n";
//    const double i0 = 0.4;
//    t_mc(i0);
//    
//    
//    bool t_cm = true;
//    std::cout << " ... copy data" << std::endl;
//    timer.start();
//    if (t_cm) {
//        nr_v = (int)m_tvertices.size();
//        nr_t = (int)m_ttriangles.size();
//        
//        std::cout << "tot. nr. of triangles: " << nr_t << std::endl;
//        std::cout << "tot. nr. of vertices:  " << nr_v << std::endl;
//        
//        *vertices  = new float[4*nr_v];
//        *normals   = new float[4*nr_v];
//        *triangles = new int[4*nr_t];
//        std::vector<int> t_vt(nr_v,0);
//        // global vertex id in the triangle mesh was already set
//        int vcount = 0;
//        for (auto v : m_tvertices) {
//            auto id = v.g_idx;
//            (*vertices)[4*id]     = (float)v.p[0];
//            (*vertices)[4*id + 1] = (float)v.p[1];
//            (*vertices)[4*id + 2] = (float)v.p[2];
//            (*vertices)[4*id + 3] = 1.0f;
//            (*normals)[4*id]     = (float)v.n[0];
//            (*normals)[4*id + 1] = (float)v.n[1];
//            (*normals)[4*id + 2] = (float)v.n[2];
//            (*normals)[4*id + 3] = 0.0f;
//            vcount++;
//            //t_vt[id]++;
//        }
//        
//        int count = 0;
//        for (auto t : m_ttriangles) {
//            (*triangles)[4*count]     = t.v[0];
//            (*triangles)[4*count + 1] = t.v[1];
//            (*triangles)[4*count + 2] = t.v[2];
//            t_vt[t.v[0]] += 1;
//            t_vt[t.v[1]] += 1;
//            t_vt[t.v[2]] += 1;
//            count++;
//        }
//        for (auto v : t_vt) {
//            if (v == 0) {
//                std::cout << "EXTREM ERROR!\n";
//            }
//        }
//    } else {
//        nr_v = (int)m_vertices.size();
//        nr_t = (int)m_triangles.size();
//        *vertices  = new float[4*nr_v];
//        *normals   = new float[4*nr_v];
//        *triangles = new int[4*nr_t];
//        std::vector<int> t_vt(nr_v,0);
//        // global vertex id in the triangle mesh was already set
//        int vcount = 0;
//        for (auto v : m_vertices) {
//            auto id = v->i;
//            (*vertices)[4*id]     = (float)v->p[0];
//            (*vertices)[4*id + 1] = (float)v->p[1];
//            (*vertices)[4*id + 2] = (float)v->p[2];
//            (*vertices)[4*id + 3] = 1.0f;
//            (*normals)[4*id]     = (float)v->n[0];
//            (*normals)[4*id + 1] = (float)v->n[1];
//            (*normals)[4*id + 2] = (float)v->n[2];
//            (*normals)[4*id + 3] = 0.0f;
//            vcount++;
//            t_vt[id]++;
//        }
//        for (auto v : t_vt) {
//            if (v != 1) {
//                std::cout << "EXTREM ERROR!\n";
//            }
//        }
//        int count = 0;
//        for (auto t : m_triangles) {
//            (*triangles)[4*count]     = t.v[0];
//            (*triangles)[4*count + 1] = t.v[1];
//            (*triangles)[4*count + 2] = t.v[2];
//            count++;
//        }
//        
//    }
//    
//    // create a halfedge mesh
//    reconstruct(true, mesh);
//    // analize mesh data
//    topology(mesh);
//}