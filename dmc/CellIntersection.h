#pragma once

// C++ libs
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"
#include "UniformGrid.h"
#include "QuadrilateralHashTable.h"
#include "Vertices.h"
#include "MarchingCubesLookupTables.h"
#include "GaussElimination.h"


namespace p_mc {
    /// <summary>
    /// Computes the MC polygons given a cell, the iso-value and the scalar values
    /// at the cell vertices.
    /// </summary>
    struct CellIntersection {
        using uchar = unsigned char;
        using ushort = unsigned short;
        using uint = unsigned int;
        using UGrid = p_mc::UniformGrid;
        /// mc polygons
        char* r_pattern{ nullptr };
        std::shared_ptr<char> r_pattern_;
        /// list of ambiguous cases
        uchar* t_ambig{ nullptr };
        std::shared_ptr<uchar> t_ambig_;
        /// constructor
        __host__ CellIntersection() {}
        __host__ CellIntersection(const std::array<char, 4352>& p, const std::array<uchar, 256>& a)
        {
            // alloc and init r_pattern
            cudaMalloc(&r_pattern, 4352 * sizeof(char));
            cudaMemcpy(r_pattern, &p[0], 4352 * sizeof(char), cudaMemcpyHostToDevice);
            cudaMalloc(&t_ambig, 256 * sizeof(uchar));
            cudaMemcpy(t_ambig, &a[0], 256 * sizeof(uchar), cudaMemcpyHostToDevice);
            r_pattern_ = std::shared_ptr<char>(r_pattern, cudaFree);
            t_ambig_ = std::shared_ptr<uchar>(t_ambig, cudaFree);
        }
        /// destructor
        __host__ ~CellIntersection()
        {
            r_pattern = nullptr;
            t_ambig = nullptr;
            r_pattern_.reset();
            t_ambig_.reset();
        }
        /// set size of mc polygon
        __device__ void set_cnt_size(const int cnt, const int size, unsigned long long& c_) {
            // unset contour size
            c_ &= ~(0xF << 4 * cnt);
            c_ |= (size << 4 * cnt);
        }
        /// get size of mc polygon
        __device__ int get_cnt_size(const int cnt, unsigned long long& c_) {
            return static_cast<int>((c_ & (0xF << 4 * cnt)) >> 4 * cnt);
        }
        /// set corresponging edge
        __device__ void set_c(const int cnt, const int pos, const int val, unsigned long long& c_) {
            const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
            const uint c_sz = c_ & mask[cnt];
            const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
            c_ &= ~(((unsigned long long)0xF) << e);
            c_ |= (((unsigned long long)val) << e);
        }
        /// read edge from polygon
        __device__ int get_c(const int cnt, const int pos, unsigned long long c_) {
            const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
            const uint c_sz = (uint)(c_ & mask[cnt]);
            const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
            return static_cast<int>((c_ >> e) & 0xF);
        }
        /// intersect a cell with iso-surface for unambiguous cases
        __device__ int mc_polygon(const int i_case, unsigned long long& c_)
        {
            int cnt_{ 0 };
            const char* c_lt = &r_pattern[17 * i_case];
            unsigned char pos2 = c_lt[0] + 1;
            for (int c = 1; c <= c_lt[0]; c++) // loop over contours
            {
                // set polygon size
                set_cnt_size(cnt_, c_lt[c], c_);
                // for all this edges save the vertex
                for (int i = 0; i < c_lt[c]; i++)
                {
                    const uint e = c_lt[pos2++];
                    set_c(cnt_, i, e, c_);
                }
                cnt_++;
            }
            return cnt_;
        }
        /// intersect cell for ambiguous cases
        __device__ int mc_polygon(const float i0, const float F[8], unsigned long long& c_)
        {
            // compute oriented contours
            // 1. build segments
            // 2. connect segments
            // build up segments
            // set segments map
            unsigned char segm_[12] = { 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };
            auto set_segm = [](const int ei, const int eo, unsigned char segm_[12]) {
                segm_[ei] &= 0xF0;
                segm_[ei] |= ((unsigned char)eo) & 0xF;
                segm_[eo] &= 0xF;
                segm_[eo] |= ((unsigned char)ei) << 4;
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
            unsigned short e_face_[6]{ (ushort)291, (ushort)18277, (ushort)18696, (ushort)10859, (ushort)33719, (ushort)38305 };
            // code vertices at face
            unsigned short v_face_[6]{ (ushort)12576, (ushort)25717, (ushort)5380, (ushort)29538, (ushort)8292, (ushort)30001 };

            // reading edge from face
            auto get_face_e = [e_face_](const int f, const int e) { return ((e_face_[f] >> (4 * e)) & 0xF); };
            auto get_face_v = [v_face_](const int f, const int e) { return ((v_face_[f] >> (4 * e)) & 0xF); };
            // compute oriented segments using the isoline scheme at the faces
            auto asymptotic_decider = [](const float f0, const float f1, const float f2, const float f3) {
                return (f0 * f3 - f1 * f2) / (f0 + f3 - f1 - f2);
            };
            for (int f = 0; f < 6; f++) {
                // classify face
                unsigned int f_case{ 0 };
                const int v0 = get_face_v(f, 0);
                const int v1 = get_face_v(f, 1);
                const int v2 = get_face_v(f, 2);
                const int v3 = get_face_v(f, 3);
                const int e0 = get_face_e(f, 0);
                const int e1 = get_face_e(f, 1);
                const int e2 = get_face_e(f, 2);
                const int e3 = get_face_e(f, 3);
                const float f0 = F[v0];
                const float f1 = F[v1];
                const float f2 = F[v2];
                const float f3 = F[v3];
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
                    set_segm(e0, e3, segm_);
                    break;
                case 2:
                    set_segm(e1, e0, segm_);
                    break;
                case 3:
                    set_segm(e1, e3, segm_);
                    break;
                case 4:
                    set_segm(e3, e2, segm_);
                    break;
                case 5:
                    set_segm(e0, e2, segm_);
                    break;
                case 6:
                {
                    const float val = asymptotic_decider(f0, f1, f2, f3);
                    if (val >= i0) {
                        set_segm(e3, e0, segm_);
                        set_segm(e1, e2, segm_);
                    }
                    else if (val < i0) {
                        set_segm(e1, e0, segm_);
                        set_segm(e3, e2, segm_);
                    }
                }
                break;
                case 7:
                    set_segm(e1, e2, segm_);
                    break;
                case 8:
                    set_segm(e2, e1, segm_);
                    break;
                case 9:
                {
                    const float val = asymptotic_decider(f0, f1, f2, f3);
                    if (val >= i0) {
                        set_segm(e0, e1, segm_);
                        set_segm(e2, e3, segm_);
                    }
                    else if (val < i0) {
                        set_segm(e0, e3, segm_);
                        set_segm(e2, e1, segm_);
                    }
                }
                break;
                case 10:
                    set_segm(e2, e0, segm_);
                    break;
                case 11:
                    set_segm(e2, e3, segm_);
                    break;
                case 12:
                    set_segm(e3, e1, segm_);
                    break;
                case 13:
                    set_segm(e0, e1, segm_);
                    break;
                case 14:
                    set_segm(e3, e0, segm_);
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
            //unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            // in the 4 first bits store size of contours
            // connect oriented contours
            int cnt_{ 0 };
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

            // return number of countures
            return cnt_;
        }
        /// intersect cell, compute vertices, set edges
        __device__ void slice(const float i0, const int i_case, const int i_index, const int j_index, const int k_index,
            float u[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            auto e_glIndex = [](const int e, const int i_idx, const int j_idx, const int k_idx, UGrid& ugrid)
            {
                const unsigned long long gei_pattern_ = 670526590282893600ull;
                const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
                const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
                const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
                const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
                return (3 * ugrid.gl_index(i, j, k) + offs);
            };
            // edge table unsigned long long e_table = 240177437832960;
            // input e and offset, which direction the normal has to point
            auto get_vertex_pos = [](const int e, const int offset)
            {
                const unsigned long long e_talbe = 240177437832960;
                return (e_talbe >> (4 * e + 2 * offset)) & 3;
            };
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mc_polygon(i0, u, c_);
            }
            else {
                cnt_ = mc_polygon(i_case, c_);
            }
            const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
            ushort e_{ 0 };
            // compute normals at vertices
            float3 n[8];
            ugrid.gradient(n, u, i_index, j_index, k_index);
            for (int t = 0; t < cnt_; t++)
            {
                float3 ei = make_float3(0, 0, 0);
                const int cnt_size = get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = get_c(t, i, c_);
                    const int v0 = (l_edges_[e] & 0xF);
                    const int v1 = (l_edges_[e] >> 4) & 0xF;
                    const float l = (i0 - u[v0]) / (u[v1] - u[v0]);
                    const float3 vt0{ ugrid.x0 + (i_index + (v0 & 0x1)) * ugrid.dx, ugrid.y0 + (j_index + ((v0 & 0x2) >> 1)) * ugrid.dy, ugrid.z0 + (k_index + ((v0 & 0x4) >> 2)) * ugrid.dz };
                    const float3 vt1{ ugrid.x0 + (i_index + (v1 & 0x1)) * ugrid.dx, ugrid.y0 + (j_index + ((v1 & 0x2) >> 1)) * ugrid.dy, ugrid.z0 + (k_index + ((v1 & 0x4) >> 2)) * ugrid.dz };
                    ei.x += vt0.x + l * (vt1.x - vt0.x);
                    ei.y += vt0.y + l * (vt1.y - vt0.y);
                    ei.z += vt0.z + l * (vt1.z - vt0.z);
                    //set edge case to construct oriented quadrilateral
                    if (u[v0] < i0) e_ |= (1 << e);
                }
                // normalize
                ei.x /= cnt_size;
                ei.y /= cnt_size;
                ei.z /= cnt_size;
                //float3 err{ 0,0,0 };
                float3 ni = make_float3(0, 0, 0);
                levelsetIntersection(i0, i_index, j_index, k_index, ugrid, u, n, ei, ni);
                // save unique vertex
                const int v_addr = v_.addVertex(ei, ni);
                // for all this edges save the vertex
                for (int i = 0; i < get_cnt_size(t, c_); i++)
                {
                    const uint e = get_c(t, i, c_);
                    // compute unique edges id
                    const int e_glId = e_glIndex(e, i_index, j_index, k_index, ugrid);
                    const int pos = get_vertex_pos(e, (e_ >> e) & 1);
                    //printf("add edge\n");
                    ht_.addVertex(e_glId, pos, v_addr);
                }
            }
        }
        /// check if case is ambiguous
        __device__ bool isAmbiguous(const int i_case)
        {
            if (t_ambig[i_case] == MC_AMBIGUOUS) return true;
            else return false;
        }
        /// project point back to surface
        __device__ void projectPoint(const float i0, const int i, const int j, const int k, UGrid& ugrid, float f[8], float3 n[8], float3& p, float3& n_)
        {
            const int max_iter = 5;
            int iter = 1;
            //const float delta = 0.01;
            float3 grad;
            float u, v, w;
            u = (p.x - ugrid.x0 - i * ugrid.dx) / ugrid.dx;
            v = (p.y - ugrid.y0 - j * ugrid.dy) / ugrid.dy;
            w = (p.z - ugrid.z0 - k * ugrid.dz) / ugrid.dz;
            float ii = ugrid.trilinear(f, u, v, w);
            const float err0 = fabsf(i0 - ii);
            const float delta = 0.1 * err0;
            float err1 = err0;
            while (iter < max_iter && err1 > delta)
            {

                grad.x = (1 - w) * ((1 - v) * (f[1] - f[0]) + v * (f[3] - f[2])) + w * ((1 - v) * (f[5] - f[4]) + v * (f[7] - f[6]));
                grad.y = (1 - w) * ((1 - u) * (f[2] - f[0]) + u * (f[3] - f[1])) + w * ((1 - u) * (f[6] - f[4]) + u * (f[7] - f[5]));
                grad.z = (1 - v) * ((1 - u) * (f[4] - f[0]) + u * (f[5] - f[1])) + v * ((1 - u) * (f[6] - f[2]) + u * (f[7] - f[3]));
                float l = 0.5 * (i0 - ii) / (grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
                u = u + l * grad.x;
                v = v + l * grad.y;
                w = w + l * grad.z;

                // check if we are within the cell
                if ((u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1))
                {
                    err1 = err0 + 0.1;
                    break;
                }
                ii = ugrid.trilinear(f, u, v, w);
                err1 = fabsf(i0 - ii);
                iter++;
            }

            if (err1 < err0) { // there were some improvement in the approx.
                p.x = ugrid.x0 + (i + u) * ugrid.dx;
                p.y = ugrid.y0 + (j + v) * ugrid.dy;
                p.z = ugrid.z0 + (k + w) * ugrid.dz;
                n_ = ugrid.trilinear(n, u, v, w);
                const float factor = sqrtf(n_.x * n_.x + n_.y * n_.y + n_.z * n_.z);
                n_.x /= factor;
                n_.y /= factor;
                n_.z /= factor;
            }
            else
            {
                u = (p.x - ugrid.x0 - i * ugrid.dx) / ugrid.dx;
                v = (p.y - ugrid.y0 - j * ugrid.dy) / ugrid.dy;
                w = (p.z - ugrid.z0 - k * ugrid.dz) / ugrid.dz;
                n_ = ugrid.trilinear(n, u, v, w);
                const float factor = sqrtf(n_.x * n_.x + n_.y * n_.y + n_.z * n_.z);
                n_.x /= factor;
                n_.y /= factor;
                n_.z /= factor;
            }
        }
        /// Starting at an input point move until level set intersection
        __device__ void levelsetIntersection(const float i0, const int i, const int j, const int k, UGrid& ugrid, float f[8], float3 n[8], float3& p, float3& n_)
        {
            const int max_iter = 20;
            const float step = 0.05;
            float3 grad;
            float u1, v1, w1;
            u1 = (p.x - ugrid.x0 - i * ugrid.dx) / ugrid.dx;
            v1 = (p.y - ugrid.y0 - j * ugrid.dy) / ugrid.dy;
            w1 = (p.z - ugrid.z0 - k * ugrid.dz) / ugrid.dz;
            float val1 = ugrid.trilinear(f, u1, v1, w1);
            for (int iter = 0; iter < max_iter; iter++)
            {

                grad.x = (1 - w1) * ((1 - v1) * (f[1] - f[0]) + v1 * (f[3] - f[2])) + w1 * ((1 - v1) * (f[5] - f[4]) + v1 * (f[7] - f[6]));
                grad.y = (1 - w1) * ((1 - u1) * (f[2] - f[0]) + u1 * (f[3] - f[1])) + w1 * ((1 - u1) * (f[6] - f[4]) + u1 * (f[7] - f[5]));
                grad.z = (1 - v1) * ((1 - u1) * (f[4] - f[0]) + u1 * (f[5] - f[1])) + v1 * ((1 - u1) * (f[6] - f[2]) + u1 * (f[7] - f[3]));
                if (val1 > i0)
                {
                    // invert gradient
                    grad.x = -grad.x;
                    grad.y = -grad.y;
                    grad.z = -grad.z;
                }
                // normalize
                const float sz = sqrt(grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
                grad.x /= sz;
                grad.y /= sz;
                grad.z /= sz;
                float u2, v2, w2;
                u2 = u1 + step * grad.x;
                v2 = v1 + step * grad.y;
                w2 = w1 + step * grad.z;

                // check if we are within the cell
                if ((u2 < 0 || u2 > 1) || (v2 < 0 || v2 > 1) || (w2 < 0 || w2 > 1))
                {
                    u1 = (p.x - ugrid.x0 - i * ugrid.dx) / ugrid.dx;
                    v1 = (p.y - ugrid.y0 - j * ugrid.dy) / ugrid.dy;
                    w1 = (p.z - ugrid.z0 - k * ugrid.dz) / ugrid.dz;
                    break;
                }
                const float val2 = ugrid.trilinear(f, u2, v2, w2);
                if ((val1 < i0 && i0 < val2) || (val2 < i0 && i0 < val1))
                {
                    const float e = (i0 - val1) / (val2 - val1);
                    u1 = u1 + e * (u2 - u1);
                    v1 = v1 + e * (v2 - v1);
                    w1 = w1 + e * (w2 - w1);
                    break;
                }
                // update
                u1 = u2;
                v1 = v2;
                w1 = w2;
                val1 = val2;
            }

            p.x = ugrid.x0 + (i + u1) * ugrid.dx;
            p.y = ugrid.y0 + (j + v1) * ugrid.dy;
            p.z = ugrid.z0 + (k + w1) * ugrid.dz;
            n_ = ugrid.trilinear(n, u1, v1, w1);
            const float factor = sqrt(n_.x * n_.x + n_.y * n_.y + n_.z * n_.z);
            n_.x /= factor;
            n_.y /= factor;
            n_.z /= factor;
        }
        /// Revised implementation of the slice method, including Quadric based repositioning of vertex
        __device__ void sliceP(const float i0, const int i_case, const int i_index, const int j_index, const int k_index,
            float f[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            auto e_glIndex = [](const int e, const int i_idx, const int j_idx, const int k_idx, UGrid& ugrid)
            {
                const unsigned long long gei_pattern_ = 670526590282893600ull;
                const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
                const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
                const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
                const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
                return (3 * ugrid.gl_index(i, j, k) + offs);
            };
            // edge table unsigned long long e_table = 240177437832960;
            // input e and offset, which direction the normal has to point
            auto get_vertex_pos = [](const int e, const int offset)
            {
                const unsigned long long e_talbe = 240177437832960;
                return (e_talbe >> (4 * e + 2 * offset)) & 3;
            };
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mc_polygon(i0, f, c_);
            }
            else {
                cnt_ = mc_polygon(i_case, c_);
            }
            const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
            ushort e_{ 0 };
            // compute normals at vertices
            float3 n[8];
            ugrid.gradient(n, f, i_index, j_index, k_index);
            for (int t = 0; t < cnt_; t++)
            {
                float3 ei = make_float3(0, 0, 0);
                const int cnt_size = get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = get_c(t, i, c_);
                    getLocalCoordinates(l_edges_[e], e, f, i0, ei);
                    //set edge case to construct oriented quadrilateral
                    if (f[(l_edges_[e] & 0xF)] < i0) e_ |= (1 << e);
                }
                // normalize
                ei.x /= cnt_size;
                ei.y /= cnt_size;
                ei.z /= cnt_size;

                movePointToSurface(f, i0, ei);
                // compute normal at mesh vertex
                float3 ni = trilinear(n, ei.x, ei.y, ei.z);
                normalize(ni);

                // compute point in space (voxel grid)
                ei.x = ugrid.x0 + (i_index + ei.x) * ugrid.dx;
                ei.y = ugrid.y0 + (j_index + ei.y) * ugrid.dy;
                ei.z = ugrid.z0 + (k_index + ei.z) * ugrid.dz;

                const int v_addr = v_.addVertex(ei, ni);

                // for all this edges save the vertex
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = get_c(t, i, c_);
                    int color{ -1 };
                    if (e == 0) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2);
                    if (e == 3) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2) + 1;
                    if (e == 8) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2) + 2;
                    // compute unique edges id
                    const int e_glId = e_glIndex(e, i_index, j_index, k_index, ugrid);
                    const int pos = get_vertex_pos(e, (e_ >> e) & 1);
                    if (color > -1) ht_.addVertex(e_glId, pos, v_addr, color);
                    else ht_.addVertex(e_glId, pos, v_addr);
                }
            }
        }
        /// Revised implementation of the slice method, including Quadric based repositioning of vertex
        __device__ void sliceQ(const float i0, const int i_case, const int i_index, const int j_index, const int k_index,
            float f[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            auto e_glIndex = [](const int e, const int i_idx, const int j_idx, const int k_idx, UGrid& ugrid)
            {
                const unsigned long long gei_pattern_ = 670526590282893600ull;
                const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
                const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
                const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
                const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
                return (3 * ugrid.gl_index(i, j, k) + offs);
            };
            // edge table unsigned long long e_table = 240177437832960;
            // input e and offset, which direction the normal has to point
            auto get_vertex_pos = [](const int e, const int offset)
            {
                const unsigned long long e_talbe = 240177437832960;
                return (e_talbe >> (4 * e + 2 * offset)) & 3;
            };
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mc_polygon(i0, f, c_);
            }
            else {
                cnt_ = mc_polygon(i_case, c_);
            }
            const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
            ushort e_{ 0 };
            // compute normals at vertices
            float3 n[8];
            ugrid.gradient(n, f, i_index, j_index, k_index);
            p_mc::GaussElimination::Matrix A;
            p_mc::GaussElimination::Vector X;
            p_mc::GaussElimination g;
            for (int t = 0; t < cnt_; t++)
            {
                float3 ei = make_float3(0, 0, 0);
                const int cnt_size = get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = get_c(t, i, c_);
                    getLocalCoordinates(l_edges_[e], e, f, i0, ei);
                    //set edge case to construct oriented quadrilateral
                    if (f[(l_edges_[e] & 0xF)] < i0) e_ |= (1 << e);
                }
                // normalize
                ei.x /= cnt_size;
                ei.y /= cnt_size;
                ei.z /= cnt_size;

                movePointToSurface(f, i0, ei);
                float3 ni = gradient(f, ei.x, ei.y, ei.z);
                normalize(ni);
                A(0, 0) = ni.x * ni.x;
                A(0, 1) = ni.x * ni.y;
                A(0, 2) = ni.x * ni.z;
                A(1, 0) = ni.y * ni.x;
                A(1, 1) = ni.y * ni.y;
                A(1, 2) = ni.y * ni.z;
                A(2, 0) = ni.z * ni.x;
                A(2, 1) = ni.z * ni.y;
                A(2, 2) = ni.z * ni.z;
                float d = ni.x * ei.x + ni.y * ei.y + ni.z * ei.z;
                A(0, 3) = ni.x * d;
                A(1, 3) = ni.y * d;
                A(2, 3) = ni.z * d;
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = get_c(t, i, c_);
                    float3 u{ 0,0,0 };
                    getLocalCoordinates(l_edges_[e], e, f, i0, u);
                    ni = gradient(f, u.x, u.y, u.z);
                    normalize(ni);
                    A(0, 0) += ni.x * ni.x;
                    A(0, 1) += ni.x * ni.y;
                    A(0, 2) += ni.x * ni.z;
                    A(1, 0) += ni.y * ni.x;
                    A(1, 1) += ni.y * ni.y;
                    A(1, 2) += ni.y * ni.z;
                    A(2, 0) += ni.z * ni.x;
                    A(2, 1) += ni.z * ni.y;
                    A(2, 2) += ni.z * ni.z;
                    d = ni.x * u.x + ni.y * u.y + ni.z * u.z;
                    A(0, 3) += ni.x * d;
                    A(1, 3) += ni.y * d;
                    A(2, 3) += ni.z * d;
                    u.x = 0.5 * (u.x + ei.x);
                    u.y = 0.5 * (u.y + ei.y);
                    u.z = 0.5 * (u.y + ei.z);
                    movePointToSurface(f, i0, u);
                    ni = gradient(f, u.x, u.y, u.z);
                    normalize(ni);
                    A(0, 0) += ni.x * ni.x;
                    A(0, 1) += ni.x * ni.y;
                    A(0, 2) += ni.x * ni.z;
                    A(1, 0) += ni.y * ni.x;
                    A(1, 1) += ni.y * ni.y;
                    A(1, 2) += ni.y * ni.z;
                    A(2, 0) += ni.z * ni.x;
                    A(2, 1) += ni.z * ni.y;
                    A(2, 2) += ni.z * ni.z;
                    d = ni.x * u.x + ni.y * u.y + ni.z * u.z;
                    A(0, 3) += ni.x * d;
                    A(1, 3) += ni.y * d;
                    A(2, 3) += ni.z * d;
                }
                X(0) = 0;
                X(1) = 0;
                X(2) = 0;
                if (g.gaussianElimination(A, X))
                {
                    // test if in domain
                    if ((X(0) >= 0.01 && X(0) <= 0.99) && (X(1) >= 0.01 && X(1) <= 0.99) && (X(2) >= 0.01 && X(2) <= 0.99))
                    {
                        ei.x = X(0);
                        ei.y = X(1);
                        ei.z = X(2);
                    }
                }

                // compute point in space (voxel grid)
                ni = trilinear(n, ei.x, ei.y, ei.z);
                ei.x = ugrid.x0 + (i_index + ei.x) * ugrid.dx;
                ei.y = ugrid.y0 + (j_index + ei.y) * ugrid.dy;
                ei.z = ugrid.z0 + (k_index + ei.z) * ugrid.dz;
                const int v_addr = v_.addVertex(ei, ni);

                // for all this edges save the vertex
                for (int i = 0; i < get_cnt_size(t, c_); i++)
                {
                    const uint e = get_c(t, i, c_);
                    // compute unique edges id
                    const int e_glId = e_glIndex(e, i_index, j_index, k_index, ugrid);
                    const int pos = get_vertex_pos(e, (e_ >> e) & 1);
                    //printf("add edge\n");
                    ht_.addVertex(e_glId, pos, v_addr);
                }
            }
        }
        __device__ bool movePointToSurface(const float f[8], const float i0, float3& u, float3& error)
        {
            const int max_iter = 20;
            float u1 = u.x, v1 = u.y, w1 = u.z;
            float val1 = trilinear(f, u1, v1, w1), val2;

            for (int iter = 0; iter < max_iter; iter++)
            {
                float3 grad = gradient(f, u1, v1, w1);
                if (val1 > i0)
                {
                    // invert gradient
                    grad.x = -grad.x;
                    grad.y = -grad.y;
                    grad.z = -grad.z;
                }
                // normalize
                normalize(grad);
                const float step = 0.05f;
                float u2 = u1 + step * grad.x;
                float v2 = v1 + step * grad.y;
                float w2 = w1 + step * grad.z;

                // check if we are within the cell
                bool c = clamp(u2, v2, w2);
                //const float val2 = trilinear(f, u2, v2, w2);
                val2 = trilinear(f, u2, v2, w2);
                if ((val1 <= i0 && i0 < val2) || (val2 < i0 && i0 <= val1))
                {
                    const float e = (i0 - val1) / (val2 - val1);
                    u.x = u1 + e * (u2 - u1);
                    u.y = v1 + e * (v2 - v1);
                    u.z = w1 + e * (w2 - w1);
                    error.x = 1;
                    break;
                }
                if (c) break;
                // update
                u1 = u2;
                v1 = v2;
                w1 = w2;
                val1 = val2;
            }
            return (error.x > 0);
        }
        /// Move point to surface based on gradient of scalar field
        /// @param [in] f values of scalar field at cell Vertices
        /// @param [in] i0 iso-values
        /// @param [in/out] starting position of cell representative, at end of
        ///                 iteration contains the point on the trilinear surface
        /// This is the original method proposed in the first version of the paper
        /// This method does not hit the surface under certain conditions.
        __device__ void movePointToSurfaceOld(const float f[8], const float i0, float3& u)
        {
            const int max_iter = 20;
            float u1 = u.x, v1 = u.y, w1 = u.z;
            float val1 = trilinear(f, u1, v1, w1), val2;
            for (int iter = 0; iter < max_iter; iter++)
            {
                float3 grad = gradient(f, u1, v1, w1);
                if (val1 > i0)
                {
                    // invert gradient
                    grad.x = -grad.x;
                    grad.y = -grad.y;
                    grad.z = -grad.z;
                }
                // normalize
                normalize(grad);
                const float step = 0.05f;
                float u2 = u1 + step * grad.x;
                float v2 = v1 + step * grad.y;
                float w2 = w1 + step * grad.z;

                // check if we are within the cell
                bool c = clamp(u2, v2, w2);
                //const float val2 = trilinear(f, u2, v2, w2);
                val2 = trilinear(f, u2, v2, w2);
                if ((val1 <= i0 && i0 < val2) || (val2 < i0 && i0 <= val1))
                {
                    const float e = (i0 - val1) / (val2 - val1);
                    u.x = u1 + e * (u2 - u1);
                    u.y = v1 + e * (v2 - v1);
                    u.z = w1 + e * (w2 - w1);
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
        /// <summary>
        /// Move point to surface using the gradient of the scalar field
        /// to set direction.
        /// </summary>
        /// <param name="f">scalar values at cell vertices</param>
        /// <param name="i0">iso-value</param>
        /// <param name="u">position of point, being moving to trilinear surface</param>
        /// <returns></returns>
        __device__ void movePointToSurface(const float f[8], const float i0, float3& u)
        {
            const int max_iter = 30;
            float u1 = u.x, v1 = u.y, w1 = u.z;
            float val1 = trilinear(f, u1, v1, w1), val2;
            for (int iter = 0; iter < max_iter; iter++)
            {
                float3 grad = gradient(f, u1, v1, w1);
                if (val1 > i0)
                {
                    // invert gradient
                    grad.x = -grad.x;
                    grad.y = -grad.y;
                    grad.z = -grad.z;
                }
                // normalize
                normalize(grad);
                const float step = 0.05f;
                float u2 = u1 + step * grad.x;
                float v2 = v1 + step * grad.y;
                float w2 = w1 + step * grad.z;

                // check if we are within the cell
                int c = clampDim(u2, v2, w2);
                val2 = trilinear(f, u2, v2, w2);
                if ((val1 <= i0 && i0 <= val2) || (val2 <= i0 && i0 <= val1))
                {
                    float e = val1;
                    if (val1 != val2)
                        e = (i0 - val1) / (val2 - val1);
                    u.x = u1 + e * (u2 - u1);
                    u.y = v1 + e * (v2 - v1);
                    u.z = w1 + e * (w2 - w1);
                    break;
                }

                if (c != -1) {
                    if (c == 0) {
                        grad.x = 0;
                    }
                    else if (c == 1) {
                        grad.y = 0;
                    }
                    else if (c == 2) {
                        grad.z = 0;
                    }
                    normalize(grad);
                    u2 = u1 + step * grad.x;
                    v2 = v1 + step * grad.y;
                    w2 = w1 + step * grad.z;
                    clampDim(u2, v2, w2);

                    val2 = trilinear(f, u2, v2, w2);
                    if ((val1 <= i0 && i0 <= val2) || (val2 <= i0 && i0 <= val1))
                    {
                        float e = val1;
                        if (val1 != val2)
                            e = (i0 - val1) / (val2 - val1);
                        u.x = u1 + e * (u2 - u1);
                        u.y = v1 + e * (v2 - v1);
                        u.z = w1 + e * (w2 - w1);
                        break;
                    }
                }
                // update
                u1 = u2;
                v1 = v2;
                w1 = w2;
                val1 = val2;
            }
        }
        /// trilinear interpolation of scalar
        __device__ float trilinear(const float f[8], const float u, const float v, const float w)
        {
            return (1 - w) * ((1 - v) * ((1 - u) * f[0] + u * f[1]) + v * ((1 - u) * f[2] + u * f[3])) + w * ((1 - v) * ((1 - u) * f[4] + u * f[5]) + v * ((1 - u) * f[6] + u * f[7]));
        }
        /// <summary>Trilinar Interpolation of a vector</summary>
        /// <param name="p">Points at the vertices</param>
        /// <param name="u">local coordinate</param>
        /// <param name="v">local coordinate</param>
        /// <param name="w">local coordinate</param>
        /// <returns></returns>
        __device__ float3 trilinear(const float3 p[8], const float u, const float v, const float w)
        {
            float3 po;
            po.x = (1 - w) * ((1 - v) * (p[0].x + u * (p[1].x - p[0].x)) + v * (p[2].x + u * (p[3].x - p[2].x))) + w * ((1 - v) * (p[4].x + u * (p[5].x - p[4].x)) + v * (p[6].x + u * (p[7].x - p[6].x)));
            po.y = (1 - w) * ((1 - v) * (p[0].y + u * (p[1].y - p[0].y)) + v * (p[2].y + u * (p[3].y - p[2].y))) + w * ((1 - v) * (p[4].y + u * (p[5].y - p[4].y)) + v * (p[6].y + u * (p[7].y - p[6].y)));
            po.z = (1 - w) * ((1 - v) * (p[0].z + u * (p[1].z - p[0].z)) + v * (p[2].z + u * (p[3].z - p[2].z))) + w * ((1 - v) * (p[4].z + u * (p[5].z - p[4].z)) + v * (p[6].z + u * (p[7].z - p[6].z)));
            return po;
        }
        /// <summary>
        /// compute the gradient of the trilinear interpolant within a cell
        /// </summary>
        /// <param name="f">scalar function values at the cell vertices</param>
        /// <param name="u">first local coordinate</param>
        /// <param name="v">second local coordinate</param>
        /// <param name="w">third local coordinate</param>
        /// <returns>gradient of a at point (u,v,w)</returns>
        __device__ float3 gradient(const float f[8], const float u, const float v, const float w)
        {
            float3 grad;
            grad.x = (1 - w) * ((1 - v) * (f[1] - f[0]) + v * (f[3] - f[2])) + w * ((1 - v) * (f[5] - f[4]) + v * (f[7] - f[6]));
            grad.y = (1 - w) * ((1 - u) * (f[2] - f[0]) + u * (f[3] - f[1])) + w * ((1 - u) * (f[6] - f[4]) + u * (f[7] - f[5]));
            grad.z = (1 - v) * ((1 - u) * (f[4] - f[0]) + u * (f[5] - f[1])) + v * ((1 - u) * (f[6] - f[2]) + u * (f[7] - f[3]));
            return grad;
        }
        /// <summary>
        /// Normalize a 3D vector
        /// </summary>
        /// <param name="v">inpute 3D vector</param>
        /// <returns></returns>
        __device__ void normalize(float3& v)
        {
            const float sz = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            v.x /= sz;
            v.y /= sz;
            v.z /= sz;
        }
        /// <summary>Clamp local coordinates to unit cube</summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <param name="w"></param>
        /// <returns>true is coordinated clamped</returns>
        __device__ bool clamp(float& u, float& v, float& w)
        {
            bool flag{ false };
            (u <= 0 ? (u = 0, flag = true) : u <= 1 ? u : (u = 1, flag = true));
            (v <= 0 ? (v = 0, flag = true) : v <= 1 ? v : (v = 1, flag = true));
            (w <= 0 ? (w = 0, flag = true) : w <= 1 ? w : (w = 1, flag = true));
            return flag;
        }
        /// <summary>
        /// Clamp local coordinates to unit cube and returns
        /// the dimension being clumped.
        /// </summary>
        /// <param name="u"></param>
        /// <param name="v"></param>
        /// <param name="w"></param>
        /// <returns>true is coordinated clamped</returns>
        __device__ int clampDim(float& u, float& v, float& w) {
            int flag{ -1 };
            (u <= 0 ? (u = 0, flag = 0) : u <= 1 ? u : (u = 1, flag = 0));
            (v <= 0 ? (v = 0, flag = 1) : v <= 1 ? v : (v = 1, flag = 1));
            (w <= 0 ? (w = 0, flag = 2) : w <= 1 ? w : (w = 1, flag = 2));
            return flag;
        }
        /// <summary>
        /// Compute and return local coordinates of input point u with respect to edge
        /// </summary>
        /// <param name="l_edge">encode vertices ids of edge</param>
        /// <param name="e">edge being processed</param>
        /// <param name="f">scalar values at cell vertices</param>
        /// <param name="i0">iso-value</param>
        /// <param name="u">local coordinates being accumulated</param>
        /// <returns></returns>
        __device__ float3 getLocalCoordinates(const unsigned char l_edge, const int e, const float f[8], const float i0, float3& u)
        {
            const int v0 = (l_edge & 0xF);
            const int v1 = (l_edge >> 4) & 0xF;
            const float l = (i0 - f[v0]) / (f[v1] - f[v0]);
            unsigned long long int e_{ 75059404284193ULL };
            unsigned long long c_{ 38552806359568ULL };
            float val = (e_ & 1ULL << (4 * e)) >> (4 * e);
            float con = (c_ & 1ULL << (4 * e)) >> (4 * e);
            u.x += l * val + con;
            val = (e_ & 1ULL << (4 * e + 1)) >> (4 * e + 1);
            con = (c_ & 1ULL << (4 * e + 1)) >> (4 * e + 1);
            u.y += l * val + con;
            val = (e_ & 1ULL << (4 * e + 2)) >> (4 * e + 2);
            con = (c_ & 1ULL << (4 * e + 2)) >> (4 * e + 2);
            u.z += l * val + con;
            return u;
        }

    };
} // namespace p_mc
