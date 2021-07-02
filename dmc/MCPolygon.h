#pragma once

// CUDA
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project
#include "helper_cuda.h"

namespace p_mc {
    /// <summary>
    /// Computes the MC polygon, i.e. the intersection of the iso-surface with the cell
    /// where the hyperbolic arcs are approximated with straight segment
    /// </summary>
    struct MCPolygon {
        //
        using ulong = unsigned long long;
        using uchar = unsigned char;
        using ushort = unsigned short;
        using uint = unsigned int;
        /// data
        const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };

        /// set size of mc polygon
        __device__ void set_cnt_size(const int cnt, const int size, ulong& c_) {
            // unset contour size
            c_ &= ~(0xFll << 4 * cnt);
            c_ |= (static_cast<ulong>(size) << 4 * cnt);
        }
        /// get size of mc polygon
        __device__ int get_cnt_size(const int cnt, unsigned long long& c_) {
            return static_cast<int>((c_ & (0xFll << 4 * cnt)) >> 4 * cnt);
        }
        /// set corresponding edge
        __device__ void set_c(const int cnt, const int pos, const int val, unsigned long long& c_) {
            //const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
            const uint c_sz = c_ & mask[cnt];
            const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
            c_ &= ~(((unsigned long long)0xF) << e);
            c_ |= (((unsigned long long)val) << e);
        }
        /// read edge from polygon
        __device__ int get_c(const int cnt, const int pos, unsigned long long c_) {
            //const uint mask[4] = { 0x0, 0xF, 0xFF, 0xFFF };
            const uint c_sz = (uint)(c_ & mask[cnt]);
            const uint e = 16 + 4 * ((c_sz & 0xF) + ((c_sz & 0xF0) >> 4) + ((c_sz & 0xF00) >> 8) + pos);
            return static_cast<int>((c_ >> e) & 0xF);
        }
        /// intersect a cell with iso-surface for unambiguous cases
        __device__ int mc_polygon(const int i_case, ulong& c_, const char* r_pattern)
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
        /// <summary>
        /// For each edge in the reference unit cell set to which edge the MC polygon is going. The variable
        /// segm_ will at the end encode for each of the 12 edges the edge from which the MC polygon is coming
        /// and the edge towards the MC polygon is going. This way, the MC polygon can be reconstructed
        /// starting at any edge.
        /// </summary>
        /// <param name="ei">edge id where MC polygon starts at the face being processed</param>
        /// <param name="eo">edge id toward the MC polygon is moving at the face being processed</param>
        /// <param name="segm_">encodes for each edge the index of the edge from where the polygon is coming
        /// and the index of the edge towards the MC polygon is moving. </param>
        __device__ void set_segm(const int ei, const int eo, unsigned char segm_[12]) {
            segm_[ei] &= 0xF0;
            segm_[ei] |= ((unsigned char)eo) & 0xF;
            segm_[eo] &= 0xF;
            segm_[eo] |= ((unsigned char)ei) << 4;
        }
        /// <summary>
        /// Get edge index of the edge toward or from which the MC polygon is going or coming.
        /// </summary>
        /// <param name="e">index of edge being processed</param>
        /// <param name="pos">pos=0 returns index of the edge toward the MC polygon moves, pos=1 the index of the
        /// edge from where the polygon is coming</param>
        /// <param name="segm_">encodes in and out for the given edge</param>
        /// <returns>edge index</returns>
        __device__ int get_segm(const int e, const int pos, unsigned char segm_[12]) {
            if (pos == 0)
                return (int)(segm_[e] & 0xF);
            else
                return (int)((segm_[e] >> 4) & 0xF);
        };
        /// <summary>
        /// Checks if MC polygon intersects the edge
        /// </summary>
        /// <param name="e">edge index</param>
        /// <param name="segm_">array of unsigned chars encoding edges state</param>
        /// <returns>true if edge is set</returns>
        __device__ bool is_segm_set(const int e, unsigned char segm_[12]) {
            return (segm_[e] != 0xFF);
        }
        /// <summary>
        /// Initialize the entries for a given edge
        /// </summary>
        /// <param name="e">edge index</param>
        /// <param name="segm_">array of unsigned chars encoding the state of the edges</param>
        /// <returns></returns>
        __device__ void unset_segm(const int e, unsigned char segm_[12]) {
            segm_[e] = 0xFF;
        }
        /// <summary>
        /// Given four scalar values computes the value of the function at the intersection of the asymptotes
        /// </summary>
        /// <param name="f0">function value at vertex 0</param>
        /// <param name="f1">function value at vertex 1</param>
        /// <param name="f2">function value at vertex 2</param>
        /// <param name="f3">function value at vertex 3</param>
        /// <returns>function value at intersection of asymptotes</returns>
        __device__ float asymptotic_decider(const float f0, const float f1, const float f2, const float f3) {
            return (f0 * f3 - f1 * f2) / (f0 + f3 - f1 - f2);
        }
        /// compute cell edge index for edge of input face
        __device__ int get_face_e(const int f, const int e, ushort e_face_[6]) { return ((e_face_[f] >> (4 * e)) & 0xF); };
        // compute cell vertex index for vertex of input face
        __device__ int get_face_v(const int f, const int e, ushort v_face_[6]) { return ((v_face_[f] >> (4 * e)) & 0xF); };
        /// intersect cell for ambiguous cases
        __device__ int mc_polygon(const float i0, const float F[8], unsigned long long& c_)
        {
            // compute oriented contours
            // 1. build segments
            // 2. connect segments
            // build up segments
            // set segments map
            unsigned char segm_[12] = { 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF };

            // In order to compute oriented segments, the hexahedron has to be flatten.
            // The insides of the faces of the hexahedron have to be all at the same
            // side of the flattened hexa. This requires changing the order of the
            // edges when reading from the faces
            // encode edges at face
            unsigned short e_face_[6]{ (ushort)291, (ushort)18277, (ushort)18696, (ushort)10859, (ushort)33719, (ushort)38305 };
            // encode vertices at face
            unsigned short v_face_[6]{ (ushort)12576, (ushort)25717, (ushort)5380, (ushort)29538, (ushort)8292, (ushort)30001 };

            // compute oriented segments using the isoline scheme at the faces
            for (int f = 0; f < 6; f++) {
                // classify face
                unsigned int f_case{ 0 };
                const int v0 = get_face_v(f, 0, v_face_);
                const int v1 = get_face_v(f, 1, v_face_);
                const int v2 = get_face_v(f, 2, v_face_);
                const int v3 = get_face_v(f, 3, v_face_);
                const int e0 = get_face_e(f, 0, e_face_);
                const int e1 = get_face_e(f, 1, e_face_);
                const int e2 = get_face_e(f, 2, e_face_);
                const int e3 = get_face_e(f, 3, e_face_);
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
            // 3) The next 12 entries are the indices of the edges constituting the contours
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

            // return number of contours
            return cnt_;
        }
    };

} // namespace p_mc
