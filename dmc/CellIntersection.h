#pragma once

// C++ libs
#include <memory>

// CUDA stuff
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// Project files
#include "helper_cuda.h"
#include "MCPolygon.h"
#include "VertexRepresentative.h"
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
        using ulong = unsigned long long;
        using UGrid = p_mc::UniformGrid;
        /// <summary>
        /// Lookup table with the MC polygons, which produces
        /// the correct result for non-ambiguous cases.
        /// </summary>
        char* r_pattern{ nullptr };
        std::shared_ptr<char> r_pattern_;
        /// <summary>
        /// Lookup table listing all ambiguous cases.
        /// </summary>
        uchar* t_ambig{ nullptr };
        std::shared_ptr<uchar> t_ambig_;
        /// <summary>
        /// Structure which contains methods to compute the MC polygon for a unit cell. These
        /// methods are only for ambiguous cases used.
        /// </summary>
        MCPolygon mcP;
        /// <summary>
        /// Structure with methods to compute the vertex representative of each branch of the
        /// iso-surface intersecting the cell.
        /// </summary>
        VertexRepresentative vRep;
        /// <summary>
        /// Atomic integer to count the number of cases, where the point projection onto the
        /// iso-surface fail.
        /// </summary>
        uint* ptProjection1{ nullptr };
        std::shared_ptr<uint> ptProjection1_;
        uint* ptProjection2{ nullptr };
        std::shared_ptr<uint> ptProjection2_;
        /// constructor
        __host__ CellIntersection() {}
        __host__ CellIntersection(const std::array<char, 4352>& p, const std::array<uchar, 256>& a)
        {
            // alloc and init r_pattern
            cudaMalloc(&r_pattern, 4352 * sizeof(char));
            cudaMemcpy(r_pattern, &p[0], 4352 * sizeof(char), cudaMemcpyHostToDevice);
            cudaCheckError();
            cudaMalloc(&t_ambig, 256 * sizeof(uchar));
            cudaMemcpy(t_ambig, &a[0], 256 * sizeof(uchar), cudaMemcpyHostToDevice);
            cudaCheckError();
            cudaMalloc(&ptProjection1, sizeof(uint));
            cudaMemset(ptProjection1, 0, sizeof(uint));
            cudaCheckError();
            cudaMalloc(&ptProjection2, sizeof(uint));
            cudaMemset(ptProjection2, 0, sizeof(uint));
            cudaCheckError();
            // set shared pointer to be sure memory is being freed
            r_pattern_ = std::shared_ptr<char>(r_pattern, cudaFree);
            t_ambig_ = std::shared_ptr<uchar>(t_ambig, cudaFree);
            ptProjection1_ = std::shared_ptr<uint>(ptProjection1, cudaFree);
            ptProjection2_ = std::shared_ptr<uint>(ptProjection2, cudaFree);
        }
        /// destructor
        __host__ ~CellIntersection()
        {
            r_pattern = nullptr;
            t_ambig = nullptr;
            ptProjection1 = nullptr;
            ptProjection2 = nullptr;
            r_pattern_.reset();
            t_ambig_.reset();
            ptProjection1_.reset();
            ptProjection2_.reset();
        }
        __host__ uint failedProjections1()
        {
            uint nrP{ 0 };
            cudaMemcpy(&nrP, ptProjection1, sizeof(uint), cudaMemcpyDeviceToHost);
            return nrP;
        }
        __host__ uint failedProjections2()
        {
            uint nrP{ 0 };
            cudaMemcpy(&nrP, ptProjection2, sizeof(uint), cudaMemcpyDeviceToHost);
            return nrP;
        }
        /// <summary>
        /// Compute edge global index
        /// </summary>
        /// <param name="e">edge index within the unit cell</param>
        /// <param name="i_idx">i-index of unit cell, i.e. index of origin vertex</param>
        /// <param name="j_idx">j-index of unit cell</param>
        /// <param name="k_idx">k-index of unit cell</param>
        /// <param name="ugrid">uniform grid, contains grid size information</param>
        /// <returns>edge global index</returns>
        __device__ int e_glIndex(const int e, const int i_idx, const int j_idx, const int k_idx, UGrid& ugrid)
        {
            const unsigned long long gei_pattern_ = 670526590282893600ull;
            const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
            const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
            const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
            const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
            return (3 * ugrid.gl_index(i, j, k) + offs);
        }
        /// <summary>
        /// Determines which is the position of the vertex index in the quadrilateral. It guaranties
        /// the generation of consistently oriented quads, where vertices with function values larger
        /// than the iso-value are outside the surface.
        /// </summary>
        /// <param name="e">edge index in the unit cell</param>
        /// <param name="offset">indicates which is the vertex case, i.e. larger than iso-value</param>
        /// <returns>vertex index position in quadrilateral</returns>
        __device__ int get_vertex_pos(const int e, const int offset)
        {
            const unsigned long long e_talbe = 240177437832960;
            return (e_talbe >> (4 * e + 2 * offset)) & 3;
        }
        /// count the nr. of MC polygons intersecting the cell
        __device__ int countMCPolygons(const float i0, const int i_case, float f[8])
        {
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mcP.mc_polygon(i0, f, c_);
            }
            else {
                cnt_ = mcP.mc_polygon(i_case, c_, r_pattern);
            }
            return static_cast<int>(cnt_);
        }
        /// <summary>
        /// Computes the surface vertices representing the intersection of the iso-surface
        /// with a cell. The MC polygon is computed using a modified MC lookup table for unambiguous cases
        /// or the asymptotic decider for ambiguous cases. The vertex representing the iso-surface for
        /// each branch intersecting the cell is computed by sampling the trilinear interpolant.
        /// </summary>
        /// <param name="i0">iso-value</param>
        /// <param name="i_case">MC case</param>
        /// <param name="i_index">i-index of cell in the volume</param>
        /// <param name="j_index">j-index of cell in the volume</param>
        /// <param name="k_index">k-index of cell in the volume</param>
        /// <param name="f">array containing the eight scalar values at the cell vertices</param>
        /// <param name="ugrid">uniform grid</param>
        /// <param name="ht_">quadrilateral hash table</param>
        /// <param name="v_">vertices</param>
        /// <returns></returns>
        __device__ void slice(const float i0, const int i_case, const int i_index, const int j_index, const int k_index,
            float f[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mcP.mc_polygon(i0, f, c_);
            }
            else {
                cnt_ = mcP.mc_polygon(i_case, c_, r_pattern);
            }
            vertexRepresentatives(i0, cnt_, c_, i_index, j_index, k_index, f, ugrid, ht_, v_);
        }
        /// <summary>
        /// Compute for each branch of the iso-surface intersecting the cell, the vertex representative.
        /// The trilinear interpolant is sampled until a good candidate is found. If the sampling process
        /// fails, the vertex representative is set to the mean value of the vertices of the MC polygon.
        /// </summary>
        /// <param name="i0"></param>
        /// <param name="cnt_"></param>
        /// <param name="c_"></param>
        /// <param name="i_index"></param>
        /// <param name="j_index"></param>
        /// <param name="k_index"></param>
        /// <param name="f"></param>
        /// <param name="ugrid"></param>
        /// <param name="ht_"></param>
        /// <param name="v_"></param>
        /// <returns></returns>
        __device__ void vertexRepresentatives(const float i0, uint cnt_, ulong& c_,
            const int i_index, const int j_index, const int k_index,
            float f[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            // compute normals at cell vertices
            float3 n[8];
            ugrid.gradient(n, f, i_index, j_index, k_index);
            // compute MC polygons
            const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
            // set flag at edge to compute oriented quads, it depends on gradient along edge
            // use the fact that the iso-surface in the trilinear sense intersects each edge only once.
            ushort e_{ 0 };
            // bounding box of MP polygon
            float2 bboxU[4]; // there can be at most four MC polygons
            float2 bboxV[4]; // there can be at most four MC polygons
            float2 bboxW[4]; // there can be at most four MC polygons
            // initial values for vertex representative
            //float3 p[4];
            for (int t = 0; t < cnt_; t++)
            {
                // init bbox for this MC polygon
                bboxU[t].x = 1; // min
                bboxU[t].y = 0; // max
                bboxV[t].x = 1; // min
                bboxV[t].y = 0; // max
                bboxW[t].x = 1; // min
                bboxW[t].y = 0; // max
                // init vertex representative
                //p[t].x = 0;
                //p[t].y = 0;
                //p[t].z = 0;
                const int cnt_size = mcP.get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    // id of edge being intersected
                    const int e = mcP.get_c(t, i, c_);
                    // coordinates of edge intersection
                    float3 ei = getLocalCoordinates(l_edges_[e], e, f, i0);
                    // compute bounding box of MC polygon
                    bboxU[t].x = bboxU[t].x > ei.x ? ei.x : bboxU[t].x;
                    bboxU[t].y = bboxU[t].y < ei.x ? ei.x : bboxU[t].y;
                    bboxV[t].x = bboxV[t].x > ei.y ? ei.y : bboxV[t].x;
                    bboxV[t].y = bboxV[t].y < ei.y ? ei.y : bboxV[t].y;
                    bboxW[t].x = bboxW[t].x > ei.z ? ei.z : bboxW[t].x;
                    bboxW[t].y = bboxW[t].y < ei.z ? ei.z : bboxW[t].y;

                    //set edge case to construct oriented quadrilateral
                    if (f[(l_edges_[e] & 0xF)] < i0) e_ |= (1 << e);
                }
            } // loop over MC polygons

            for (int t = 0; t < cnt_; t++)
            {
                float3 p = make_float3(0, 0, 0);
                const int cnt_size = mcP.get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    // id of edge being intersected
                    const int e = mcP.get_c(t, i, c_);
                    // coordinates of edge intersection
                    float3 ei = getLocalCoordinates(l_edges_[e], e, f, i0);
                    // compute initial vertex representative
                    p.x += ei.x;
                    p.y += ei.y;
                    p.z += ei.z;
                }
                p.x /= cnt_size;
                p.y /= cnt_size;
                p.z /= cnt_size;
                uint nrProjectionsFailed1{ 0 };
                uint nrProjectionsFailed2{ 0 };
                vRep.representative(i0, f, cnt_, c_, t, p, bboxU, bboxV, bboxW, nrProjectionsFailed1, nrProjectionsFailed2);
                if (nrProjectionsFailed1 > 0)
                {
                    // count how many times back projection failed
                    atomicAdd(ptProjection1, nrProjectionsFailed1);
                    if (nrProjectionsFailed2 > 0)
                    {
                        // count how many times back projection failed
                        atomicAdd(ptProjection2, nrProjectionsFailed2);
                    }
                }
                // compute normal
                float3 ni = trilinear(n, p.x, p.y, p.z);
                normalize(ni);
                // compute point in space (voxel grid)
                p.x = ugrid.x0 + (i_index + p.x) * ugrid.dx;
                p.y = ugrid.y0 + (j_index + p.y) * ugrid.dy;
                p.z = ugrid.z0 + (k_index + p.z) * ugrid.dz;
                // save unique vertex and normal
                const int v_addr = v_.addVertex(p, ni);
                // for all this edges save the vertex
                for (int i = 0; i < mcP.get_cnt_size(t, c_); i++)
                {
                    // compute element
                    const uint e = mcP.get_c(t, i, c_);
                    int color{ -1 };
                    if (e == 0) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2);
                    if (e == 3) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2) + 1;
                    if (e == 8) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2) + 2;
                    // compute unique edges id
                    const int e_glId = e_glIndex(e, i_index, j_index, k_index, ugrid);
                    const int pos = get_vertex_pos(e, (e_ >> e) & 1);
                    if (color == -1) ht_.addVertex(e_glId, pos, v_addr);
                    else  ht_.addVertex(e_glId, pos, v_addr, color);
                }
            }
        }
        /// <summary>
        /// check if MC case is ambiguous.
        /// </summary>
        /// <param name="i_case">MC case</param>
        /// <returns>true if MC case is ambiguous</returns>
        __device__ bool isAmbiguous(const int i_case)
        {
            if (t_ambig[i_case] == MC_AMBIGUOUS) return true;
            else return false;
        }
        /// <summary>
        /// Compute surface representative by moving the mean value of the vertices of the MC
        /// polygon towards the iso-surface along the direction of the gratien.
        /// </summary>
        /// <param name="i0">iso-value</param>
        /// <param name="i_case">MC case</param>
        /// <param name="i_index">i-index of cell in volume</param>
        /// <param name="j_index">j-index of cell in volume</param>
        /// <param name="k_index">k-index of cell in volume</param>
        /// <param name="f">scalar values at cell vertices</param>
        /// <param name="ugrid">uniform grid</param>
        /// <param name="ht_">quadrilateral hash table</param>
        /// <param name="v_">vertices</param>
        /// <returns></returns>
        __device__ void sliceP(const float i0, const int i_case, const int i_index, const int j_index, const int k_index,
            float f[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mcP.mc_polygon(i0, f, c_);
            }
            else {
                cnt_ = mcP.mc_polygon(i_case, c_, r_pattern);
            }
            const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
            ushort e_{ 0 };
            // compute normals at vertices
            float3 n[8];
            ugrid.gradient(n, f, i_index, j_index, k_index);
            for (int t = 0; t < cnt_; t++)
            {
                float3 ei = make_float3(0, 0, 0);
                const int cnt_size = mcP.get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = mcP.get_c(t, i, c_);
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
                    const uint e = mcP.get_c(t, i, c_);
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
        /// <summary>
        /// Intersect cell with iso-surface. The point repositioning is based on error quadrics.
        /// The cell intersection if computed as follows
        /// <list>
        /// <item>
        /// <term>MC polygon</term>
        /// <description>Compute first the MC polygon. For non-ambiguous cases use the lookup table.
        /// Ambiguous case compute the MC polygon based on the asymptotic decider.</description>
        /// </item>
        /// <item>
        /// <term>Representative</term>
        /// <description>
        /// <para>Compute intersection of iso-surface with cell edges</param>
        /// <para>Use the points computed in the previous step and the interpolated
        ///      normals to computed the surface representative using error quadrics.</para>
        /// </desctiption>
        /// </item>
        /// <term>Element indices</term>
        /// <description>Computed elements indices based on function gradient to obtain
        ///              a consistently oriented quad only mesh.</description>
        /// </list>
        /// </summary>
        /// <param name="i0">iso-value</param>
        /// <param name="i_case">MC case</param>
        /// <param name="i_index">i-index of origin vertex of the cell</param>
        /// <param name="j_index">j-index of origin vertex of the cell</param>
        /// <param name="k_index">k-index of origin vertex of the cell</param>
        /// <param name="f">function values at cell vertices</param>
        /// <param name="ugrid">structure with methods for the underlying uniform grid</param>
        /// <param name="ht_">hash table to compute element global indices</param>
        /// <param name="v_">array containing the vertices</param>
        /// <returns>void</returns>
        __device__ void sliceQ(const float i0, const int i_case, const int i_index, const int j_index, const int k_index,
            float f[8], UGrid& ugrid, QuadrilateralHashTable& ht_, Vertices& v_)
        {
            /// compute mc polygons
            unsigned long long c_ = 0xFFFFFFFFFFFF0000;
            uint cnt_{ 0 };
            if (t_ambig[i_case] == MC_AMBIGUOUS)
            {
                cnt_ = mcP.mc_polygon(i0, f, c_);
            }
            else {
                cnt_ = mcP.mc_polygon(i_case, c_, r_pattern);
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
                const int cnt_size = mcP.get_cnt_size(t, c_);
                for (int i = 0; i < cnt_size; i++)
                {
                    const uint e = mcP.get_c(t, i, c_);
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
                    const uint e = mcP.get_c(t, i, c_);
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
                for (int i = 0; i < mcP.get_cnt_size(t, c_); i++)
                {
                    const uint e = mcP.get_c(t, i, c_);
                    // compute unique edges id
                    const int e_glId = e_glIndex(e, i_index, j_index, k_index, ugrid);
                    const int pos = get_vertex_pos(e, (e_ >> e) & 1);
                    //printf("add edge\n");
                    ht_.addVertex(e_glId, pos, v_addr);
                }
            }
        }
        /// <summary>
        /// Move representative towards iso-surface following the direction of the gradient. If point
        /// steps outside cell, clamp movement in this direction by setting gradient component to zero.
        /// </summary>
        /// <param name="f">array function values at cell vertices</param>
        /// <param name="i0">iso-value</param>
        /// <param name="u">local coordinates of point representative</param>
        /// <returns></returns>
        __device__ void movePointToSurface(const float f[8], const float i0, float3& u)
        {
            const int max_iter = 30;
            float u1 = u.x, v1 = u.y, w1 = u.z;
            float val1 = trilinear(f, u1, v1, w1), val2;
            //float c0 = val1;
            //float u0 = u1, v0 = v1, w0 = w1;
            //int flag = 0;
            //float3 grad;
            bool foundSurface = false;
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
                //const float step = l * (i0 - val1) / (grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
                float u2 = u1 + step * grad.x;
                float v2 = v1 + step * grad.y;
                float w2 = w1 + step * grad.z;

                // check if we are within the cell
                int c = clampDim(u2, v2, w2);
                //const float val2 = trilinear(f, u2, v2, w2);
                val2 = trilinear(f, u2, v2, w2);
                if ((val1 <= i0 && i0 <= val2) || (val2 <= i0 && i0 <= val1))
                {
                    float e = val1;
                    if (val1 != val2)
                        e = (i0 - val1) / (val2 - val1);
                    u.x = u1 + e * (u2 - u1);
                    u.y = v1 + e * (v2 - v1);
                    u.z = w1 + e * (w2 - w1);
                    foundSurface = true;
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

                    //const float val2 = trilinear(f, u2, v2, w2);
                    val2 = trilinear(f, u2, v2, w2);
                    if ((val1 <= i0 && i0 <= val2) || (val2 <= i0 && i0 <= val1))
                    {
                        float e = val1;
                        if (val1 != val2)
                            e = (i0 - val1) / (val2 - val1);
                        u.x = u1 + e * (u2 - u1);
                        u.y = v1 + e * (v2 - v1);
                        u.z = w1 + e * (w2 - w1);
                        foundSurface = true;
                        break;
                    }
                }


                // update
                u1 = u2;
                v1 = v2;
                w1 = w2;
                val1 = val2;
            }
            if (!foundSurface) {
                // count number of case, where back projection to iso-surface failed
                atomicAdd(ptProjection1, 1);
                //printf("Surface not found\n");
                //printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", u0, v0, w0, u1, v1, w1, f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7], c0, val1);
            }
        }
        /// <summary>
        /// Trilinear interpolation of a scalar value at local coordinates
        /// </summary>
        /// <param name="f">array with function values at vertices</param>
        /// <param name="u">local coordinate</param>
        /// <param name="v">local coordinate</param>
        /// <param name="w">local coordinate</param>
        /// <returns>interpolated function value</returns>
        __device__ float trilinear(const float f[8], const float u, const float v, const float w)
        {
            return (1 - w) * ((1 - v) * ((1 - u) * f[0] + u * f[1]) + v * ((1 - u) * f[2] + u * f[3])) + w * ((1 - v) * ((1 - u) * f[4] + u * f[5]) + v * ((1 - u) * f[6] + u * f[7]));
        }
        /// <summary>
        /// Trilinar Interpolation of a vector using positions or gradients at cell vertices. Input the local
        /// coordinates for the trilinear interpolation.
        /// </summary>
        /// <param name="p">Point or gradient at cell vertices</param>
        /// <param name="u">local coordinate</param>
        /// <param name="v">local coordinate</param>
        /// <param name="w">local coordinate</param>
        /// <returns>interpolated point</returns>
        __device__ float3 trilinear(const float3 p[8], const float u, const float v, const float w)
        {
            float3 po;
            po.x = (1 - w) * ((1 - v) * (p[0].x + u * (p[1].x - p[0].x)) + v * (p[2].x + u * (p[3].x - p[2].x))) + w * ((1 - v) * (p[4].x + u * (p[5].x - p[4].x)) + v * (p[6].x + u * (p[7].x - p[6].x)));
            po.y = (1 - w) * ((1 - v) * (p[0].y + u * (p[1].y - p[0].y)) + v * (p[2].y + u * (p[3].y - p[2].y))) + w * ((1 - v) * (p[4].y + u * (p[5].y - p[4].y)) + v * (p[6].y + u * (p[7].y - p[6].y)));
            po.z = (1 - w) * ((1 - v) * (p[0].z + u * (p[1].z - p[0].z)) + v * (p[2].z + u * (p[3].z - p[2].z))) + w * ((1 - v) * (p[4].z + u * (p[5].z - p[4].z)) + v * (p[6].z + u * (p[7].z - p[6].z)));
            return po;
        }
        /// <summary>
        /// Compute the gradient of the trilinear interpolant within a cell at the given local coordinates.
        /// </summary>
        /// <param name="f">scalar function values at the cell vertices</param>
        /// <param name="u">first local coordinate</param>
        /// <param name="v">second local coordinate</param>
        /// <param name="w">third local coordinate</param>
        /// <returns>gradient of trilinear interpolant at point (u,v,w)</returns>
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
        /// <param name="v">input vector</param>
        /// <returns>void</returns>
        __device__ void normalize(float3& v)
        {
            const float sz = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            v.x /= sz;
            v.y /= sz;
            v.z /= sz;
        }
        /// <summary>
        /// Clamp local coordinates to unit cube.
        /// </summary>
        /// <param name="u">local coordinate</param>
        /// <param name="v">local coordinate</param>
        /// <param name="w">local coordinate</param>
        /// <returns>true if coordinated clamped</returns>
        __device__ bool clamp(float& u, float& v, float& w)
        {
            bool flag{ false };
            (u <= 0 ? (u = 0, flag = true) : u <= 1 ? u : (u = 1, flag = true));
            (v <= 0 ? (v = 0, flag = true) : v <= 1 ? v : (v = 1, flag = true));
            (w <= 0 ? (w = 0, flag = true) : w <= 1 ? w : (w = 1, flag = true));
            return flag;
        }
        /// <summary>
        /// Clamp local coordinates to unit cube and return the dimension being clamped.
        /// </summary>
        /// <param name="u">u local coord.</param>
        /// <param name="v">v local coord.</param>
        /// <param name="w">w local coord.</param>
        /// <returns>dimension clampoed</returns>
        __device__ int clampDim(float& u, float& v, float& w) {
            int flag{ -1 };
            (u <= 0 ? (u = 0, flag = 0) : u <= 1 ? u : (u = 1, flag = 0));
            (v <= 0 ? (v = 0, flag = 1) : v <= 1 ? v : (v = 1, flag = 1));
            (w <= 0 ? (w = 0, flag = 2) : w <= 1 ? w : (w = 1, flag = 2));
            return flag;
        }
        /// <summary>
        /// Computes the intersection of the mc polygon with the edge and add the position
        /// to the vertex representative for this polygon. This way, the vertex representative
        /// can be computed as the mean value of the intersections of the mc polygon with the cell.
        /// The vertex representative is given in local coordinates.
        /// </summary>
        /// <param name="l_edge">lookup table to get indices of end vertices of edge</param>
        /// <param name="e">edge id</param>
        /// <param name="f">function values at cell vertices</param>
        /// <param name="i0">iso-value</param>
        /// <param name="u">vertex representative</param>
        /// <returns>void</returns>
        __device__ void getLocalCoordinates(const unsigned char l_edge, const int e, const float f[8], const float i0, float3& u)
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
            //return u;
        }
        /// <summary>
        /// Computes the local coordinates of the intersection of the MC polygon with the edge.
        /// </summary>
        /// <param name="l_edge">lookup table to get indices of end vertices of edge</param>
        /// <param name="e">edge id</param>
        /// <param name="f">function values at cell vertices</param>
        /// <param name="i0">iso-value</param>
        /// <returns>local coordinates of iso-surface intersection with the edge</returns>
        __device__ float3 getLocalCoordinates(const unsigned char l_edge, const int e, const float f[8], const float i0)
        {
            const int v0 = (l_edge & 0xF);
            const int v1 = (l_edge >> 4) & 0xF;
            const float l = (i0 - f[v0]) / (f[v1] - f[v0]);
            unsigned long long int e_{ 75059404284193ULL };
            unsigned long long c_{ 38552806359568ULL };
            float val = (e_ & 1ULL << (4 * e)) >> (4 * e);
            float con = (c_ & 1ULL << (4 * e)) >> (4 * e);
            float3 u;
            u.x = l * val + con;
            val = (e_ & 1ULL << (4 * e + 1)) >> (4 * e + 1);
            con = (c_ & 1ULL << (4 * e + 1)) >> (4 * e + 1);
            u.y = l * val + con;
            val = (e_ & 1ULL << (4 * e + 2)) >> (4 * e + 2);
            con = (c_ & 1ULL << (4 * e + 2)) >> (4 * e + 2);
            u.z = l * val + con;
            return u;
        }

    };
} // namespace p_mc
