#include "DualMarchingCubes.h"
#include "LookupTables.h"


void cpp_mc::DualMarchingCubes::dual_mc(const double i0, UGrid& ugrid,
	std::vector<Vertex>& v, std::vector<Normal>& n,
	std::vector<int>& tris, std::vector<int>& quads)
{
	// collect information about the uniform grid
	const int idim = ugrid.x_size();
	const int jdim = ugrid.y_size();
	const int kdim = ugrid.z_size();
	std::map<int, std::array<int, 5>> m_;
    std::cout << " ... compute iso-surface" << std::endl;
    for (int k = 0; k < (kdim - 1); k++)
	{
		for (int j = 0; j < (jdim - 1); j++)
		{
			for (int i = 0; i < (idim - 1); i++)
			{
				double u[8];
				u[0] = ugrid.scalar(i, j, k);
				u[1] = ugrid.scalar(i + 1, j, k);
				u[2] = ugrid.scalar(i, j + 1, k);
				u[3] = ugrid.scalar(i + 1, j + 1, k);
				u[4] = ugrid.scalar(i, j, k + 1);
				u[5] = ugrid.scalar(i + 1, j, k + 1);
				u[6] = ugrid.scalar(i, j + 1, k + 1);
				u[7] = ugrid.scalar(i + 1, j + 1, k + 1);

				//
				uint i_case{ 0 };
				i_case = i_case + ((uint)(u[0] >= i0));
				i_case = i_case + ((uint)(u[1] >= i0)) * 2;
				i_case = i_case + ((uint)(u[2] >= i0)) * 4;
				i_case = i_case + ((uint)(u[3] >= i0)) * 8;
				i_case = i_case + ((uint)(u[4] >= i0)) * 16;
				i_case = i_case + ((uint)(u[5] >= i0)) * 32;
				i_case = i_case + ((uint)(u[6] >= i0)) * 64;
				i_case = i_case + ((uint)(u[7] >= i0)) * 128;

				if (i_case == 0 || i_case == 255)
					continue;
				else
					slice(i0, i_case, i, j, k, u, ugrid, v, n, m_);
			}
		}
	}

	// collect quadrilaterals
    std::vector<int> colors;
    std::cout << " ... collect quadrilaterals" << std::endl;
	collectQuadrilaterals(quads, colors, m_);
    // simplify 3X3Y
    std::cout << " ... mesh simplification" << std::endl;
    simplify3X3Y(v, n, quads, colors);
	// compute triangles
    std::cout << " ... compute triangles" << std::endl;
	collectTriangles(tris, quads, v);
}

void cpp_mc::DualMarchingCubes::slice(const double i0, const uint i_case, const int i_index, const int j_index, const int k_index,
	double f[8], UGrid& ugrid, std::vector<Vertex>& v, std::vector<Normal>& n, std::map<int, std::array<int, 5>>& m_)
{
	auto e_glIndex = [](const int e, const int i_idx, const int j_idx, const int k_idx, UGrid& ugrid)
	{
		const unsigned long long gei_pattern_ = 670526590282893600ull;
		const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
		const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
		const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
		const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
		return (3 * ugrid.global_index(i, j, k) + offs);
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
	// compute vertex representative for each mc polygon
	for (uint t = 0; t < cnt_; t++)
	{
		Vertex ei;
		const int cnt_size = get_cnt_size(t, c_);
		for (int i = 0; i < cnt_size; i++)
		{
			const uint e = get_c(t, i, c_);
			getLocalCoordinates(l_edges_[e], e, f, i0, ei);
			//set edge case to construct oriented quadrilateral
			if (f[(l_edges_[e] & 0xF)] < i0) e_ |= (1 << e);
		}
		// normalize
		ei /= cnt_size;

		movePointToSurface(f, i0, ei);
		// compute normal at mesh vertex
		Normal ni;
		ugrid.normal(ni, i_index, j_index, k_index, ei[0], ei[1], ei[2]);
		// compute euclidean coordinates of new vertex
		Point pi;
		ugrid.position(pi, i_index, j_index, k_index, ei[0], ei[1], ei[2]);
		// add vertex and normal to list
		const int v_addr = static_cast<int>(v.size());
		v.push_back(pi);
		n.push_back(ni);

		// for all this edges save the vertex
		for (int i = 0; i < cnt_size; i++)
		{
			const uint e = get_c(t, i, c_);
			// compute unique edges id
			const int e_glId = e_glIndex(e, i_index, j_index, k_index, ugrid);
			const int pos = get_vertex_pos(e, (e_ >> e) & 1);
            // compute color
            int color{ -1 };
            if (e == 0) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2);
            if (e == 3) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2) + 1;
            if (e == 8) color = 3 * ((i_index & 1) | (j_index & 1) << 1 | (k_index & 1) << 2) + 2;
			//add vertex to has table
            if (color == -1)
                addVertex(e_glId, pos, v_addr, m_);
            else
                addVertex(e_glId, pos, v_addr, color, m_);
		}
	}
}


unsigned int cpp_mc::DualMarchingCubes::mc_polygon(const double i0, const double F[8], ulong& c_)
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
    auto asymptotic_decider = [](const double f0, const double f1, const double f2, const double f3) {
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
        const double f0 = F[v0];
        const double f1 = F[v1];
        const double f2 = F[v2];
        const double f3 = F[v3];
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
            const double val = asymptotic_decider(f0, f1, f2, f3);
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
            const double val = asymptotic_decider(f0, f1, f2, f3);
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

unsigned int cpp_mc::DualMarchingCubes::mc_polygon(const int i_case, ulong& c_)
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

// Mesh simplification
void cpp_mc::DualMarchingCubes::halfedges(const int nr_v, std::vector<int>& quads, std::vector<Halfedge>& he, std::vector<int>& he_v, std::vector<int>& he_f)
{
    const int nr_q = static_cast<int>(quads.size()) / 4;
    he_v.resize(nr_v);
    he_f.resize(nr_q);
    he.resize(4 * nr_q);
    std::map<ulong, std::array<int, 5>> m_;
    auto setKey = [](const int v0, const int v1)
    {
        if (v0 < v1)
            return (static_cast<ulong>(v0) << 32) | (v1 & 0xffffffffL);
        else
            return (static_cast<ulong>(v1) << 32) | (v0 & 0xffffffffL);
    };
    auto addHalfedge = [](const ulong key, const int he, std::map<ulong, std::array<int,5>>& m)
    {
        auto e = m.find(key);
        if (e != m.end())
        {
            // the element alread exist
            const int pos = e->second[4];
            e->second[pos] = he;
            e->second[4] = pos + 1;
        }
        else
        {
            // the element has to be created
            std::array<int, 5> h{ INVALID_INDEX,INVALID_INDEX,INVALID_INDEX,INVALID_INDEX, 1 };
            h[0] = he;
            m[key] = h;
        }
    };
    for (int i = 0; i < nr_q; i++)
    {
        // quad vertices
        const int v0 = quads[4 * i];
        const int v1 = quads[4 * i + 1];
        const int v2 = quads[4 * i + 2];
        const int v3 = quads[4 * i + 3];
        // he[0] = origin vertex
        // he[1] = face
        // he[2] = next
        // he[3] = tween
        // he[4] = 0 - manifold, 1 non-manifold
        std::array<int, 5> e0_{ INVALID_INDEX, i, INVALID_INDEX, INVALID_INDEX, MANIFOLD };
        std::array<int, 5> e1_{ INVALID_INDEX, i, INVALID_INDEX, INVALID_INDEX, MANIFOLD };
        std::array<int, 5> e2_{ INVALID_INDEX, i, INVALID_INDEX, INVALID_INDEX, MANIFOLD };
        std::array<int, 5> e3_{ INVALID_INDEX, i, INVALID_INDEX, INVALID_INDEX, MANIFOLD };
        // origin vertex
        e0_[0] = v0;
        e1_[0] = v1;
        e2_[0] = v2;
        e3_[0] = v3;
        // next
        e0_[2] = 4 * i + 1;
        e1_[2] = 4 * i + 2;
        e2_[2] = 4 * i + 3;
        e3_[2] = 4 * i;
        // set
        he[4 * i] = e0_;
        he[4 * i + 1] = e1_;
        he[4 * i + 2] = e2_;
        he[4 * i + 3] = e3_;
        // collect faces
        he_f[i] = 4 * i; // this face is pointing to the "first" halfedge
        // collect vertices
        he_v[v0] = 4 * i;
        he_v[v1] = 4 * i + 1;
        he_v[v2] = 4 * i + 2;
        he_v[v3] = 4 * i + 3;
        // add halfedges to hash table
        addHalfedge(setKey(v0, v1), 4 * i, m_);
        addHalfedge(setKey(v1, v2), 4 * i + 1, m_);
        addHalfedge(setKey(v2, v3), 4 * i + 2, m_);
        addHalfedge(setKey(v3, v0), 4 * i + 3, m_);
    }
    // connect halfedge twins
    int nr_nonmanifold{ 0 };
    for (auto e : m_)
    {
        const int c = e.second[4];
        switch (c)
        {
        case 1:
        {
            // this is a boundary edge
            const int e0 = e.second[0];
            he[e0][3] = INVALID_INDEX;
            break;
        }
        case 2:
        {
            const int e0 = e.second[0];
            const int e1 = e.second[1];
            he[e0][3] = e1;
            he[e1][3] = e0;
            break;
        }
        case 4:
        {
            nr_nonmanifold++;
            const int e0 = e.second[0];
            const int e1 = e.second[1];
            const int e2 = e.second[2];
            const int e3 = e.second[3];
            // collect vertices
            const int v0 = he[e0][0];
            const int v1 = he[e1][0];
            const int v2 = he[e2][0];
            const int v3 = he[e3][0];
            if (v0 != v1)
            {
                he[e0][3] = e1;
                he[e1][3] = e0;
                he[e2][3] = e3;
                he[e3][3] = e2;
            }
            else
            {
                he[e0][3] = e2;
                he[e2][3] = e0;
                he[e1][3] = e3;
                he[e3][3] = e1;
            }
            he[e0][4] = NON_MANIFOLD;
            he[e1][4] = NON_MANIFOLD;
            he[e2][4] = NON_MANIFOLD;
            he[e3][4] = NON_MANIFOLD;
            break;
        }
        default:
            std::cout << "ERROR: wrong nr. of faces sharing an edge: " << c << std::endl;
            break;
        }
    }
    std::cout << " ... nr. of non-manifold edges: " << nr_nonmanifold << std::endl;
}

std::array<int, 4> cpp_mc::DualMarchingCubes::collectNeighbors(const int quad, std::vector<Halfedge>& he, std::vector<int>& he_f)
{
    // he[0] = origin vertex
    // he[1] = face
    // he[2] = next
    // he[3] = tween
    // he[4] = 0 - manifold, 1 non-manifold
    // collect the four halfedges
    const int e0 = he_f[quad];
    const int e1 = he[e0][2];
    const int e2 = he[e1][2];
    const int e3 = he[e2][2];
    // collect neighbors
    const int t0 = he[e0][3];
    const int t1 = he[e1][3];
    const int t2 = he[e2][3];
    const int t3 = he[e3][3];
    // collect faces
    int f0{ -1 };
    int f1{ -1 };
    int f2{ -1 };
    int f3{ -1 };
    // set neighbors
    if (t0 != INVALID_INDEX) f0 = he[t0][1];
    if (t1 != INVALID_INDEX) f1 = he[t1][1];
    if (t2 != INVALID_INDEX) f2 = he[t2][1];
    if (t3 != INVALID_INDEX) f3 = he[t3][1];

    return { f0,f1,f2,f3 };

}

/// <summary>
/// Faces hinterit from uniform grid a coloring with 24 colors, the gird coloring is extended to an grid edge coloring
/// with 24 colors. The number of colors has to be reduced to 4, if possible. Proceed in two steps
/// 1. reduce the number of colors to 5: 0, 1, 2, 3, 4
/// 2. re-color as many faces as possible that have color 4 (the fifth color)
/// </summary>
/// <param name="he">halfedges</param>
/// <param name="he_f">halfedge faces</param>
/// <param name="colors">face colors</param>
void cpp_mc::DualMarchingCubes::colorFaces(std::vector<Halfedge>& he, std::vector<int>& he_f, std::vector<int>& colors)
{
    const int nr_q{ static_cast<int>(he_f.size()) };
    // classify faces with colors larger than 5
    std::array<std::list<int>, 19> fC;
    std::list<int> fColor4;
    for (int f = 0; f < nr_q; f++)
    {
        if (colors[f] > 4) // color larger than fifth color
        {
            fC[colors[f] - 5].push_back(f);
        }
        else if (colors[f] == 4)
        {
            fColor4.push_back(f);
        }
    }
    // for each color simplify
    for (auto c : fC)
    {
        for (auto f : c)
        {
            // collect neighbors
            std::array<int, 4> n = collectNeighbors(f, he, he_f);
            // collect collors
            int c0{ INVALID_COLOR };
            int c1{ INVALID_COLOR };
            int c2{ INVALID_COLOR };
            int c3{ INVALID_COLOR };
            if (n[0] != INVALID_INDEX) c0 = colors[n[0]];
            if (n[1] != INVALID_INDEX) c1 = colors[n[1]];
            if (n[2] != INVALID_INDEX) c2 = colors[n[2]];
            if (n[3] != INVALID_INDEX) c3 = colors[n[3]];
            // this colors must be all larger than 4
            if (c0 != 0 && c1 != 0 && c2 != 0 && c3 != 0)
            {
                colors[f] = 0;
                continue;
            }
            if (c0 != 1 && c1 != 1 && c2 != 1 && c3 != 1)
            {
                colors[f] = 1;
                continue;
            }
            if (c0 != 2 && c1 != 2 && c2 != 2 && c3 != 2)
            {
                colors[f] = 2;
                continue;
            }
            if (c0 != 3 && c1 != 3 && c2 != 3 && c3 != 3)
            {
                colors[f] = 3;
                continue;
            }
            if (c0 != 4 && c1 != 4 && c2 != 4 && c2 != 4)
            {
                colors[f] = 4;
                continue;
            }
        }
    }
    // opimize colors, reduce the number of faces with color 4 as much as possible
    //for (int f = 0; f < nr_q; f++)
    for (auto f : fColor4)
    {
        // collect neighbors
        std::array<int, 4> n = collectNeighbors(f, he, he_f);
        // collect collors
        int c0{ INVALID_COLOR };
        int c1{ INVALID_COLOR };
        int c2{ INVALID_COLOR };
        int c3{ INVALID_COLOR };
        if (n[0] != INVALID_INDEX) c0 = colors[n[0]];
        if (n[1] != INVALID_INDEX) c1 = colors[n[1]];
        if (n[2] != INVALID_INDEX) c2 = colors[n[2]];
        if (n[3] != INVALID_INDEX) c3 = colors[n[3]];
        // this colors must be all larger than 4
        if (c0 != 0 && c1 != 0 && c2 != 0 && c3 != 0)
        {
            colors[f] = 0;
            continue;
        }
        if (c0 != 1 && c1 != 1 && c2 != 1 && c3 != 1)
        {
            colors[f] = 1;
            continue;
        }
        if (c0 != 2 && c1 != 2 && c2 != 2 && c3 != 2)
        {
            colors[f] = 2;
            continue;
        }
        if (c0 != 3 && c1 != 3 && c2 != 3 && c3 != 3)
        {
            colors[f] = 3;
            continue;
        }
    }
    // check which colors remain
    std::array<int, 6> c_{ 0,0,0,0,0,0 };
    for (auto c : colors)
    {
        switch (c)
        {
        case 0:
            c_[0]++;
            break;
        case 1:
            c_[1]++;
            break;
        case 2:
            c_[2]++;
            break;
        case 3:
            c_[3]++;
            break;
        case 4:
            c_[4]++;
            break;
        default:
            c_[5]++;
            break;
        }
    }
    for (auto c : c_)
    {
        std::cout << " ... colors: " << c << std::endl;
    }
}

void cpp_mc::DualMarchingCubes::vertexValence(const int nr_v, std::vector<Halfedge>& he_, std::vector<int>& vV)
{
    // he[0] = origin vertex
    // he[1] = face
    // he[2] = next
    // he[3] = tween
    // he[4] = 0 - manifold, 1 non-manifold
    const int nr_he{ static_cast<int>(he_.size()) };
    vV.resize(nr_v);
    std::fill(vV.begin(), vV.end(), 0);
    for (int e = 0; e < nr_he; e++)
    {
        vV[he_[e][0]]++;
        // check if boundary edge
        if (he_[e][3] == INVALID_INDEX)
        {
            const int ne{ he_[e][2] };
            vV[he_[ne][0]]++;
        }
    }
}

bool cpp_mc::DualMarchingCubes::isNonManifold(const int f, std::vector<Halfedge>& he, std::vector<int>& he_f)
{
    const int e0 = he_f[f];
    const int e1 = he[e0][2];
    const int e2 = he[e1][2];
    const int e3 = he[e2][2];
    if (he[e0][4] == NON_MANIFOLD || he[e1][4] == NON_MANIFOLD || he[e2][4] == NON_MANIFOLD || he[e3][4] == NON_MANIFOLD)
        return true;
    else
        return false;
}
void cpp_mc::DualMarchingCubes::mark3X3Y(std::vector<int>& quads, std::vector<int>& vV, std::vector<bool>& p3X3Y)
{
    const int nr_q{ static_cast<int>(quads.size()) / 4 };
    p3X3Y.resize(nr_q);
    std::fill(p3X3Y.begin(), p3X3Y.end(), false);
    for (int f = 0; f < nr_q; f++)
    {
        const int v0{ quads[4 * f] };
        const int v1{ quads[4 * f + 1] };
        const int v2{ quads[4 * f + 2] };
        const int v3{ quads[4 * f + 3] };
        const int valence0 = vV[v0];
        const int valence1 = vV[v1];
        const int valence2 = vV[v2];
        const int valence3 = vV[v3];

        bool flag1 = (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5);
        bool flag2 = (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3);
        if (flag1 || flag2)
        {
            p3X3Y[f] = true;
        }
    }
}

// Color based simplification of elements with vertex valence pattern 3X3Y
void cpp_mc::DualMarchingCubes::mergeVertices3X3Y(std::vector<Vertex>& v, std::vector<Normal>& normals, std::vector<bool> p3X3Y,
    std::vector<int>& vV, std::vector<int>& colors, std::vector<Halfedge>& he, std::vector<int>& he_f,
    std::vector<std::pair<bool, int>>& vm_, std::vector<bool>& em_)
{
    const int nr_q{ static_cast<int>(he_f.size()) };
    const int nr_v{ static_cast<int>(v.size()) };
    vm_.resize(nr_v);
    em_.resize(nr_q, false);
    std::fill(vm_.begin(), vm_.end(), std::make_pair(false, INVALID_INDEX));
    std::fill(em_.begin(), em_.end(), false);
    for (int f = 0; f < nr_q; f++)
    {
        const int c = colors[f];
        if (!p3X3Y[f]) continue;

        // if one edge in non-manifold, do not remove
        // he[0] = origin vertex
        // he[1] = face
        // he[2] = next
        // he[3] = tween
        // he[4] = 0 - manifold, 1 non-manifold
        const int e0 = he_f[f];
        const int e1 = he[e0][2];
        const int e2 = he[e1][2];
        const int e3 = he[e2][2];
        //if (he[e0][4] == NON_MANIFOLD || he[e1][4] == NON_MANIFOLD || he[e2][4] == NON_MANIFOLD || he[e3][4] == NON_MANIFOLD) continue;
        // collect neghibors and check pattern and colors, be careful not to be at the boundary
        std::array<int, 4> n = collectNeighbors(f, he, he_f);
        bool nonManifoldNeighbor{ false };
        for (auto e : n)
        {
            if (isNonManifold(e,he,he_f))
                nonManifoldNeighbor = true;
        }
        if (nonManifoldNeighbor) continue;
        int n_color[4]{ INVALID_COLOR,INVALID_COLOR, INVALID_COLOR, INVALID_COLOR };
        bool n_pattern[4]{ false,false,false,false };
        for (int i = 0; i < 4; i++)
        {
            if (n[i] != INVALID_INDEX)
            {
                n_color[i] = colors[n[i]];
                n_pattern[i] = p3X3Y[n[i]];
            }
        }
        // check if element can be removed
        bool flag{ true };
        for (int i = 0; i < 4; i++)
        {
            if (n_pattern[i] && (n_color[i] <= c)) flag = false;
        }
        if (!flag) continue;
        // the element can be removed
        const int v0 = he[e0][0];
        const int v1 = he[e1][0];
        const int v2 = he[e2][0];
        const int v3 = he[e3][0];

        const int valence0 = vV[v0];
        const int valence1 = vV[v1];
        const int valence2 = vV[v2];
        const int valence3 = vV[v3];

        if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5)
        {
            // compute new position and normals of v0
            v[v0] = v[v1] + 0.5 * (v[v3] - v[v1]);
            normals[v0] = normals[v1] + 0.5 * (normals[v3] - normals[v1]);
            normals[v0].normalize();

            // mark v2 to be removed
            vm_[v2].first = true;

            // set twins of v2 to be v0, to be able to remove element later
            vm_[v2].second = v0;

            // element has to be removed
            em_[f] = true;
        }
        else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3)
        {
            // compute new position and normal of v1
            v[v1] = v[v0] + 0.5 * (v[v2] - v[v0]);
            normals[v1] = normals[v0] + 0.5 * (normals[v2] - normals[v0]);
            normals[v1].normalize();

            // mark v3 to be removed
            vm_[v3].first = true;
            // set twins, remove v3, use addres of v1
            vm_[v3].second = v1;

            // element has to be removed
            em_[f] = true;
        }
    }
}
// Simplification of elements with vertex valence pattern 3X3Y, only consider neighbors
void cpp_mc::DualMarchingCubes::mergeVertices3X3Y(std::vector<Vertex>& v, std::vector<Normal>& normals, std::vector<bool> p3X3Y,
    std::vector<int>& vV, std::vector<Halfedge>& he, std::vector<int>& he_f,
    std::vector<std::pair<bool, int>>& vm_, std::vector<bool>& em_)
{
    const int nr_q{ static_cast<int>(he_f.size()) };
    const int nr_v{ static_cast<int>(v.size()) };
    vm_.resize(nr_v);
    em_.resize(nr_q, false);
    std::fill(vm_.begin(), vm_.end(), std::make_pair(false, INVALID_INDEX));
    std::fill(em_.begin(), em_.end(), false);
    for (int f = 0; f < nr_q; f++)
    {
        if (!p3X3Y[f]) continue;

        // if one edge in non-manifold, do not remove
        // he[0] = origin vertex
        // he[1] = face
        // he[2] = next
        // he[3] = tween
        // he[4] = 0 - manifold, 1 non-manifold
        const int e0 = he_f[f];
        const int e1 = he[e0][2];
        const int e2 = he[e1][2];
        const int e3 = he[e2][2];
        // collect neghibors and check pattern and colors, be careful not to be at the boundary
        std::array<int, 4> n = collectNeighbors(f, he, he_f);
        if (n[0] != INVALID_INDEX && isNonManifold(n[0], he, he_f)) continue;
        if (n[1] != INVALID_INDEX && isNonManifold(n[1], he, he_f)) continue;
        if (n[2] != INVALID_INDEX && isNonManifold(n[2], he, he_f)) continue;
        if (n[3] != INVALID_INDEX && isNonManifold(n[3], he, he_f)) continue;

        // check if element has neighbor with same vertex valence pattern
        if (n[0] != INVALID_INDEX && p3X3Y[n[0]]) continue;
        if (n[1] != INVALID_INDEX && p3X3Y[n[1]]) continue;
        if (n[2] != INVALID_INDEX && p3X3Y[n[2]]) continue;
        if (n[3] != INVALID_INDEX && p3X3Y[n[3]]) continue;

        // the element can be removed
        const int v0 = he[e0][0];
        const int v1 = he[e1][0];
        const int v2 = he[e2][0];
        const int v3 = he[e3][0];

        const int valence0 = vV[v0];
        const int valence1 = vV[v1];
        const int valence2 = vV[v2];
        const int valence3 = vV[v3];

        if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5)
        {
            // compute new position and normals of v0
            v[v0] = v[v1] + 0.5 * (v[v3] - v[v1]);
            normals[v0] = normals[v1] + 0.5 * (normals[v3] - normals[v1]);
            normals[v0].normalize();

            // mark v2 to be removed
            vm_[v2].first = true;

            // set twins of v2 to be v0, to be able to remove element later
            vm_[v2].second = v0;

            // element has to be removed
            em_[f] = true;
        }
        else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3)
        {
            // compute new position and normal of v1
            v[v1] = v[v0] + 0.5 * (v[v2] - v[v0]);
            normals[v1] = normals[v0] + 0.5 * (normals[v2] - normals[v0]);
            normals[v1].normalize();

            // mark v3 to be removed
            vm_[v3].first = true;
            // set twins, remove v3, use addres of v1
            vm_[v3].second = v1;

            // element has to be removed
            em_[f] = true;
        }
    }
}

void cpp_mc::DualMarchingCubes::removeVertices3X3Y(std::vector<Vertex>& v, std::vector<Normal>& n, std::vector<std::pair<bool, int>>& vm_,
    std::vector<Vertex>& nv, std::vector<Normal>& nn)
{
    const int nr_v{ static_cast<int>(v.size()) };
    nv.reserve(nr_v);
    nn.reserve(nr_v);
    for (int i = 0; i < nr_v; i++)
    {
        if (!vm_[i].first)
        {
            const int addr{ static_cast<int>(nv.size()) };
            nv.push_back(v[i]);
            nn.push_back(n[i]);
            vm_[i].second = addr;
        }
    }
}

void cpp_mc::DualMarchingCubes::removeQuadrilaterals3X3Y(std::vector<int>& q, std::vector<bool>& em,
    std::vector<std::pair<bool, int>>& vm, std::vector<int>& nq)
{
    const int nr_q{ static_cast<int>(q.size()) / 4 };
    nq.reserve(nr_q);
    for (int f = 0; f < nr_q; f++)
    {
        if (!em[f])
        {
            const int v0 = q[4 * f];
            const int v1 = q[4 * f + 1];
            const int v2 = q[4 * f + 2];
            const int v3 = q[4 * f + 3];
            if (vm[v0].first) nq.push_back(vm[vm[v0].second].second);
            else nq.push_back(vm[v0].second);
            if (vm[v1].first) nq.push_back(vm[vm[v1].second].second);
            else nq.push_back(vm[v1].second);
            if (vm[v2].first) nq.push_back(vm[vm[v2].second].second);
            else nq.push_back(vm[v2].second);
            if (vm[v3].first) nq.push_back(vm[vm[v3].second].second);
            else nq.push_back(vm[v3].second);
        }
    }
}

void cpp_mc::DualMarchingCubes::simplify3X3Y(std::vector<Vertex>& v, std::vector<Normal>& n, std::vector<int>& quads, std::vector<int>& colors)
{
    int nr_v = static_cast<int>(v.size());
    std::vector<std::array<int, 5>> he_; // mark manifold and non-manifold halfedges
    std::vector<int> v_;
    std::vector<int> f_;
    std::vector<int> vV_;
    std::vector<bool> p3X3Y;
    std::vector<std::pair<bool, int>> vm_;
    std::vector<bool> em_;
    std::vector<Vertex> nv;
    std::vector<Normal> nn;
    std::vector<int> nq;
    halfedges(nr_v, quads, he_, v_, f_);
    colorFaces(he_, f_, colors);
    vertexValence(nr_v, he_, vV_);
    mark3X3Y(quads, vV_, p3X3Y);
    mergeVertices3X3Y(v, n, p3X3Y, vV_, colors, he_, f_, vm_, em_);
    removeVertices3X3Y(v, n, vm_, nv, nn);
    removeQuadrilaterals3X3Y(quads, em_, vm_, nq);
    // copy elements back to input arrays
    v.resize(nv.size());
    std::copy(nv.begin(), nv.end(), v.begin());
    n.resize(nn.size());
    std::copy(nn.begin(), nn.end(), n.begin());
    quads.resize(nq.size());
    std::copy(nq.begin(), nq.end(), quads.begin());
    nv.clear();
    nn.clear();
    nq.clear();

    // Remove elements with vertex valence pattern 3X3Y if they do not have a
    // neighbor with the same vertex valence pattern
    nr_v = static_cast<int>(v.size());
    halfedges(nr_v, quads, he_, v_, f_);
    vertexValence(nr_v, he_, vV_);
    mark3X3Y(quads, vV_, p3X3Y);
    mergeVertices3X3Y(v, n, p3X3Y, vV_, he_, f_, vm_, em_);
    removeVertices3X3Y(v, n, vm_, nv, nn);
    removeQuadrilaterals3X3Y(quads, em_, vm_, nq);
    // copy elements back to input arrays
    v.resize(nv.size());
    std::copy(nv.begin(), nv.end(), v.begin());
    n.resize(nn.size());
    std::copy(nn.begin(), nn.end(), n.begin());
    quads.resize(nq.size());
    std::copy(nq.begin(), nq.end(), quads.begin());
}
