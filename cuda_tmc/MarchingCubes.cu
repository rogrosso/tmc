#include "MarchingCubes.h"

// Lauch bounds
#define THREADS_PER_BLOCK 256
#if __CUDA_ARCH__ >= 200
	#define MY_KERNEL_MAX_THREADS (2 * THREADS_PER_BLOCK)
	#define MY_KERNEL_MIN_BLOCKS 3
#else
	#define MY_KERNEL_MAX_THREADS THREADS_PER_BLOCK
	#define MY_KERNEL_MIN_BLOCKS 2
#endif

// defines
#define BIT_1 0x1
#define BIT_2 0x2
#define BIT_3 0x4
#define BIT_4 0x8
#define BIT_5 0x10
#define BIT_6 0x20
#define BIT_7 0x40
#define BIT_8 0x80
#define BIT_16 0x8000

// Empty Bucket
#define EMPTY_BUCKET_32 -1
#define EMPTY_BUCKET_64 0ull


// Shared memory experiments
#define AMB_BLOCKSIZE 64
#define MC_BLOCKSIZE 512

// type aliases
// Introduce convenient aliases here
using uint = unsigned int;
using uchar = unsigned char;
using ushort = unsigned short;
using ullong = unsigned long long;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// error handling
#define cudaCheckError() {                                          \
  cudaError_t e=cudaGetLastError();                                  \
  if(e!=cudaSuccess) {                                               \
  printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
  exit(0); \
      }                                                                  \
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Convenience function
template<typename T>
int size(T t) {
    int s{ 0 };
    cudaMemcpy(&s, t.t_size, sizeof(int), cudaMemcpyDeviceToHost);
    return s;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The hash table to obtain a unique vertex index
struct VertexHashTable {
    int* key{ nullptr };
    int* addr{ nullptr };
    int t_size;
};

void initVertexHashTable(VertexHashTable& ht, const int size) {
    ht.t_size = size;
    cudaMalloc(&ht.key, size * sizeof(int));
    cudaMalloc(&ht.addr, size * sizeof(int));
    //cudaMemset(key, EMPTY_BUCKET_32, size * sizeof(int));
}
void freeVertexHashTable(VertexHashTable& h) {
    if (h.key != nullptr) {
        cudaFree(h.key);
    }
    if (h.addr != nullptr) {
        cudaFree(h.addr);
    }
    h.t_size = 0;
    h.key = nullptr;
    h.addr = nullptr;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vertex array
struct Vertices {
    float4* vertices{ nullptr };
    float4* normals{ nullptr };
    int* t_size{ nullptr }; // atomic counter to compute nr. of vertices
};

void initVertices(Vertices& v, const int size) {
    cudaMalloc(&v.vertices, size * sizeof(float4));
    cudaMalloc(&v.normals, size * sizeof(float4));
    cudaMalloc(&v.t_size, sizeof(int));
    cudaMemset(v.t_size, 0, sizeof(int));
}
void freeVertices(Vertices& v) {
    if (v.vertices != nullptr) {
        cudaFree(v.vertices);
    }
    if (v.normals != nullptr) {
        cudaFree(v.normals);
    }
    if (v.t_size != nullptr) {
        cudaFree(v.t_size);
    }
    v.t_size = nullptr;
    v.vertices = nullptr;
    v.normals = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a triangle consist of three indices
struct Triangles {
    int a_size; // size of buffer
    int4* triangles{ nullptr };
    int* t_size{ nullptr }; // atomic counter to compute nr. of triangles
};


void initTriangles(Triangles& t, const int size) {
    t.a_size = size;
    cudaMalloc(&t.triangles, size * sizeof(int4));
    cudaMalloc(&t.t_size, sizeof(int));
    cudaMemset(t.t_size, 0, sizeof(int));
}

void freeTriangles(Triangles& t) {
    if (t.triangles != nullptr) {
        cudaFree(t.triangles);
    }
    if (t.t_size != nullptr) {
        cudaFree(t.t_size);
    }
    t.a_size = 0;
    t.triangles = nullptr;
    t.t_size = nullptr;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Track cell cases and ids
struct CellsIds {
    int* cells_{ nullptr };
    int* t_size{ nullptr }; // atomic counter to get address of ambiguous cell in a_cells array
};

void initCells(CellsIds& c, const int size) {
    cudaMalloc(&c.cells_, size * sizeof(int));
    cudaMalloc(&c.t_size, sizeof(int));
    cudaMemset(c.t_size, 0, sizeof(int));
}
void freeCells(CellsIds& c) {
    if (c.cells_ != nullptr) {
        cudaFree(c.cells_);
    }
    if (c.t_size != nullptr) {
        cudaFree(c.t_size);
    }
    c.t_size = nullptr;
    c.cells_ = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Track ambiguous cases
struct AmbiguousCells {
    int* a_cells{ nullptr };
    int* t_size{ nullptr }; // atomic counter to get address of ambiguous cell in a_cells array
};

void initACells(AmbiguousCells& ac, const int size) {
    cudaMalloc(&ac.a_cells, size * sizeof(int));
    cudaMalloc(&ac.t_size, sizeof(int));
    cudaMemset(ac.t_size, 0, sizeof(int));
}
void freeACells(AmbiguousCells& ac) {
    if (ac.a_cells != nullptr) {
        cudaFree(ac.a_cells);
    }
    if (ac.t_size != nullptr) {
        cudaFree(ac.t_size);
    }
    ac.t_size = nullptr;
    ac.a_cells = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hash table to map halfedge twins
struct HalfedgeHashTable {
    int t_size{ 0 };
    unsigned long long* key{ nullptr };
    int2* he_ids{ nullptr };
};

void initHalfedgeHashTable(HalfedgeHashTable& t, const int size) {
    t.t_size = size;
    cudaMalloc(&t.key, size * sizeof(unsigned long long));
    //cudaCheckError();
    cudaMemset(t.key, 0, size * sizeof(unsigned long long));
    //cudaCheckError();
    cudaMalloc(&t.he_ids, size * sizeof(int2));
}


__device__ bool addHalfedgeToHashTable (HalfedgeHashTable t, const int addr, const int v0, const int v1) {
    unsigned long long x = (unsigned long long)v0;
    unsigned long long y = (unsigned long long)v1;
    unsigned long long key = (x < y) ? y : x;
    key = key + (x + y) * (x + y + 1) / 2ull;
    
    {
        key = (~key) + (key << 21); // key = (key << 21) - key - 1;
        key = key ^ (key >> 24);
        key = (key + (key << 3)) + (key << 8); // key * 265
        key = key ^ (key >> 14);
        key = (key + (key << 2)) + (key << 4); // key * 21
        key = key ^ (key >> 28);
        key = key + (key << 31);
    }
    // open hashing
    int h = int(key % (unsigned long long)t.t_size);
    int e = 1;
    for (int loop = 0; loop < 128; loop++) {
        unsigned long long old = atomicCAS(&t.key[h], EMPTY_BUCKET_64, key);
        if (old == EMPTY_BUCKET_64 || old == key) {
            if (v0 < v1) {
                t.he_ids[h].x = addr;
            } 
            else {
                t.he_ids[h].y = addr;
            }
            return true;
        }
        else {
            // step with linear probing
            h = (h + e*e) % t.t_size;
            e = e + 1;
        }
    }
    //printf("ERROR: can't add halfedge\n");
    return false;
}

void freeHalfedgeHashTable(HalfedgeHashTable& t) {
    if (t.key != nullptr) {
        cudaFree(t.key);
    }
    if (t.he_ids != nullptr) {
        cudaFree(t.he_ids);
    }
    // 
    t.t_size = 0;
    t.key = nullptr;
    t.he_ids = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Halfedge data structure
struct Halfedges {
    int* t_size{ nullptr };
    int buffSize{ 0 };
    // halfedge int4: 
    //    he.x = origin vertex
    //    he.y = face
    //    he.z = next
    //   he.w = tween
    int4* he_e{ nullptr };
    int*  he_v{ nullptr }; // helfedge id
    int*  he_f{ nullptr };
};


void initHalfedges(Halfedges& h, const int nr_he, const int nr_v, const int nr_t) {
    h.buffSize = nr_he;
    cudaMalloc(&h.he_e, nr_he * sizeof(int4));
    //cudaCheckError();
    cudaMalloc(&h.he_v, nr_v * sizeof(int));
    //cudaCheckError();
    cudaMalloc(&h.he_f, nr_t * sizeof(int));
    cudaMalloc(&h.t_size, sizeof(int));
    cudaMemset(&h.t_size, 0, sizeof(int));
}
//__device__ void add(HalfedgeHashTable het_, const int v0, const int v1, const int v2) {
//    const int a_ = atomicAdd(t_size, 3);
//    const int f_ = a_ / 3;
//    // he 0
//    he_e[a_].x = v0;
//    he_e[a_].y = f_; // there are three halfedges for each face (triangle)
//    he_e[a_].z = a_ + 1; // next
//    he_e[a_].w = -1; // default is boundary edge
//
//                     // he 1
//    he_e[a_ + 1].x = v1;
//    he_e[a_ + 1].y = f_; // there are three halfedges for each face (triangle)
//    he_e[a_ + 1].z = a_ + 2;
//    he_e[a_ + 1].w = -1; // default is boundary edge
//
//                         // he 2
//    he_e[a_ + 2].x = v2;
//    he_e[a_ + 2].y = f_; // there are three halfedges for each face (triangle)
//    he_e[a_ + 2].z = a_;
//    he_e[a_ + 2].w = -1; // default is boundary edge
//
//                         // add halfedges ids to hash table
//    het_.add(a_, v0, v1);
//    het_.add(a_ + 1, v1, v2);
//    het_.add(a_ + 2, v2, v0);
//}

void freeHalfedges(Halfedges& h) {
    if (h.t_size != nullptr) {
        cudaFree(h.t_size);
    }
    if (h.he_e != nullptr) {
        cudaFree(h.he_e);
    }
    if (h.he_v != nullptr) {
        cudaFree(h.he_v);
    }
    if (h.he_f != nullptr) {
        cudaFree(h.he_f);
    }
    //
    h.buffSize = 0;
    h.t_size = nullptr;
    h.he_e = nullptr;
    h.he_v = nullptr;
    h.he_f = nullptr;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MC lookup tables
struct MC_lookup {
    ushort* e_{ nullptr };
    unsigned long long* t_{ nullptr };
};

void initMC_lookup(MC_lookup& l, const std::array<unsigned short, 256>& ep_, const std::array<int, 4096>& tp_, const std::array<unsigned char, 256>& ta_) {
    ushort* le_ = new ushort[256];
    for (int i = 0; i < 256; i++) {
        le_[i] = ep_[i];
        le_[i] |= (ta_[i] == 105) ? BIT_16 : 0x0;
    }
    cudaMalloc(&l.e_, 256 * sizeof(ushort));
    cudaMemcpy(l.e_, &le_[0], 256 * sizeof(ushort), cudaMemcpyHostToDevice);
    cudaCheckError();
    // create MC loolup table
    unsigned long long* l_ = new unsigned long long[256];
    unsigned long long flg = 0xFull;
    for (int i = 0; i < 256; i++) {
        int i_case = i * 16;
        unsigned long long f = 0ull;

        for (int t = 0; t < 16; t++) {
            int mcval = tp_[i_case + t];
            unsigned long long lmcval = (unsigned long long)mcval;
            if (mcval == -1) {
                f |= (flg << (t * 4));
            }
            else {
                f |= (lmcval << (t * 4));
            }
        }
        l_[i] = f;
    }
    cudaMalloc(&l.t_, 256 * sizeof(unsigned long long));
    cudaMemcpy(l.t_, &l_[0], 256 * sizeof(unsigned long long), cudaMemcpyHostToDevice);
    cudaCheckError();
    delete[] le_;
    delete[] l_;
}


void freeMC_lookup(MC_lookup& l)
{
    if (l.e_ != nullptr) {
        cudaFree(l.e_);
    }
    if (l.t_ != nullptr) {
        cudaFree(l.t_);
    }
    l.e_ = nullptr;
    l.t_ = nullptr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Uniform grid
struct UniformGrid {
    int idim{ 0 };
    int jdim{ 0 };
    int kdim{ 0 };
    float x0{ 0 };
    float y0{ 0 };
    float z0{ 0 };
    float dx{ 0 };
    float dy{ 0 };
    float dz{ 0 };

    __device__ int gl_index(const int i, const int j, const int k) {
        return (k * jdim * idim + j * idim + i);
    }
    __device__ int i_index(const int gl_index) {
        return (gl_index % idim);
    }
    __device__ int j_index(const int gl_index) {
        return ((gl_index / idim) % jdim);
    }
    __device__ int k_index(const int gl_index) {
        return (gl_index / (idim * jdim));
    }
    __host__ void size(const int x_size, const int y_size, const int z_size) {
        idim = x_size;
        jdim = y_size;
        kdim = z_size;
    }
    __host__ void origin(const float x, const float y, const float z) {
        x0 = x;
        y0 = y;
        z0 = z;
    }
    __host__ void spacing(const float x, const float y, const float z) {
        dx = x;
        dy = y;
        dz = z;
    }

};

using UGrid = UniformGrid;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// timer based on CUDA timen routines
struct CTimer {
    float e_milliseconds;
    cudaEvent_t c_start;
    cudaEvent_t c_stop;
    CTimer() {
        cudaEventCreate(&c_start);
        cudaEventCreate(&c_stop);
    }

    void __host__ start() {
        cudaEventRecord(c_start);
    }
    void __host__ stop() {
        cudaEventRecord(c_stop);
        cudaEventSynchronize(c_stop);
        cudaEventElapsedTime(&e_milliseconds, c_start, c_stop);
    }
    void __host__ print() {
        std::cout << std::setprecision(7) << " ... time in ms: " << e_milliseconds << std::endl;
    }
    void __host__ print(std::string& m) {
        std::cout << std::setprecision(7) << " ... " << m << " time in ms: " << e_milliseconds << std::endl;
    }

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
using namespace p_mc;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute number of vertices computed for this cell
// compute only the intersection of the iso-surface with
// the cell edges
template<typename T>
__device__ uint numberOfSetBits(T n) {
    // C or C++: use uint32_t
    uint b = (uint)n;
    b = b - ((b >> 1) & 0x55555555);
    b = (b & 0x33333333) + ((b >> 2) & 0x33333333);
    return (((b + (b >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reduce 
// warp reduce based on __shfl_down
template<typename T>
__device__ int warpReduceSum(T val) {
	for (int offset = warpSize / 2; offset > 0; offset /= 2) {
		val += __shfl_down(val, offset);
	}
	return val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// warp reduce kernel
template<typename T>
__global__ void warp_reduce_kernel(T *in, T* out, int N) {
	T sum = T(0);
	for (int i = blockIdx.x*blockDim.x + threadIdx.x; i<N; i += blockDim.x*gridDim.x) {
		sum += in[i];
	}
	sum = warpReduceSum(sum);
	if (threadIdx.x%warpSize == 0)
		atomicAdd(out, sum);
}

// host function for warp reduce
template<typename T>
int warpReduce(T *i_data, const int N) {
	int threads = 256;
	int blocks = std::min((N + threads - 1) / threads, 2048);

	T* d_sum{ nullptr };
	cudaMalloc(&d_sum, sizeof(int));
	cudaMemsetAsync(d_sum, 0, sizeof(int));
	warp_reduce_kernel<typename T> << <blocks, threads >> >(i_data, d_sum, N);
    cudaCheckError();
    // return sum
	T h_sum{ 0 };
	cudaMemcpy(&h_sum, d_sum, sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(d_sum);

	return h_sum;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Index computations
//__host__ __device__ int global_index(const int i, const int j, const int k, const int idim, const int jdim, const int kdim)
//{
//	return (k * jdim * idim + j * idim + i);
//}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// trilinear interpolation
__device__ void trilinear(float4& po, const float3 p[8], const float u, const float v, const float w)
{
    po.x = (1 - w) * ((1 - v) * (p[0].x + u * (p[1].x - p[0].x)) + v * (p[2].x + u * (p[3].x - p[2].x))) + w * ((1 - v) * (p[4].x + u * (p[5].x - p[4].x)) + v * (p[6].x + u * (p[7].x - p[6].x)));
    po.y = (1 - w) * ((1 - v) * (p[0].y + u * (p[1].y - p[0].y)) + v * (p[2].y + u * (p[3].y - p[2].y))) + w * ((1 - v) * (p[4].y + u * (p[5].y - p[4].y)) + v * (p[6].y + u * (p[7].y - p[6].y)));
    po.z = (1 - w) * ((1 - v) * (p[0].z + u * (p[1].z - p[0].z)) + v * (p[2].z + u * (p[3].z - p[2].z))) + w * ((1 - v) * (p[4].z + u * (p[5].z - p[4].z)) + v * (p[6].z + u * (p[7].z - p[6].z)));
    //po.y = (1 - w) * ((1 - v) * (p[0].y * (1 - u) + p[1].y * u) + v * (p[2].y * (1 - u) + p[3].y * u)) + w * ((1 - v) * (p[4].y * (1 - u) + p[5].y * u) + v * (p[6].y * (1 - u) + p[7].y * u));
    //po.z = (1 - w) * ((1 - v) * (p[0].z * (1 - u) + p[1].z * u) + v * (p[2].z * (1 - u) + p[3].z * u)) + w * ((1 - v) * (p[4].z * (1 - u) + p[5].z * u) + v * (p[6].z * (1 - u) + p[7].z * u));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hash tables
// device hash function
__device__ uint hash_function( uint key )
{
    key = (key ^ 61) ^ (key >> 16);
    key = key + (key << 3);
    key = key ^ (key >> 4);
    key = key * 0x27d4eb2d;
    key = key ^ (key >> 15);
    return key;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Insert a key assigned to a vertex index in the hash table
// return the position in the array, where the key was inserted
// if the key was already in the array, e.g. other kernel was processing
// the same edge in the uniform grid, return false. This way, the kernel will
// not generate a new vertex. If the vertex was not created jet, the address in
// the array is returned so that the calling kernel can save this position
// of the vertex in the hash table
// v_gindex key is a unique global index assigned to the vertex
// Hash table:
//   struct HashTable {
//       int* key;
//       int* addr;
//       int t_size;
//   };
// v_addr contains the position in the key array, where the key = v_gindex was stored
__device__ bool insert_vertex_key(const int v_gindex, VertexHashTable ht_, int& v_addr)
{
    //const int start_address = int(hash_function((uint)v_gindex) % ht_.t_size);
    // open hashing
    //int h = start_address;
    int h = int(hash_function((uint)v_gindex) % ht_.t_size);
    int e = 1;
    for (int loop = 0; loop < 128; loop++) {
        int old = atomicCAS(&ht_.key[h], EMPTY_BUCKET_32, v_gindex);
        if (old == EMPTY_BUCKET_32) {
            v_addr = h;
            return true;
        }
        else if (v_gindex == old) {
            // vertex key already in table
            return false;
        }
        else {
            // step with linear probing
            h = (h + e*e) % ht_.t_size;
            e = e + 1;
            //if (h == start_address) {
              //  printf("ERROR: can't find free bucket for %d\n", v_gindex);
               // return false;
            //}
        }
    }

    return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a probably faster strategy to reduce number of glabal memory access
// This function sets vertex and normal using an atomic counter.
// Keep the address where vertex was stored, therefore the hast table knows the address of vertices and normals in vertex and normal arrays
// It resturn the address in the hash table where the address of vertex and normal are stored
__device__ int insert_vertex(const int v_gindex, VertexHashTable ht_, Vertices v_, const float4 vc, const float4 vn)
{
    const int start_address = int(hash_function((uint)v_gindex) % ht_.t_size);
    // open hashing
    int h = start_address;
    int e = 1;
    for (int loop = 0; loop < 128; loop++) {
        int old = atomicCAS(&ht_.key[h], EMPTY_BUCKET_32, v_gindex);
        if (old == EMPTY_BUCKET_32) {
            const int a_ = atomicAdd(v_.t_size, 1);
            v_.vertices[a_] = vc;
            v_.normals[a_] = vn;
            ht_.addr[h] = a_;
            return h;
        }
        else if (v_gindex == old) {
            // vertex key already in table
            return h;
        }
        else {
            // step with linear probing
            h = (h + e*e) % ht_.t_size;
            e = e + 1;
            if (h == start_address) {
                printf("ERROR: can't find free bucket for %d\n", v_gindex);
                return -1;
            }
        }
    }

    return -1;
}

__device__ int insert_vertex_fast(const int v_gindex, VertexHashTable ht_, Vertices v_,int& address)
{
    const int start_address = int(hash_function((uint)v_gindex) % ht_.t_size);
    // open hashing
    int h = start_address;
    int e = 1;
    for (int loop = 0; loop < 128; loop++) {
        int old = atomicCAS(&ht_.key[h], EMPTY_BUCKET_32, v_gindex);
        if (old == EMPTY_BUCKET_32) {
            const int a_ = atomicAdd(v_.t_size, 1);
            //v_.vertices[a_] = vc;
            //v_.normals[a_] = vn;
            ht_.addr[h] = a_;
            address = a_;
            return h;
        }
        else if (v_gindex == old) {
            // vertex key already in table
            address = -1;
            return h;
        }
        else {
            // step with linear probing
            h = (h + e*e) % ht_.t_size;
            e = e + 1;
            if (h == start_address) {
                printf("ERROR: can't find free bucket for %d\n", v_gindex);
                return -1;
            }
        }
    }

    return -1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// find vertex global index in hash table
// Values were store at the hash address with open hashing and
__device__ int find_vertex(const int gl_index, VertexHashTable ht_)
{
	// compute hash for global index
	const int pos = int(hash_function((uint)gl_index) % ht_.t_size);

	// open hashing with quadratic probing
	int h = pos;
	int e = 1;
	for (int loop = 0; loop < 128; loop++) {
		if (ht_.key[h] == gl_index) {
			return ht_.addr[h];
		}
		else {
			// step with linear probing
			h = (h + e*e) % ht_.t_size;
			e = e + 1;
		}
	}
	printf("ERROR: can't find gl_index in hash table: gl_index %d at %d\n",gl_index, pos);
	return -1;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 64 bit hash table to use Cantor's pairing function
// 64 bit mix function
__device__ unsigned long long hash64shift(unsigned long long key)
{
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Insert the unique id of a halfedge obtained using the bijective Cantor's pairing function
// into the hast table. At the same position in array, save actuall have edge address, which
// will be used later to connect twin edges.
__device__ bool insert_halfedge_id(const int t_size, unsigned long long *he_table, int2* he_ids, int he_addr, int v0, int v1)
{
    //const unsigned long long EMPTY_BUCKET = 0ull;
    // compute pairing function value
    unsigned long long x = (unsigned long long)v0;
    unsigned long long y = (unsigned long long)v1;
    unsigned long long he_id = (x < y) ? y : x;
    he_id = he_id + (x + y) * (x + y + 1) / 2ull;
    // evalue hash function
    unsigned long long l_size = (unsigned long long)t_size;
    //unsigned long long he_id = (v0 < v1) ? (unsigned long long)v0 | ((unsigned long long)v1 << 32) : (unsigned long long)v1 | ((unsigned long long)v0 << 32);
    //const int start_address = int(hash64shift(he_id) % l_size);
    const int start_address = int( he_id % l_size );

    // open hashing
    int h = start_address;
    int e = 1;
    for (int loop = 0; loop < 128; loop++) {
        unsigned long long old = atomicCAS(&he_table[h], EMPTY_BUCKET_64, he_id);
        if (old == EMPTY_BUCKET_64 || old == he_id) {
            if (v0 < v1) {
                he_ids[h].x = he_addr;
            }
            else {
                he_ids[h].y = he_addr;
            }

            return true;
        }
        else {
            // step with linear probing
            h = (h + e*e) % t_size;
            e = e + 1;
            if (h == start_address) {
                printf("ERROR: he can't find free bucket for %d\n", he_id);
                return false;
            }
        }
    }

    return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute the cell vertices from uniform grid and cell indices
// Use spacing to compute vertex position
__device__ void cell_vertices(float3 v[8], const int i, const int j, const int k, UGrid ugrid)
{
	v[0].x = ugrid.x0 + i * ugrid.dx;
	v[0].y = ugrid.y0 + j * ugrid.dy;
	v[0].z = ugrid.z0 + k * ugrid.dz;

	v[1].x = v[0].x + ugrid.dx;
	v[1].y = v[0].y;
	v[1].z = v[0].z;

	v[2].x = v[0].x;
	v[2].y = v[0].y + ugrid.dy;
	v[2].z = v[0].z;

	v[3].x = v[0].x + ugrid.dx;
	v[3].y = v[0].y + ugrid.dy;
	v[3].z = v[0].z;

	v[4].x = v[0].x;
	v[4].y = v[0].y;
	v[4].z = v[0].z + ugrid.dz;

	v[5].x = v[0].x + ugrid.dx;
	v[5].y = v[0].y;
	v[5].z = v[0].z + ugrid.dz;

	v[6].x = v[0].x;
	v[6].y = v[0].y + ugrid.dy;
	v[6].z = v[0].z + ugrid.dz;

	v[7].x = v[0].x + ugrid.dx;
	v[7].y = v[0].y + ugrid.dy;
	v[7].z = v[0].z + ugrid.dz;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute gradient of scalar field at the vertices
// Use central differences, at the boundaries use forward
// or backward differences correspondigly
__device__ void gradient(float3 n[8], cudaTextureObject_t u, UGrid ugrid, const int i_index, const int j_index, const int k_index)
{
    const int idim = ugrid.idim;
    const int jdim = ugrid.jdim;
    const int kdim = ugrid.kdim;
    const float dx = ugrid.dx;
    const float dy = ugrid.dy;
    const float dz = ugrid.dz;

    int v0, v1;
    float f = 2.f;
    auto u_index = [](const int dim, int i, float& f) {
        f = (i<0) || (i >= dim) ? 1 : 2;
        i = (i<0) ? 0 : i;
        i = (i >= dim) ? dim - 1 : i;
        return i;
    };
    // 8 vertices
    // v0, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[0].x = (tex3D<float>(u, v1, j_index, k_index) - tex3D<float>(u, v0, j_index, k_index)) / (f * dx);
    // v0, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[0].y = (tex3D<float>(u, i_index, v1, k_index) - tex3D<float>(u, i_index, v0, k_index)) / (f * dy);
    // v0, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1;
    n[0].z = (tex3D<float>(u, i_index, j_index, v1) - tex3D<float>(u, i_index, j_index, v0)) / (f * dz);

    // v1, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[1].x = (tex3D<float>(u, v1, j_index, k_index) - tex3D<float>(u, v0, j_index, k_index)) / (f * dx);
    // v1, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[1].y = (tex3D<float>(u, i_index+1, v1, k_index) - tex3D<float>(u, i_index+1, v0, k_index)) / (f * dy);
    // v1, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1;
    n[1].z = (tex3D<float>(u, i_index+1, j_index, v1) - tex3D<float>(u, i_index+1, j_index, v0)) / (f * dz);

    // v2, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[2].x = (tex3D<float>(u, v1, j_index+1, k_index) - tex3D<float>(u, v0, j_index+1, k_index)) / (f * dx);
    // v2, y
    v0 = j_index;
    v1 = u_index(jdim,j_index + 2, f);
    n[2].y = (tex3D<float>(u, i_index, v1, k_index) - tex3D<float>(u, i_index, v0, k_index)) / (f * dy);
    // v2, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1, f;
    n[2].z = (tex3D<float>(u, i_index, j_index+1, v1) - tex3D<float>(u, i_index, j_index+1, v0)) / (f * dz);

    // v3, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[3].x = (tex3D<float>(u, v1, j_index+1, k_index) - tex3D<float>(u, v0, j_index+1, k_index)) / (f * dx);
    // v3, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[3].y = (tex3D<float>(u, i_index+1, v1, k_index) - tex3D<float>(u, i_index+1, v0, k_index)) / (f * dy);
    // v3, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1;
    n[3].z = (tex3D<float>(u, i_index+1, j_index+1, v1) - tex3D<float>(u, i_index+1, j_index+1, v0)) / (f * dz);

    // v4, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[4].x = (tex3D<float>(u, v1, j_index, k_index+1) - tex3D<float>(u, v0, j_index, k_index+1)) / (f * dx);
    // v4, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[4].y = (tex3D<float>(u, i_index, v1, k_index+1) - tex3D<float>(u, i_index, v0, k_index+1)) / (f * dy);
    // v4, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[4].z = (tex3D<float>(u, i_index, j_index, v1) - tex3D<float>(u, i_index, j_index, v0)) / (f * dz);

    // v5, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[5].x = (tex3D<float>(u, v1, j_index, k_index+1) - tex3D<float>(u, v0, j_index, k_index+1)) / (f * dx);
    // v5, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[5].y = (tex3D<float>(u, i_index+1, v1, k_index+1) - tex3D<float>(u, i_index+1, v0, k_index+1)) / (f * dy);
    // v5, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[5].z = (tex3D<float>(u, i_index+1, j_index, v1) - tex3D<float>(u, i_index+1, j_index, v0)) / (f * dz);

    // v6, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[6].x = (tex3D<float>(u, v1, j_index+1, k_index+1) - tex3D<float>(u, v0, j_index+1, k_index+1)) / (f * dx);
    // v6, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[6].y = (tex3D<float>(u, i_index, v1, k_index+1) - tex3D<float>(u, i_index, v0, k_index+1)) / (f * dy);
    // v6, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[6].z = (tex3D<float>(u, i_index, j_index+1, v1) - tex3D<float>(u, i_index, j_index+1, v0)) / (f * dz);

    // v7, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[7].x = (tex3D<float>(u, v1, j_index + 1, k_index + 1) - tex3D<float>(u, v0, j_index + 1, k_index + 1)) / (f * dx);
    // v7, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[7].y = (tex3D<float>(u, i_index + 1, v1, k_index + 1) - tex3D<float>(u, i_index + 1, v0, k_index + 1)) / (f * dy);
    // v7, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[7].z = (tex3D<float>(u, i_index+1, j_index+1, v1) - tex3D<float>(u, i_index + 1, j_index + 1, v0)) / (f * dz);
}

__device__ void gradientShared(const int tr, float3 n[8], cudaTextureObject_t u, UGrid ugrid, const int i_index, const int j_index, const int k_index)
{
    const int idim = ugrid.idim;
    const int jdim = ugrid.jdim;
    const int kdim = ugrid.kdim;
    const float dx = ugrid.dx;
    const float dy = ugrid.dy;
    const float dz = ugrid.dz;

    int v0, v1;
    float f = 2.f;
    auto u_index = [](const int dim, int i, float& f) {
        f = (i<0) || (i >= dim) ? 1 : 2;
        i = (i<0) ? 0 : i;
        i = (i >= dim) ? dim - 1 : i;
        return i;
    };
    // 8 vertices
    // v0, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[tr].x = (tex3D<float>(u, v1, j_index, k_index) - tex3D<float>(u, v0, j_index, k_index)) / (f * dx);
    // v0, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[tr].y = (tex3D<float>(u, i_index, v1, k_index) - tex3D<float>(u, i_index, v0, k_index)) / (f * dy);
    // v0, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1;
    n[tr].z = (tex3D<float>(u, i_index, j_index, v1) - tex3D<float>(u, i_index, j_index, v0)) / (f * dz);

    // v1, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[tr+1].x = (tex3D<float>(u, v1, j_index, k_index) - tex3D<float>(u, v0, j_index, k_index)) / (f * dx);
    // v1, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[tr + 1].y = (tex3D<float>(u, i_index + 1, v1, k_index) - tex3D<float>(u, i_index + 1, v0, k_index)) / (f * dy);
    // v1, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1;
    n[tr + 1].z = (tex3D<float>(u, i_index + 1, j_index, v1) - tex3D<float>(u, i_index + 1, j_index, v0)) / (f * dz);

    // v2, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[tr + 2].x = (tex3D<float>(u, v1, j_index + 1, k_index) - tex3D<float>(u, v0, j_index + 1, k_index)) / (f * dx);
    // v2, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[tr + 2].y = (tex3D<float>(u, i_index, v1, k_index) - tex3D<float>(u, i_index, v0, k_index)) / (f * dy);
    // v2, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1, f;
    n[tr + 2].z = (tex3D<float>(u, i_index, j_index + 1, v1) - tex3D<float>(u, i_index, j_index + 1, v0)) / (f * dz);

    // v3, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[tr + 3].x = (tex3D<float>(u, v1, j_index + 1, k_index) - tex3D<float>(u, v0, j_index + 1, k_index)) / (f * dx);
    // v3, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[tr + 3].y = (tex3D<float>(u, i_index + 1, v1, k_index) - tex3D<float>(u, i_index + 1, v0, k_index)) / (f * dy);
    // v3, z
    v0 = u_index(kdim, k_index - 1, f);
    v1 = k_index + 1;
    n[tr + 3].z = (tex3D<float>(u, i_index + 1, j_index + 1, v1) - tex3D<float>(u, i_index + 1, j_index + 1, v0)) / (f * dz);

    // v4, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[tr + 4].x = (tex3D<float>(u, v1, j_index, k_index + 1) - tex3D<float>(u, v0, j_index, k_index + 1)) / (f * dx);
    // v4, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[tr + 4].y = (tex3D<float>(u, i_index, v1, k_index + 1) - tex3D<float>(u, i_index, v0, k_index + 1)) / (f * dy);
    // v4, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[tr + 4].z = (tex3D<float>(u, i_index, j_index, v1) - tex3D<float>(u, i_index, j_index, v0)) / (f * dz);

    // v5, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[tr + 5].x = (tex3D<float>(u, v1, j_index, k_index + 1) - tex3D<float>(u, v0, j_index, k_index + 1)) / (f * dx);
    // v5, y
    v0 = u_index(jdim, j_index - 1, f);
    v1 = j_index + 1;
    n[tr + 5].y = (tex3D<float>(u, i_index + 1, v1, k_index + 1) - tex3D<float>(u, i_index + 1, v0, k_index + 1)) / (f * dy);
    // v5, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[tr + 5].z = (tex3D<float>(u, i_index + 1, j_index, v1) - tex3D<float>(u, i_index + 1, j_index, v0)) / (f * dz);

    // v6, x
    v0 = u_index(idim, i_index - 1, f);
    v1 = i_index + 1;
    n[tr + 6].x = (tex3D<float>(u, v1, j_index + 1, k_index + 1) - tex3D<float>(u, v0, j_index + 1, k_index + 1)) / (f * dx);
    // v6, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[tr + 6].y = (tex3D<float>(u, i_index, v1, k_index + 1) - tex3D<float>(u, i_index, v0, k_index + 1)) / (f * dy);
    // v6, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[tr + 6].z = (tex3D<float>(u, i_index, j_index + 1, v1) - tex3D<float>(u, i_index, j_index + 1, v0)) / (f * dz);

    // v7, x
    v0 = i_index;
    v1 = u_index(idim, i_index + 2, f);
    n[tr + 7].x = (tex3D<float>(u, v1, j_index + 1, k_index + 1) - tex3D<float>(u, v0, j_index + 1, k_index + 1)) / (f * dx);
    // v7, y
    v0 = j_index;
    v1 = u_index(jdim, j_index + 2, f);
    n[tr + 7].y = (tex3D<float>(u, i_index + 1, v1, k_index + 1) - tex3D<float>(u, i_index + 1, v0, k_index + 1)) / (f * dy);
    // v7, z
    v0 = k_index;
    v1 = u_index(kdim, k_index + 2, f);
    n[tr + 7].z = (tex3D<float>(u, i_index + 1, j_index + 1, v1) - tex3D<float>(u, i_index + 1, j_index + 1, v0)) / (f * dz);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// old fashined way to compute the gradient from a uniform grid
__device__ void gradient2(float3 n[8], cudaTextureObject_t u, const UGrid ugrid, const int i_index, const int j_index, const int k_index)
{
    for (int _k = 0; _k <= 1; _k++) {
        int k = k_index + _k;
        for (int _j = 0; _j <= 1; _j++) {
            int j = j_index + _j;
            for (int _i = 0; _i <= 1; _i++) {
                int i = i_index + _i;
                // set gradient at vertex
                unsigned int v_index{ 0 };
                v_index |= (_i) & 1;
                v_index |= (_j << 1) & 2;
                v_index |= (_k << 2) & 4;
                // x-component
                float factor{ 1.f };
                int i1{ 0 };
                int i2{ 0 };
                if (i == 0) {
                    i1 = i;
                    i2 = i + 1;
                }
                else if (i == ugrid.idim - 1) {
                    i1 = i - 1;
                    i2 = i;
                }
                else {
                    i1 = i - 1;
                    i2 = i + 1;
                    factor = 2.f;
                }
                n[v_index].x = (tex3D<float>(u, i2, j, k) - tex3D<float>(u, i1, j, k)) / (factor * ugrid.dx);

                // y-component
                factor = 1.f;
                if (j == 0) {
                    i1 = j;
                    i2 = j + 1;
                }
                else if (j == ugrid.jdim - 1) {
                    i1 = j - 1;
                    i2 = j;
                }
                else {
                    i1 = j - 1;
                    i2 = j + 1;
                    factor = 2.f;
                }
                n[v_index].y = (tex3D<float>(u, i, i2, k) - tex3D<float>(u, i, i1, k)) / (factor * ugrid.dy);

                // z-component
                factor = 1.f;
                if (k == 0) {
                    i1 = k;
                    i2 = k + 1;
                }
                else if (k == ugrid.kdim - 1) {
                    i1 = k - 1;
                    i2 = k;
                }
                else {
                    i1 = k - 1;
                    i2 = k + 1;
                    factor = 2.f;
                }
                n[v_index].z = (tex3D<float>(u, i, j, i2) - tex3D<float>(u, i, j, i1)) / (factor * ugrid.dz);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      CUDA GLOBAL FUNCTIONS
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table
__global__ void init_hash_table(VertexHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size)
        return;
    ht_.key[tid] = EMPTY_BUCKET_32;
    //ht_.addr[tid] = -1;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      TOPLOGICALLY CORRECT MARCHING CUBES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Marching Cubes
// Count cells being intersected by isosurface
// Coompute a  nr. of vertices
__global__ void mc_count(CellsIds cells, AmbiguousCells acells, int* v_count, int* a_count, float i0, cudaTextureObject_t v_data, UGrid ugrid, MC_lookup l_tables)
{
	const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    const int i = ugrid.i_index(gl_index);
    const int j = ugrid.j_index(gl_index);
    const int k = ugrid.k_index(gl_index);

	if (i >= ugrid.idim - 1 || j >= ugrid.jdim - 1 || k >= ugrid.kdim - 1) {
		return;
	}

	// scalar values at vertices
	float u[8];
	u[0] = tex3D<float>(v_data, i, j, k);
	u[1] = tex3D<float>(v_data, i + 1, j, k);
	u[2] = tex3D<float>(v_data, i, j + 1, k);
	u[3] = tex3D<float>(v_data, i + 1, j + 1, k);
	u[4] = tex3D<float>(v_data, i, j, k + 1);
	u[5] = tex3D<float>(v_data, i + 1, j, k + 1);
	u[6] = tex3D<float>(v_data, i, j + 1, k + 1);
	u[7] = tex3D<float>(v_data, i + 1, j + 1, k + 1);

	// compute case
	uint i_case{ 0 };
	i_case = i_case + ((uint)(u[0] >= i0));
	i_case = i_case + ((uint)(u[1] >= i0)) * 2;
	i_case = i_case + ((uint)(u[2] >= i0)) * 4;
	i_case = i_case + ((uint)(u[3] >= i0)) * 8;
	i_case = i_case + ((uint)(u[4] >= i0)) * 16;
	i_case = i_case + ((uint)(u[5] >= i0)) * 32;
	i_case = i_case + ((uint)(u[6] >= i0)) * 64;
	i_case = i_case + ((uint)(u[7] >= i0)) * 128;

	// compute number of vertices computed for this cell
    const ushort e_ = l_tables.e_[i_case];
    int nr_vertices = numberOfSetBits<ushort>(e_);
 
    // copy into global memory
	if (nr_vertices > 0) {
		// get an address
        if (e_ & BIT_16) {
            atomicAdd(a_count, nr_vertices);
            //acells.add(gl_index);
            acells.a_cells[atomicAdd(acells.t_size, 1)] = gl_index;
        }
        else {
            atomicAdd(v_count, nr_vertices);
            cells.cells_[atomicAdd(cells.t_size,1)] = gl_index;
        }
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Based on standard Marching Cubes compute cell triangulation for non ambiguous cases
// Save cell global id to process ambiguous cases later
// Parameters
// @param i0 is the iso-value
// @param v_data is a 3D texture with the volume data
// @param ugird contains information describing the uniform grid
// @param l_tables is a structure with pointers to the lookup table for MC
//        struct MC_lookup {
//            unsigned short* e_pattern;
//            int* t_pattern;
//            uint* t_ambig;
//		  };
// @param nr_cells is to the total nr.of cells intersected by the iso - surfece
// @param cellid is a field with the global id of the cells in the uniform grid
// @param c_addr is a atomic counter to compute the index in the array a_cells with ids of the cells which have an ambiguous case
// @param a_cell a point to a field containing the ids of the ambiguous cases to be processed later
// @param ht_ hash table to compute unique vertex index
// @param v_ a structure containing all pointer required for vertex processing
//        struct Vertices {
//            float4* vertices{ nullptr };
//            float4* normals{ nullptr };
//            int* t_size{ nullptr };
//         };
// @param t_  a structure containing all pointer required for triangle processing
//        struct Triangles {
//            int3* triangles{ nullptr };
//            int* t_size{ nullptr };
//        };
__global__ void mc_slice(const float i0, cudaTextureObject_t v_data, UGrid ugrid, MC_lookup l_tables, const int nr_cells, CellsIds cells, VertexHashTable ht_, Vertices v_, Triangles t_)
{
    __shared__ int4 tris[5 * MC_BLOCKSIZE];

	// get thread id
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int bz = blockDim.x;
    const int tr = threadIdx.x;
    if (nr_cells <= tid)
		return;
    
    for (int i = tr; i < 5 * bz; i += bz) {
        tris[i].x = -1;
    }
    //__syncthreads();

	// compute grid indices from global index
	const int gl_index = cells.cells_[tid];
    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
	
    // get address to store vertices
    //const int vtpos = cells.vtpos_[tid];

	// construct 8 cell vertices
	float3 v[8];
	cell_vertices(v, i_index, j_index, k_index, ugrid);

	// scalar values at vertices
	float u[8];
	u[0] = tex3D<float>(v_data, i_index, j_index, k_index);
	u[1] = tex3D<float>(v_data, i_index + 1, j_index, k_index);
	u[2] = tex3D<float>(v_data, i_index, j_index + 1, k_index);
	u[3] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index);
	u[4] = tex3D<float>(v_data, i_index, j_index, k_index + 1);
	u[5] = tex3D<float>(v_data, i_index + 1, j_index, k_index + 1);
	u[6] = tex3D<float>(v_data, i_index, j_index + 1, k_index + 1);
	u[7] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index + 1);

	// compute normals at vertices
	float3 n[8];
	gradient(n, v_data, ugrid, i_index, j_index, k_index);

	// compute case
	uint i_case{ 0 };
	i_case = i_case + ((uint)(u[0] >= i0));
	i_case = i_case + ((uint)(u[1] >= i0)) * 2;
	i_case = i_case + ((uint)(u[2] >= i0)) * 4;
	i_case = i_case + ((uint)(u[3] >= i0)) * 8;
	i_case = i_case + ((uint)(u[4] >= i0)) * 16;
	i_case = i_case + ((uint)(u[5] >= i0)) * 32;
	i_case = i_case + ((uint)(u[6] >= i0)) * 64;
	i_case = i_case + ((uint)(u[7] >= i0)) * 128;

	// ambiguous cases are processed in the next pass
    const ushort e_ = l_tables.e_[i_case];
	
	// Compute intersection with edges
	const unsigned long long gei_pattern_ = 670526590282893600ull;
	const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
	int v_gindex[12];
	ushort flag{ 1 };
    //const ushort e_pattern = (ushort)l_tables.ePattern(i_case); //  l_tables.e_pattern[i_case];
	for (int e = 0; e < 12; e++) {
		v_gindex[e] = -1;
		if (flag & e_) {
            // get unique vertex index
            // compute vertex global index
            const int ix = i_index + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
            const int iy = j_index + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
            const int iz = k_index + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
            const int off_val = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
            int address{ -1 };
            v_gindex[e] = insert_vertex_fast(int(9 * ugrid.gl_index(ix, iy, iz) + off_val), ht_, v_,address);
            if (address > -1) {
                // compute edge inersection
                // compute local coordinate along edge
                const int v0 = (l_edges_[e] & 0xF);
                const int v1 = (l_edges_[e] >> 4) & 0xF;
                const float l = (i0 - u[v0]) / (u[v1] - u[v0]);
                float4 vp = make_float4(v[v0].x + l*(v[v1].x - v[v0].x), v[v0].y + l*(v[v1].y - v[v0].y), v[v0].z + l*(v[v1].z - v[v0].z), 1.f);
                float4 np = make_float4(n[v0].x + l*(n[v1].x - n[v0].x), n[v0].y + l*(n[v1].y - n[v0].y), n[v0].z + l*(n[v1].z - n[v0].z), 0.f);
                const float length = std::sqrt(np.x * np.x + np.y * np.y + np.z * np.z);
                np.x = np.x / length;
                np.y = np.y / length;
                np.z = np.z / length;
                //np.w = 0.f;
                v_.vertices[address] = vp;
                v_.normals[address] = np;
            }
		}
		flag <<= 1;
	}

	// compute triangles
	//const unsigned char* t_ambig = l_tables.t_ambig;
	unsigned long long tl_ = l_tables.t_[i_case];
	for (int t = 0; t < 16; t += 3) {
        const int v0 = (int)((tl_ >> (4 * t)) & 0xFull); 
        //if (((tl_ >> (4*t))& 0xFull) == 0xF) {
        if (v0 == 0xF) {
			// there are no more triangles
			break;
		}
		// save tirangle
		//const int v0 = (int)((tl_ >> (4 * t)) & 0xFull);
        const int v1 = (int)((tl_ >> (4 * (t + 1))) & 0xFull);
        const int v2 = (int)((tl_ >> (4 * (t + 2))) & 0xFull);
        //t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[v0], v_gindex[v1], v_gindex[v2], 0);
        tris[tr + (t/3) * bz] = make_int4(v_gindex[v0], v_gindex[v1], v_gindex[v2], 0);
        //t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[((tl_ >> (4 * t)) & 0xFull)], v_gindex[((tl_ >> (4 * (t + 1))) & 0xFull)], v_gindex[((tl_ >> (4 * (t + 2))) & 0xFull)], 0);
	}
    __syncthreads();

    // write tris
    for (int i = tr; i < 5 * bz; i += bz) {
        if (tris[i].x > -1) {
            t_.triangles[atomicAdd(t_.t_size, 1)] = tris[i];
        }
        //__syncthreads();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute the triangulation of a cell with an ambiguous case
__global__ void t_slice(const float i0, cudaTextureObject_t v_data, UGrid ugrid, MC_lookup l_tables, const int nr_cells, AmbiguousCells acells, VertexHashTable ht_, Vertices v_, Triangles t_)
{
    __shared__ float3 n[8 * AMB_BLOCKSIZE];
    __shared__ float3 p[8 * AMB_BLOCKSIZE];
    // get cell id
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (nr_cells <= tid)
		return;
    
	// compute grid indices from global index
	const int gl_index = acells.a_cells[tid];
    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    const int tr = threadIdx.x * 8;

	// construct 8 cell vertices
	//float3 p[8];
	//cell_vertices(p, i_index, j_index, k_index, ugrid);
    p[tr].x = ugrid.x0 + i_index * ugrid.dx;
    p[tr].y = ugrid.y0 + j_index * ugrid.dy;
    p[tr].z = ugrid.z0 + k_index * ugrid.dz;

    p[tr + 1].x = p[tr].x + ugrid.dx;
    p[tr + 1].y = p[tr].y;
    p[tr + 1].z = p[tr].z;

    p[tr + 2].x = p[tr].x;
    p[tr + 2].y = p[tr].y + ugrid.dy;
    p[tr + 2].z = p[tr].z;

    p[tr + 3].x = p[tr].x + ugrid.dx;
    p[tr + 3].y = p[tr].y + ugrid.dy;
    p[tr + 3].z = p[tr].z;

    p[tr + 4].x = p[tr].x;
    p[tr + 4].y = p[tr].y;
    p[tr + 4].z = p[tr].z + ugrid.dz;

    p[tr + 5].x = p[tr].x + ugrid.dx;
    p[tr + 5].y = p[tr].y;
    p[tr + 5].z = p[tr].z + ugrid.dz;

    p[tr + 6].x = p[tr].x;
    p[tr + 6].y = p[tr].y + ugrid.dy;
    p[tr + 6].z = p[tr].z + ugrid.dz;

    p[tr + 7].x = p[tr].x + ugrid.dx;
    p[tr + 7].y = p[tr].y + ugrid.dy;
    p[tr + 7].z = p[tr].z + ugrid.dz;

	// scalar values at vertices
	float F[8];
	F[0] = tex3D<float>(v_data, i_index, j_index, k_index);
	F[1] = tex3D<float>(v_data, i_index + 1, j_index, k_index);
	F[2] = tex3D<float>(v_data, i_index, j_index + 1, k_index);
	F[3] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index);
	F[4] = tex3D<float>(v_data, i_index, j_index, k_index + 1);
	F[5] = tex3D<float>(v_data, i_index + 1, j_index, k_index + 1);
	F[6] = tex3D<float>(v_data, i_index, j_index + 1, k_index + 1);
	F[7] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index + 1);

	// compute normals at vertices
	//float3 n[8];
    gradientShared(tr, n, v_data, ugrid, i_index, j_index, k_index);
    __syncthreads();

	// compute case
	uint i_case{ 0 };
	i_case = i_case + ((uint)(F[0] >= i0));
	i_case = i_case + ((uint)(F[1] >= i0)) * 2;
	i_case = i_case + ((uint)(F[2] >= i0)) * 4;
	i_case = i_case + ((uint)(F[3] >= i0)) * 8;
	i_case = i_case + ((uint)(F[4] >= i0)) * 16;
	i_case = i_case + ((uint)(F[5] >= i0)) * 32;
	i_case = i_case + ((uint)(F[6] >= i0)) * 64;
	i_case = i_case + ((uint)(F[7] >= i0)) * 128;

	// Compute intersection with edges
	const unsigned long long gei_pattern_ = 670526590282893600ull;
	const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
	// compute intersection with cell edges
	float ecoord[12]{};
	int v_gindex[12]{};
	ushort flag{ 1 };
    ushort e_ = l_tables.e_[i_case];
    for (int e = 0; e < 12; e++) {
		v_gindex[e] = -1;
		//ecoord[e] = 0.f;
		if (flag & e_) {
            // get unique vertex index
            // compute vertex global index
            const int ix = i_index + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
            const int iy = j_index + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
            const int iz = k_index + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
            const int off_val = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
            int address{ -1 };
            v_gindex[e] = insert_vertex_fast(int(9 * ugrid.gl_index(ix, iy, iz) + off_val), ht_, v_,address);
            

            // compute edge inersection
            // compute local coordinate along edge
            const int v0 = (l_edges_[e] & 0xF);
            const int v1 = (l_edges_[e] >> 4) & 0xF;
            const float l = (i0 - F[v0]) / (F[v1] - F[v0]);
            if (address > -1) {
                v_.vertices[address] = make_float4(p[tr + v0].x + l*(p[tr + v1].x - p[tr + v0].x), p[tr + v0].y + l*(p[tr + v1].y - p[tr + v0].y), p[tr + v0].z + l*(p[tr + v1].z - p[tr + v0].z), 1.f);
                float4 np = make_float4(n[tr + v0].x + l*(n[tr + v1].x - n[tr + v0].x), n[tr + v0].y + l*(n[tr + v1].y - n[tr + v0].y), n[tr + v0].z + l*(n[tr + v1].z - n[v0].z), 0.f);
                const float length = std::sqrt(np.x * np.x + np.y * np.y + np.z * np.z);
                np.x = np.x / length;
                np.y = np.y / length;
                np.z = np.z / length;
                v_.normals[address] = np;
            }
            
            //v_gindex[e] =  v_.add(ht_, int(9 * ugrid.gl_index(ix, iy, iz) + off_val), vp, np);
			// remember local coordinate along edge
			ecoord[e] = l;
		}
		flag <<= 1;
	}

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
        return (f0*f3 - f1*f2) / (f0 + f3 - f1 - f2);
    };
    uchar f_flag{ 0 };
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
            if (val > i0) {
                set_segm(e3, e0, segm_);
                set_segm(e1, e2, segm_);
            }
            else if (val < i0) {
                set_segm(e1, e0, segm_);
                set_segm(e3, e2, segm_);
            }
            else {
                // set flag for this face
                f_flag |= (1 << f);
                float ec0 = ecoord[e0];
                float ec1 = ecoord[e1];
                float ec2 = ecoord[e2];
                float ec3 = ecoord[e3];
                if ((0x218 >> (f * 2)) & BIT_1) {
                    ec0 = 1 - ec0;
                    ec2 = 1 - ec2;
                }
                if ((0x218 >> (f * 2)) & BIT_2) {
                    ec1 = 1 - ec1;
                    ec3 = 1 - ec3;
                }
                if (ec1 < ec3 && ec0 > ec2) {
                    set_segm(e1, e0, segm_);
                    set_segm(e3, e2, segm_);
                }
                else if (ec1 > ec3 && ec0 < ec2) {
                    set_segm(e3, e0, segm_);
                    set_segm(e1, e2, segm_);
                }
                else {
                    return;
                }
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
            if (val > i0) {
                set_segm(e0, e1, segm_);
                set_segm(e2, e3, segm_);
            }
            else if (val < i0) {
                set_segm(e0, e3, segm_);
                set_segm(e2, e1, segm_);
            }
            else {
                f_flag = (1 << f);
                // singular case val == i0, there are no asymptotes
                // check if there is a reasonable triangulation of the face
                float ec0 = ecoord[e0];
                float ec1 = ecoord[e1];
                float ec2 = ecoord[e2];
                float ec3 = ecoord[e3];
                if ((0x218 >> (f * 2)) & BIT_1) {
                    ec0 = 1 - ec0;
                    ec2 = 1 - ec2;
                }
                if ((0x218 >> (f * 2)) & BIT_2) {
                    ec1 = 1 - ec1;
                    ec3 = 1 - ec3;
                }
                if (ec1 < ec3 && ec0 > ec2) {
                    set_segm(e0, e1, segm_);
                    set_segm(e2, e3, segm_);
                }
                else if (ec1 > ec3 && ec0 < ec2) {
                    set_segm(e0, e3, segm_);
                    set_segm(e2, e1, segm_);
                }
                else {
                    return;
                }
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
	float ui[2]{};
	float vi[2]{};
	float wi[2]{};
	unsigned char q_sol{ 0 };
	const float a = (F[0] - F[1])*(-F[6] + F[7] + F[4] - F[5]) - (F[4] - F[5])*(-F[2] + F[3] + F[0] - F[1]);
	const float b = (i0 - F[0])*(-F[6] + F[7] + F[4] - F[5]) + (F[0] - F[1])*(F[6] - F[4]) - (i0 - F[4])*(-F[2] + F[3] + F[0] - F[1]) - (F[4] - F[5])*(F[2] - F[0]);
	const float c = (i0 - F[0])*(F[6] - F[4]) - (i0 - F[4])*(F[2] - F[0]);;
	float d = b*b - 4 * a*c;
	if (d > 0) {
		d = std::sqrt(d);
		// compute u-coord of solutions
		ui[0] = (-b - d) / (2 * a);
		ui[1] = (-b + d) / (2 * a);
		// compute v-coord of solutions
		float g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
		float g2 = F[2] * (1 - ui[0]) + F[3] * ui[0];
		vi[0] = (i0 - g1) / (g2 - g1);
		if (isnan(vi[0]) || isinf(vi[0])) {
			vi[0] = -1.f;
		}
		g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
		g2 = F[2] * (1 - ui[1]) + F[3] * ui[1];
		vi[1] = (i0 - g1) / (g2 - g1);
		if (isnan(vi[1]) || isinf(vi[1])) {
			vi[1] = -1.f;
		}
		// compute w-coordinates of solutions
		g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
		g2 = F[4] * (1 - ui[0]) + F[5] * ui[0];
		wi[0] = (i0 - g1) / (g2 - g1);
		if (isnan(wi[0]) || isinf(wi[0])) {
			wi[0] = -1.f;
		}
		g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
		g2 = F[4] * (1 - ui[1]) + F[5] * ui[1];
		wi[1] = (i0 - g1) / (g2 - g1);
		if (isnan(wi[1]) || isinf(wi[1])) {
			wi[1] = -1.f;
		}

		// correct values for roots of quadratic equations
		// in case the asymptotic decider has failed
		if (f_flag & BIT_1) { // face 1, w = 0;
			if (wi[0] < wi[1]) wi[0] = 0;
			else wi[1] = 0;
		}
		if (f_flag & BIT_2) { // face 2, w = 1
			if (wi[0] > wi[1]) wi[1] = 1;
			else wi[1] = 1;
		}
		if (f_flag & BIT_3) { // face 3, v = 0
			if (vi[0] < vi[1]) vi[0] = 0;
			else vi[1] = 0;
		}
		if (f_flag & BIT_4) { // face 4, v = 1
			if (vi[0] > vi[1]) vi[0] = 1;
			else vi[1] = 1;
		}
		if (f_flag & BIT_5) { // face 5, u = 0
			if (ui[0] < ui[1]) ui[0] = 0;
			else ui[1] = 0;
		}
		if (f_flag & BIT_6) { // face 6, u = 1
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
	// Triangles are stored in global memory starting at offset
	// count nr. of inner vertices to compute right global index
	// first inner vertex has index cell_global_index + 3;
	int v_count{ 3 };
	if (numberOfSetBits<unsigned char>(q_sol) == 6) {
		// there are at most three contours
		// Possible cases:
		//  1) a single contour with 12 vertices
		//  2) two contours which build a tunnel
		//  3) three contours, one has only 3 vertices and does not belong to the tunnel

		// construct the six vertices of the inner hexagon
		float3 hvt[6];
		hvt[0].x = ui[0]; hvt[0].y = vi[0]; hvt[0].z = wi[0];
		hvt[1].x = ui[0]; hvt[1].y = vi[0]; hvt[1].z = wi[1];
		hvt[2].x = ui[1]; hvt[2].y = vi[0]; hvt[2].z = wi[1];
		hvt[3].x = ui[1]; hvt[3].y = vi[1]; hvt[3].z = wi[1];
		hvt[4].x = ui[1]; hvt[4].y = vi[1]; hvt[4].z = wi[0];
		hvt[5].x = ui[0]; hvt[5].y = vi[1]; hvt[5].z = wi[0];

		// construct vertices at intersections with the edges
		auto e_vert = [&ecoord](const int e, const int i) {
			const unsigned int l_coord[3]{ 1324855, 5299420, 16733440 };
			unsigned char flag = (l_coord[i] >> (2 * e)) & 3;
			if (flag == 3)
				return ecoord[e];
			else
				return (float)(flag);

		};

		// if there are three contours, then there is a tunnel and one
		// of the contours is not part of it.
		unsigned char _not_tunnel = 0xF;
		if (cnt_ == 3) {
			// loop over the contorus
			// triangulate the contour which is not part of
			// the tunnel
			const float uc_min = (ui[0] < ui[1]) ? ui[0] : ui[1];
			const float uc_max = (ui[0] < ui[1]) ? ui[1] : ui[0];
			for (int t = 0; t < (int)cnt_; t++) {
				if (get_cnt_size(t, c_) == 3) {
					float umin = 2;
					float umax = -2;
					uint e0 = get_c(t, 0, c_);
					uint e1 = get_c(t, 1, c_);
					uint e2 = get_c(t, 2, c_);
					const float u_e0 = e_vert(e0, 0);
					const float u_e1 = e_vert(e1, 0);
					const float u_e2 = e_vert(e2, 0);
					umin = (u_e0 < umin) ? u_e0 : umin;
					umin = (u_e1 < umin) ? u_e1 : umin;
					umin = (u_e2 < umin) ? u_e2 : umin;
					umax = (u_e0 > umax) ? u_e0 : umax;
					umax = (u_e1 > umax) ? u_e1 : umax;
					umax = (u_e2 > umax) ? u_e1 : umax;
					if (uc_min > umax || uc_max < umin) {
						// this contour is not part of the tunnel
						_not_tunnel = t;
						// save triangle in global memory
						t_.triangles[atomicAdd(t_.t_size,1)] = make_int4(v_gindex[e0], v_gindex[e1], v_gindex[e2], 0);
					}
				}
			}
		}

		// compute vertices of inner hexagon, save new vertices in list and compute and keep
		// global vertice index to build triangle connectivity later on.
		int tg_idx[6];
        float4 po;
		for (int i = 0; i < 6; i++) {
            int address{ -1 };
            tg_idx[i] = insert_vertex_fast(int(9 * gl_index + v_count), ht_, v_, address);
            // update nr. of vertices
            v_count++;
            // create a store vertex and normal
			//float4 po;
			//float4 hn;
			// local coordinates for trilinear interpolation
			const float u = hvt[i].x; const float v = hvt[i].y; const float w = hvt[i].z;
			po.x = (1 - w)*((1 - v)*(p[tr].x + u*(p[tr + 1].x - p[tr].x)) + v*(p[tr + 2].x + u*(p[tr + 3].x - p[tr + 2].x))) + w*((1 - v)*(p[tr + 4].x + u*(p[tr + 5].x - p[tr + 4].x)) + v*(p[tr + 6].x + u*(p[tr + 7].x - p[tr + 6].x)));
		    po.y = (1 - w)*((1 - v)*(p[tr].y + u*(p[tr + 1].y - p[tr].y)) + v*(p[tr + 2].y + u*(p[tr + 3].y - p[tr + 2].y))) + w*((1 - v)*(p[tr + 4].y + u*(p[tr + 5].y - p[tr + 4].y)) + v*(p[tr + 6].y + u*(p[tr + 7].y - p[tr + 6].y)));
			po.z = (1 - w)*((1 - v)*(p[tr].z + u*(p[tr + 1].z - p[tr].z)) + v*(p[tr + 2].z + u*(p[tr + 3].z - p[tr + 2].z))) + w*((1 - v)*(p[tr + 4].z + u*(p[tr + 5].z - p[tr + 4].z)) + v*(p[tr + 6].z + u*(p[tr + 7].z - p[tr + 6].z)));
            //trilinear(po, p, hvt[i].x, hvt[i].y, hvt[i].z);
            po.w = 1.f;
            v_.vertices[address] = po;
            //trilinear(po, n, hvt[i].x, hvt[i].y, hvt[i].z);
            po.x = (1 - w)*((1 - v)*(n[tr].x + u*(n[tr + 1].x - n[tr].x)) + v*(n[tr + 2].x + u*(n[tr + 3].x - n[tr + 2].x))) + w*((1 - v)*(n[tr + 4].x + u*(n[tr + 5].x - n[tr + 4].x)) + v*(n[tr + 6].x + u*(n[tr + 7].x - n[tr + 6].x)));
            po.y = (1 - w)*((1 - v)*(n[tr].y + u*(n[tr + 1].y - n[tr].y)) + v*(n[tr + 2].y + u*(n[tr + 3].y - n[tr + 2].y))) + w*((1 - v)*(n[tr + 4].y + u*(n[tr + 5].y - n[tr + 4].y)) + v*(n[tr + 6].y + u*(n[tr + 7].y - n[tr + 6].y)));
            po.z = (1 - w)*((1 - v)*(n[tr].z + u*(n[tr + 1].z - n[tr].z)) + v*(n[tr + 2].z + u*(n[tr + 3].z - n[tr + 2].z))) + w*((1 - v)*(n[tr + 4].z + u*(n[tr + 5].z - n[tr + 4].z)) + v*(n[tr + 6].z + u*(n[tr + 7].z - n[tr + 6].z)));
            // normalize normal
			const float factor = std::sqrt(po.x * po.x + po.y * po.y + po.z * po.z);
            po.x = po.x / factor;
            po.y = po.y / factor;
            po.z = po.z / factor;
            po.w = 0.f;
            v_.normals[address] = po;
		}


		// triangulate contours with inner hexagon
		unsigned char tcon_[12];
		for (int i = 0; i < (int)cnt_; i++) {
			if (_not_tunnel != i) { // contour belongs to tunnel
				const int cnt_sz = (int)get_cnt_size(i, c_);
				for (int r = 0; r < cnt_sz; r++) {
					int index = -1;
					double dist = 1000.;
					uint ci = get_c(i, r, c_);
					const float u_edge = e_vert(ci, 0);
					const float v_edge = e_vert(ci, 1);
					const float w_edge = e_vert(ci, 2);
					for (int s = 0; s < 6; s++) {
						const float uval = u_edge - hvt[s].x;
						const float vval = v_edge - hvt[s].y;
						const float wval = w_edge - hvt[s].z;
						const float val = uval*uval + vval*vval + wval*wval;
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
                        t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[tid1], v_gindex[tid2], tg_idx[cid1],0);
					}
					break;
					case 1:
					{
						// measure diagonals
						// triangulate along shortest diagonal
						float u_edge = e_vert(tid1, 0);
						float v_edge = e_vert(tid1, 1);
						float w_edge = e_vert(tid1, 2);
						const float l1 = (u_edge - hvt[cid2].x)*(u_edge - hvt[cid2].x) + (v_edge - hvt[cid2].y)*(v_edge - hvt[cid2].y) + (w_edge - hvt[cid2].z)*(w_edge - hvt[cid2].z);
						u_edge = e_vert(tid2, 0);
						v_edge = e_vert(tid2, 1);
						w_edge = e_vert(tid2, 2);
						const double l2 = (u_edge - hvt[cid1].x)*(u_edge - hvt[cid1].x) + (v_edge - hvt[cid1].y)*(v_edge - hvt[cid1].y) + (w_edge - hvt[cid1].z)*(w_edge - hvt[cid1].z);
						const int a_ = atomicAdd(t_.t_size, 2);
						if (l1 < l2) {
                            t_.triangles[a_]   = make_int4(v_gindex[tid1], v_gindex[tid2], tg_idx[cid2], 0);
                            t_.triangles[a_+1] = make_int4(v_gindex[tid1], tg_idx[cid2], tg_idx[cid1], 0);
						}
						else {
                            t_.triangles[a_]     = make_int4(v_gindex[tid1], v_gindex[tid2], tg_idx[cid1], 0);
                            t_.triangles[a_ + 1] = make_int4(v_gindex[tid2], tg_idx[cid2], tg_idx[cid1], 0);
						}
					}
					break;
					case 2:
					{
						const int cidm = midpointRingIntModulo(cid1, cid2);
                        const int a_ = atomicAdd(t_.t_size, 3);
                        t_.triangles[a_]   = make_int4(v_gindex[tid1], v_gindex[tid2], tg_idx[cidm], 0);
                        t_.triangles[a_+1] = make_int4(v_gindex[tid1], tg_idx[cidm], tg_idx[cid1], 0);
                        t_.triangles[a_+2] = make_int4(v_gindex[tid2], tg_idx[cid2], tg_idx[cidm], 0);
					}
					break;
					} // switch
				} // for loop over the vertices of the contour
			} // if (_not_tunnel)
		} // for loop over contours
		if (cnt_ == 1) {
			// there is a single contour
			// triangulate and close inner hexagon
            const int a_ = atomicAdd(t_.t_size, 4);
            const bool s_ = asymptotic_decider(F[0], F[1], F[2], F[3]);
            const bool of_ = (wi[1] < wi[0]) ? s_ : !s_;
            if (!of_) {
                t_.triangles[a_] = make_int4(tg_idx[0], tg_idx[2], tg_idx[1], 0);
                t_.triangles[a_ + 1] = make_int4(tg_idx[2], tg_idx[4], tg_idx[3], 0);
                t_.triangles[a_ + 2] = make_int4(tg_idx[0], tg_idx[5], tg_idx[4], 0);
                t_.triangles[a_ + 3] = make_int4(tg_idx[0], tg_idx[4], tg_idx[2], 0);
            }
            else {
                t_.triangles[a_] = make_int4(tg_idx[0], tg_idx[1], tg_idx[2], 0);
                t_.triangles[a_ + 1] = make_int4(tg_idx[2], tg_idx[3], tg_idx[4], 0);
                t_.triangles[a_ + 2] = make_int4(tg_idx[0], tg_idx[4], tg_idx[5], 0);
                t_.triangles[a_ + 3] = make_int4(tg_idx[0], tg_idx[2], tg_idx[4], 0);
            }
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
                    //const int a_ = atomicAdd(t_.t_size, 1);
                    t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)], 0);
				}
				break;
				case 4:
				{
                    const int a_ = atomicAdd(t_.t_size, 2);
                    t_.triangles[a_]   = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)], 0);
                    t_.triangles[a_+1] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 2, c_)], v_gindex[get_c(i, 3, c_)], 0);
				}
				break;
				case 5:
				{
                    const int a_ = atomicAdd(t_.t_size, 3);
                    t_.triangles[a_]   = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)], 0);
                    t_.triangles[a_+1] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 2, c_)], v_gindex[get_c(i, 3, c_)], 0);
                    t_.triangles[a_+2] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 3, c_)], v_gindex[get_c(i, 4, c_)], 0);
				}
				break;
				case 6:
				{
                    const int a_ = atomicAdd(t_.t_size, 4);
                    t_.triangles[a_]   = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 3, c_)], 0);
                    t_.triangles[a_+1] = make_int4(v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)], v_gindex[get_c(i, 3, c_)], 0);
                    t_.triangles[a_+2] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 3, c_)], v_gindex[get_c(i, 4, c_)], 0);
                    t_.triangles[a_+3] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 4, c_)], v_gindex[get_c(i, 5, c_)], 0);
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
			unsigned char fs[3][2]{{(unsigned char)(q_sol & 1), (unsigned char)((q_sol >> 1) & 1)}, { (unsigned char)((q_sol >> 2) & 1), (unsigned char)((q_sol >> 3) & 1) }, { (unsigned char)((q_sol >> 4) & 1), (unsigned char)((q_sol >> 5) & 1) }};

			const unsigned char fc1 = fs[0][0] * fs[1][0] + fs[0][1] * fs[1][1];
			const unsigned char fc2 = fs[0][0] * fs[2][0] + fs[0][1] * fs[2][1];
			const unsigned char fc3 = fs[1][0] * fs[2][1] + fs[1][1] * fs[2][0];
			const unsigned char c_faces = fc1 + fc2 + fc3;
			float ucoord{};
			float vcoord{};
			float wcoord{};
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
			float4 ip;
			float4 in;
			//ip.x = (1 - wcoord)*((1 - vcoord)*(p[0].x + ucoord*(p[1].x - p[0].x)) + vcoord*(p[2].x + ucoord*(p[3].x - p[2].x))) + wcoord*((1 - vcoord)*(p[4].x + ucoord*(p[5].x - p[4].x)) + vcoord*(p[6].x + ucoord*(p[7].x - p[6].x)));
			//ip.y = (1 - wcoord)*((1 - vcoord)*(p[0].y + ucoord*(p[1].y - p[0].y)) + vcoord*(p[2].y + ucoord*(p[3].y - p[2].y))) + wcoord*((1 - vcoord)*(p[4].y + ucoord*(p[5].y - p[4].y)) + vcoord*(p[6].y + ucoord*(p[7].y - p[6].y)));
			//ip.z = (1 - wcoord)*((1 - vcoord)*(p[0].z + ucoord*(p[1].z - p[0].z)) + vcoord*(p[2].z + ucoord*(p[3].z - p[2].z))) + wcoord*((1 - vcoord)*(p[4].z + ucoord*(p[5].z - p[4].z)) + vcoord*(p[6].z + ucoord*(p[7].z - p[6].z)));
			//in.x = (1 - wcoord)*((1 - vcoord)*(n[0].x + ucoord*(n[1].x - n[0].x)) + vcoord*(n[2].x + ucoord*(n[3].x - n[2].x))) + wcoord*((1 - vcoord)*(n[4].x + ucoord*(n[5].x - n[4].x)) + vcoord*(n[6].x + ucoord*(n[7].x - n[6].x)));
			//in.y = (1 - wcoord)*((1 - vcoord)*(n[0].y + ucoord*(n[1].y - n[0].y)) + vcoord*(n[2].y + ucoord*(n[3].y - n[2].y))) + wcoord*((1 - vcoord)*(n[4].y + ucoord*(n[5].y - n[4].y)) + vcoord*(n[6].y + ucoord*(n[7].y - n[6].y)));
			//in.z = (1 - wcoord)*((1 - vcoord)*(n[0].z + ucoord*(n[1].z - n[0].z)) + vcoord*(n[2].z + ucoord*(n[3].z - n[2].z))) + wcoord*((1 - vcoord)*(n[4].z + ucoord*(n[5].z - n[4].z)) + vcoord*(n[6].z + ucoord*(n[7].z - n[6].z)));
            trilinear(ip, p, ucoord, vcoord, wcoord);
            trilinear(in, n, ucoord, vcoord, wcoord);
			// normalize normal
			const float factor = std::sqrt(in.x * in.x + in.y * in.y + in.z * in.z);
			in.x = in.x / factor;
			in.y = in.y / factor;
			in.z = in.z / factor;
            // the fourth coordinate
            ip.w = 1.f;
            in.w = 0.f;
            // global index
			//const int gidx = int(9 * gl_index + v_count);
            int gidx = int(9 * gl_index + v_count);
			// this point is only used if contours with more than three vertices
			// are present
			//bool pt_used{ false };

            // check if the vertex will be used, this happens
            // if there are contours with more than three edges
            for (int i = 0; i < (int)cnt_; i++) {
                if (get_cnt_size(i, c_) > 3) {
                    int address{ -1 };
                    gidx = insert_vertex_fast(gidx, ht_, v_, address);
                    //gidx = v_.add(ht_, gidx, ip, in);
                    v_count++;
                    float4 ip;
                    //trilinear(ip, p, ucoord, vcoord, wcoord);
                    ip.x = (1 - wcoord)*((1 - vcoord)*(p[tr].x + ucoord*(p[tr + 1].x - p[tr].x)) + vcoord*(p[tr + 2].x + ucoord*(p[tr + 3].x - p[tr + 2].x))) + wcoord*((1 - vcoord)*(p[tr + 4].x + ucoord*(p[tr + 5].x - p[tr + 4].x)) + vcoord*(p[tr + 6].x + ucoord*(p[tr + 7].x - p[tr + 6].x)));
                    ip.y = (1 - wcoord)*((1 - vcoord)*(p[tr].y + ucoord*(p[tr + 1].y - p[tr].y)) + vcoord*(p[tr + 2].y + ucoord*(p[tr + 3].y - p[tr + 2].y))) + wcoord*((1 - vcoord)*(p[tr + 4].y + ucoord*(p[tr + 5].y - p[tr + 4].y)) + vcoord*(p[tr + 6].y + ucoord*(p[tr + 7].y - p[tr + 6].y)));
                    ip.z = (1 - wcoord)*((1 - vcoord)*(p[tr].z + ucoord*(p[tr + 1].z - p[tr].z)) + vcoord*(p[tr + 2].z + ucoord*(p[tr + 3].z - p[tr + 2].z))) + wcoord*((1 - vcoord)*(p[tr + 4].z + ucoord*(p[tr + 5].z - p[tr + 4].z)) + vcoord*(p[tr + 6].z + ucoord*(p[tr + 7].z - p[tr + 6].z)));
                    ip.w = 1.f;
                    v_.vertices[address] = ip;
                    //trilinear(ip, n, ucoord, vcoord, wcoord);
                    ip.x = (1 - wcoord)*((1 - vcoord)*(n[tr].x + ucoord*(n[tr + 1].x - n[tr].x)) + vcoord*(n[tr + 2].x + ucoord*(n[tr + 3].x - n[tr + 2].x))) + wcoord*((1 - vcoord)*(n[tr + 4].x + ucoord*(n[tr + 5].x - n[tr + 4].x)) + vcoord*(n[tr + 6].x + ucoord*(n[tr + 7].x - n[tr + 6].x)));
                    ip.y = (1 - wcoord)*((1 - vcoord)*(n[tr].y + ucoord*(n[tr + 1].y - n[tr].y)) + vcoord*(n[tr + 2].y + ucoord*(n[tr + 3].y - n[tr + 2].y))) + wcoord*((1 - vcoord)*(n[tr + 4].y + ucoord*(n[tr + 5].y - n[tr + 4].y)) + vcoord*(n[tr + 6].y + ucoord*(n[tr + 7].y - n[tr + 6].y)));
                    ip.z = (1 - wcoord)*((1 - vcoord)*(n[tr].z + ucoord*(n[tr + 1].z - n[tr].z)) + vcoord*(n[tr + 2].z + ucoord*(n[tr + 3].z - n[tr + 2].z))) + wcoord*((1 - vcoord)*(n[tr + 4].z + ucoord*(n[tr + 5].z - n[tr + 4].z)) + vcoord*(n[tr + 6].z + ucoord*(n[tr + 7].z - n[tr + 6].z)));
                    // normalize normal
                    const float factor = std::sqrt(ip.x * ip.x + ip.y * ip.y + ip.z * ip.z);
                    ip.x = ip.x / factor;
                    ip.y = ip.y / factor;
                    ip.z = ip.z / factor;
                    ip.w = 0.f;
                    v_.normals[address] = ip;

                    break;
                }
            }


			// loop over the contorus
			for (int i = 0; i < (int)cnt_; i++) {
				const unsigned char cnt_sz = (unsigned char)get_cnt_size(i, c_);
				if (cnt_sz == 3) {
                    //const int a_ = atomicAdd(t_.t_size, 1);
                    t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)], 0);
				}
				else {
					//pt_used = true;
					for (int t = 0; t < cnt_sz; t++) {
						// add triangle to list
                        t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[get_c(i, t, c_)], v_gindex[get_c(i, (t + 1) % cnt_sz, c_)], gidx, 0);
					}
				}
			}
		} // else - there are saddle points
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id to position in vertex array
// to construct shared vertex list
__global__ void map_triangles(const int nr_t, VertexHashTable ht_, Triangles t_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= nr_t)
		return;

	t_.triangles[tid].x = ht_.addr[t_.triangles[tid].x]; //find_vertex(k0, ht_);
    t_.triangles[tid].y = ht_.addr[t_.triangles[tid].y]; //find_vertex(k1, ht_);
    t_.triangles[tid].z = ht_.addr[t_.triangles[tid].z]; //find_vertex(k2, ht_);
}

__global__ void map_triangles_fast(const int nr_t, VertexHashTable ht_, Triangles t_)
{
    const int tid = (blockIdx.x * blockDim.x + threadIdx.x);
    const int offset = tid % 3;
    const int t = tid / 3;
    if (t >= nr_t)
        return;

    switch (offset) {
    case 0:
        t_.triangles[t].x = ht_.addr[t_.triangles[t].x]; //find_vertex(k0, ht_);
        break;
    case 1:
        t_.triangles[t].y = ht_.addr[t_.triangles[t].y]; //find_vertex(k1, ht_);
        break;
    case 2:
        t_.triangles[t].z = ht_.addr[t_.triangles[t].z]; //find_vertex(k2, ht_);
        break;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// create three halfedges from triangle
// for the triangle store a corresponding halfe edge
// for each vertex store the starting halfedge
__global__ void create_halfedge(const int nr_he, Triangles t_, Halfedges he_, HalfedgeHashTable het_) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_he) {
        return;
    }

    const int offset = tid % 3;
    const int t = tid / 3;

    // get data from global memory
    const int4 tri = t_.triangles[tid / 3];
    // create and save halfedges
    // halfedge
    switch (offset) {
    case 0:
        he_.he_e[tid].x = tri.x; // origin vertex
        he_.he_e[tid].y = t; // face id
        he_.he_e[tid].z = tid + 1; // next halfedge
        he_.he_e[tid].w = -1; // boundary edge at initialization
        he_.he_v[tri.x] = tid;
        break;
    case 1:
        he_.he_e[tid].x = tri.y; // origin vertex
        he_.he_e[tid].y = t; // face id
        he_.he_e[tid].z = tid + 1; // next halfedge
        he_.he_e[tid].w = -1; // boundary edge at initialization
        he_.he_v[tri.y] = tid;
        break;
    case 2:
        he_.he_e[tid].x = tri.z; // origin vertex
        he_.he_e[tid].y = t; // face id
        he_.he_e[tid].z = tid - 2; // next halfedge
        he_.he_e[tid].w = -1; // boundary edge at initialization
        he_.he_v[tri.z] = tid;
        break;
    }
    // save halfedge ids in hash table to construct later twins neighborhood
    // insert_halfedge_id(const int t_size, unsigned long long *he_table, int2* he_ids, int he_addr, int v0, int v1)
    //insert_halfedge_id(he_size, he_table, he_ids, 3 * tid, tri.x, tri.y); 
    //insert_halfedge_id(he_size, he_table, he_ids, 3 * tid + 1, tri.y, tri.z);
    //insert_halfedge_id(he_size, he_table, he_ids, 3 * tid + 2, tri.z, tri.x);
   addHalfedgeToHashTable(het_, tid, tri.x, tri.y);
   //addHalfedgeToHashTable(het_,3 * tid + 1, tri.y, tri.z);
   //addHalfedgeToHashTable(het_,3 * tid + 2, tri.z, tri.x);

    // map vertex to halfedge
    //he_.he_v[tri.x] = 3 * tid;
    //he_.he_v[tri.y] = 3 * tid + 1;
    //he_.he_v[tri.z] = 3 * tid + 2;
    // map face to halfedge, this is a redundant information
    //he_f[tid] = 3 * tid;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// map vertex index of src vertex in halfedge indices to final index in vertex array
// at this point we know the total nr. of half edges, and the total nr. of vertices
__global__ void map_halfedge_vertex(const int nr_he, int4* he_e, VertexHashTable ht_, int* he_v)
{
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (gl_index >= nr_he)
        return;
    // he.x = origin vertex
    // he.y = face
    // he.z = next
    // he.w = tween

    // set vertex id
    //const int v = find_vertex(he_e[gl_index].x, ht_);
    const int v = ht_.addr[he_e[gl_index].x];
    he_e[gl_index].x = v;
    he_v[v] = gl_index; // vertex points to this halfedges
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// map index to boundary vertex
__global__ void map_halfedge_bndvertex(const int nr_he, int4* he_e, int* he_v)
{
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (gl_index >= nr_he)
        return;
    // he.x = origin vertex
    // he.y = face
    // he.z = next
    // he.w = tween, -1 if boundary edge

    if (he_e[gl_index].w == -1) {
        he_v[he_e[gl_index].x] = gl_index;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// map halfedge twins
__global__ void map_halfedge_twins(Halfedges he_, HalfedgeHashTable het_)
{
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (gl_index >= het_.t_size) {
        return;
    }
    if (het_.key[gl_index] > 0) {
        const int he0 = het_.he_ids[gl_index].x;
        const int he1 = het_.he_ids[gl_index].y;
        he_.he_e[he0].w = he1;
        he_.he_e[he1].w = he0;
    }

}

__global__ void map_halfedge_twins_fast(Halfedges he_, HalfedgeHashTable het_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int gl_index = tid / 2;
    const int offset = tid % 2;
    if (gl_index >= het_.t_size) {
        return;
    }
    if (het_.key[gl_index] > 0) {
        const int he0 = het_.he_ids[gl_index].x;
        const int he1 = het_.he_ids[gl_index].y;
        switch (offset) {
        case 0:
            he_.he_e[he0].w = he1;
            break;
        case 1:
            he_.he_e[he1].w = he0;
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// add halfedge to table
__device__ void addHalfedges(Halfedges he_, HalfedgeHashTable het_, const int v0, const int v1, const int v2)
{
    const int a_ = atomicAdd(he_.t_size, 1);
    const int f_ = a_ / 3;
    // he 0
    he_.he_e[a_].x = v0;
    he_.he_e[a_].y = f_; // there are three halfedges for each face (triangle)
    he_.he_e[a_].z = a_ + 1; // next
    he_.he_e[a_].w = -1; // default is boundary edge
   addHalfedgeToHashTable(het_,a_, v0, v1);

    // he 1
    he_.he_e[a_ + 1].x = v1;
    he_.he_e[a_ + 1].y = f_; // there are three halfedges for each face (triangle)
    he_.he_e[a_ + 1].z = a_ + 2;
    he_.he_e[a_ + 1].w = -1; // default is boundary edge
   addHalfedgeToHashTable(het_,a_ + 1, v1, v2);

    // he 2
    he_.he_e[a_ + 2].x = v2;
    he_.he_e[a_ + 2].y = f_; // there are three halfedges for each face (triangle)
    he_.he_e[a_ + 2].z = a_;
    he_.he_e[a_ + 2].w = -1; // default is boundary edge
   addHalfedgeToHashTable(het_,a_ + 2, v2, v0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Halfedge Marching cubes
__global__ void he_mcSlice(const float i0, cudaTextureObject_t v_data, UGrid ugrid, MC_lookup l_tables, int nr_cells, int* cellid, AmbiguousCells ac_, VertexHashTable ht_, Vertices v_, HalfedgeHashTable het_, Halfedges he_)
{
    // get thread id
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nr_cells <= tid)
        return;

    // compute grid indices from global index
    const int gl_index = cellid[tid];
    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    
    // construct 8 cell vertices
    float3 v[8];
    cell_vertices(v, i_index, j_index, k_index, ugrid);

    // scalar values at vertices
    float u[8];
    u[0] = tex3D<float>(v_data, i_index, j_index, k_index);
    u[1] = tex3D<float>(v_data, i_index + 1, j_index, k_index);
    u[2] = tex3D<float>(v_data, i_index, j_index + 1, k_index);
    u[3] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index);
    u[4] = tex3D<float>(v_data, i_index, j_index, k_index + 1);
    u[5] = tex3D<float>(v_data, i_index + 1, j_index, k_index + 1);
    u[6] = tex3D<float>(v_data, i_index, j_index + 1, k_index + 1);
    u[7] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index + 1);

    // compute normals at vertices
    float3 n[8];
    gradient(n, v_data, ugrid, i_index, j_index, k_index);

    // compute case
    uint i_case{ 0 };
    i_case = i_case + ((uint)(u[0] >= i0));
    i_case = i_case + ((uint)(u[1] >= i0)) * 2;
    i_case = i_case + ((uint)(u[2] >= i0)) * 4;
    i_case = i_case + ((uint)(u[3] >= i0)) * 8;
    i_case = i_case + ((uint)(u[4] >= i0)) * 16;
    i_case = i_case + ((uint)(u[5] >= i0)) * 32;
    i_case = i_case + ((uint)(u[6] >= i0)) * 64;
    i_case = i_case + ((uint)(u[7] >= i0)) * 128;

    // ambiguous cases are processed in the next pass
    //if (105 == l_tables.t_ambig[i_case]) {
    ushort e_ = l_tables.e_[i_case];
    if (e_ & BIT_16) {
        //ac_.add(gl_index);
        ac_.a_cells[atomicAdd(ac_.t_size, 1)] = gl_index;
        return; // don't process this cell with standard MC
    }

    // Compute intersection with edges
    const unsigned long long gei_pattern_ = 670526590282893600ull;
    const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
    int v_gindex[12]{};
    ushort flag{ 1 };
    //const ushort e_pattern = l_tables.ePattern(i_case); // l_tables.e_pattern[i_case];
    for (int e = 0; e < 12; e++) {
        v_gindex[e] = -1;
        if (flag & e_) {
            // compute edge inersection
            // compute local coordinate along edge
            const int v0 = (l_edges_[e] & 0xF);
            const int v1 = (l_edges_[e] >> 4) & 0xF;
            const float l = (i0 - u[v0]) / (u[v1] - u[v0]);
            float4 vp = make_float4(v[v0].x + l*(v[v1].x - v[v0].x), v[v0].y + l*(v[v1].y - v[v0].y), v[v0].z + l*(v[v1].z - v[v0].z), 1.f);
            float4 np = make_float4(n[v0].x + l*(n[v1].x - n[v0].x), n[v0].y + l*(n[v1].y - n[v0].y), n[v0].z + l*(n[v1].z - n[v0].z), 0.f);
            const float length = std::sqrt(np.x * np.x + np.y * np.y + np.z * np.z);
            np.x = np.x / length;
            np.y = np.y / length;
            np.z = np.z / length;

            // get unique vertex index
            // compute vertex global index
            const int ix = i_index + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
            const int iy = j_index + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
            const int iz = k_index + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
            const int off_val = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
            v_gindex[e] = insert_vertex(int(9 * ugrid.gl_index(ix, iy, iz) + off_val), ht_, v_, vp, np);
        }
        flag <<= 1;
    }

    // compute triangles
    //const unsigned char* t_ambig = l_tables.t_ambig;
    unsigned long long tl_ = l_tables.t_[i_case];
    for (int t = 0; t < 16; t += 3) {
        //const int t_index = i_case * 16 + t;
        //if (t_pattern[i_case * 16 + t] == -1) {
        if (((tl_ >> (4 * t)) & 0xFull) == 0xF) {
            // there are no more triangles
            break;
        }
        // save tirangle
        const int v0 = (int)((tl_ >> (4 * t)) & 0xFull);
        const int v1 = (int)((tl_ >> (4 * (t + 1))) & 0xFull);
        const int v2 = (int)((tl_ >> (4 * (t + 2))) & 0xFull);
        // create three halfedges
        addHalfedges(he_, het_, v_gindex[v0], v_gindex[v1], v_gindex[v2]);
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute ambiguous cells
__global__ void he_tSlice(const float i0, cudaTextureObject_t v_data, UGrid ugrid, MC_lookup l_tables, const int nr_cells, AmbiguousCells ac_, VertexHashTable ht_,Vertices v_, HalfedgeHashTable het_, Halfedges he_)
{
    // get cell id
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nr_cells <= tid)
        return;

    // compute grid indices from global index
    const int gl_index = ac_.a_cells[tid];
    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);

    // construct 8 cell vertices
    float3 p[8];
    cell_vertices(p, i_index, j_index, k_index, ugrid);

    // scalar values at vertices
    float F[8];
    F[0] = tex3D<float>(v_data, i_index, j_index, k_index);
    F[1] = tex3D<float>(v_data, i_index + 1, j_index, k_index);
    F[2] = tex3D<float>(v_data, i_index, j_index + 1, k_index);
    F[3] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index);
    F[4] = tex3D<float>(v_data, i_index, j_index, k_index + 1);
    F[5] = tex3D<float>(v_data, i_index + 1, j_index, k_index + 1);
    F[6] = tex3D<float>(v_data, i_index, j_index + 1, k_index + 1);
    F[7] = tex3D<float>(v_data, i_index + 1, j_index + 1, k_index + 1);

    // compute normals at vertices
    float3 n[8];
    gradient(n, v_data, ugrid, i_index, j_index, k_index);

    // compute case
    uint i_case{ 0 };
    i_case = i_case + ((uint)(F[0] >= i0));
    i_case = i_case + ((uint)(F[1] >= i0)) * 2;
    i_case = i_case + ((uint)(F[2] >= i0)) * 4;
    i_case = i_case + ((uint)(F[3] >= i0)) * 8;
    i_case = i_case + ((uint)(F[4] >= i0)) * 16;
    i_case = i_case + ((uint)(F[5] >= i0)) * 32;
    i_case = i_case + ((uint)(F[6] >= i0)) * 64;
    i_case = i_case + ((uint)(F[7] >= i0)) * 128;

    // Compute intersection with edges
    const unsigned long long gei_pattern_ = 670526590282893600ull;
    const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };

    // compute intersection with cell edges
    float ecoord[12]{};
    int v_gindex[12]{};
    ushort flag{ 1 };
    ushort e_ = l_tables.e_[i_case];
    //ushort e_pattern = l_tables.ePattern(i_case); // l_tables.e_pattern[i_case];
    for (int e = 0; e < 12; e++) {
        v_gindex[e] = -1;
        //ecoord[e] = 0.f;
        if (flag & e_) {
            // compute edge inersection
            // compute local coordinate along edge
            const int v0 = (l_edges_[e] & 0xF);
            const int v1 = (l_edges_[e] >> 4) & 0xF;
            const float l = (i0 - F[v0]) / (F[v1] - F[v0]);
            float4 vp = make_float4(p[v0].x + l*(p[v1].x - p[v0].x), p[v0].y + l*(p[v1].y - p[v0].y), p[v0].z + l*(p[v1].z - p[v0].z), 1.f);
            float4 np = make_float4(n[v0].x + l*(n[v1].x - n[v0].x), n[v0].y + l*(n[v1].y - n[v0].y), n[v0].z + l*(n[v1].z - n[v0].z), 0.f);
            const float length = std::sqrt(np.x * np.x + np.y * np.y + np.z * np.z);
            np.x = np.x / length;
            np.y = np.y / length;
            np.z = np.z / length;

            // get unique vertex index
            // compute vertex global index
            const int ix = i_index + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
            const int iy = j_index + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
            const int iz = k_index + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
            const int off_val = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
            v_gindex[e] = insert_vertex(int(9 * ugrid.gl_index(ix, iy, iz) + off_val), ht_, v_, vp, np);

            // remember local coordinate along edge
            ecoord[e] = l;
        }
        flag <<= 1;
    }

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
        return (f0*f3 - f1*f2) / (f0 + f3 - f1 - f2);
    };
    uchar f_flag{ 0 };
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
            if (val > i0) {
                set_segm(e3, e0, segm_);
                set_segm(e1, e2, segm_);
            }
            else if (val < i0) {
                set_segm(e1, e0, segm_);
                set_segm(e3, e2, segm_);
            }
            else {
                // set flag for this face
                f_flag |= (1 << f);
                float ec0 = ecoord[e0];
                float ec1 = ecoord[e1];
                float ec2 = ecoord[e2];
                float ec3 = ecoord[e3];
                if ((0x218 >> (f * 2)) & BIT_1) {
                    ec0 = 1 - ec0;
                    ec2 = 1 - ec2;
                }
                if ((0x218 >> (f * 2)) & BIT_2) {
                    ec1 = 1 - ec1;
                    ec3 = 1 - ec3;
                }
                if (ec1 < ec3 && ec0 > ec2) {
                    set_segm(e1, e0, segm_);
                    set_segm(e3, e2, segm_);
                }
                else if (ec1 > ec3 && ec0 < ec2) {
                    set_segm(e3, e0, segm_);
                    set_segm(e1, e2, segm_);
                }
                else {
                    return;
                }
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
            if (val > i0) {
                set_segm(e0, e1, segm_);
                set_segm(e2, e3, segm_);
            }
            else if (val < i0) {
                set_segm(e0, e3, segm_);
                set_segm(e2, e1, segm_);
            }
            else {
                f_flag = (1 << f);
                // singular case val == i0, there are no asymptotes
                // check if there is a reasonable triangulation of the face
                float ec0 = ecoord[e0];
                float ec1 = ecoord[e1];
                float ec2 = ecoord[e2];
                float ec3 = ecoord[e3];
                if ((0x218 >> (f * 2)) & BIT_1) {
                    ec0 = 1 - ec0;
                    ec2 = 1 - ec2;
                }
                if ((0x218 >> (f * 2)) & BIT_2) {
                    ec1 = 1 - ec1;
                    ec3 = 1 - ec3;
                }
                if (ec1 < ec3 && ec0 > ec2) {
                    set_segm(e0, e1, segm_);
                    set_segm(e2, e3, segm_);
                }
                else if (ec1 > ec3 && ec0 < ec2) {
                    set_segm(e0, e3, segm_);
                    set_segm(e2, e1, segm_);
                }
                else {
                    return;
                }
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
    float ui[2]{};
    float vi[2]{};
    float wi[2]{};
    unsigned char q_sol{ 0 };
    const float a = (F[0] - F[1])*(-F[6] + F[7] + F[4] - F[5]) - (F[4] - F[5])*(-F[2] + F[3] + F[0] - F[1]);
    const float b = (i0 - F[0])*(-F[6] + F[7] + F[4] - F[5]) + (F[0] - F[1])*(F[6] - F[4]) - (i0 - F[4])*(-F[2] + F[3] + F[0] - F[1]) - (F[4] - F[5])*(F[2] - F[0]);
    const float c = (i0 - F[0])*(F[6] - F[4]) - (i0 - F[4])*(F[2] - F[0]);;
    float d = b*b - 4 * a*c;
    if (d > 0) {
        d = std::sqrt(d);
        // compute u-coord of solutions
        ui[0] = (-b - d) / (2 * a);
        ui[1] = (-b + d) / (2 * a);
        // compute v-coord of solutions
        float g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
        float g2 = F[2] * (1 - ui[0]) + F[3] * ui[0];
        vi[0] = (i0 - g1) / (g2 - g1);
        if (isnan(vi[0]) || isinf(vi[0])) {
            vi[0] = -1.f;
        }
        g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
        g2 = F[2] * (1 - ui[1]) + F[3] * ui[1];
        vi[1] = (i0 - g1) / (g2 - g1);
        if (isnan(vi[1]) || isinf(vi[1])) {
            vi[1] = -1.f;
        }
        // compute w-coordinates of solutions
        g1 = F[0] * (1 - ui[0]) + F[1] * ui[0];
        g2 = F[4] * (1 - ui[0]) + F[5] * ui[0];
        wi[0] = (i0 - g1) / (g2 - g1);
        if (isnan(wi[0]) || isinf(wi[0])) {
            wi[0] = -1.f;
        }
        g1 = F[0] * (1 - ui[1]) + F[1] * ui[1];
        g2 = F[4] * (1 - ui[1]) + F[5] * ui[1];
        wi[1] = (i0 - g1) / (g2 - g1);
        if (isnan(wi[1]) || isinf(wi[1])) {
            wi[1] = -1.f;
        }

        // correct values for roots of quadratic equations
        // in case the asymptotic decider has failed
        if (f_flag & BIT_1) { // face 1, w = 0;
            if (wi[0] < wi[1]) wi[0] = 0;
            else wi[1] = 0;
        }
        if (f_flag & BIT_2) { // face 2, w = 1
            if (wi[0] > wi[1]) wi[1] = 1;
            else wi[1] = 1;
        }
        if (f_flag & BIT_3) { // face 3, v = 0
            if (vi[0] < vi[1]) vi[0] = 0;
            else vi[1] = 0;
        }
        if (f_flag & BIT_4) { // face 4, v = 1
            if (vi[0] > vi[1]) vi[0] = 1;
            else vi[1] = 1;
        }
        if (f_flag & BIT_5) { // face 5, u = 0
            if (ui[0] < ui[1]) ui[0] = 0;
            else ui[1] = 0;
        }
        if (f_flag & BIT_6) { // face 6, u = 1
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

    // compute the number of solutions to the quadratic equation for a given face
    auto nrQSolFace = [](const uint f, const unsigned char n) {
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
    // Triangles are stored in global memory starting at offset
    // count nr. of inner vertices to compute right global index
    // first inner vertex has index cell_global_index + 3;
    int v_count{ 3 };
    if (numberOfSetBits<unsigned char>(q_sol) == 6) {
        // there are at most three contours
        // Possible cases:
        //  1) a single contour with 12 vertices
        //  2) two contours which build a tunnel
        //  3) three contours, one has only 3 vertices and does not belong to the tunnel

        // construct the six vertices of the inner hexagon
        float3 hvt[6];
        hvt[0].x = ui[0]; hvt[0].y = vi[0]; hvt[0].z = wi[0];
        hvt[1].x = ui[0]; hvt[1].y = vi[0]; hvt[1].z = wi[1];
        hvt[2].x = ui[1]; hvt[2].y = vi[0]; hvt[2].z = wi[1];
        hvt[3].x = ui[1]; hvt[3].y = vi[1]; hvt[3].z = wi[1];
        hvt[4].x = ui[1]; hvt[4].y = vi[1]; hvt[4].z = wi[0];
        hvt[5].x = ui[0]; hvt[5].y = vi[1]; hvt[5].z = wi[0];

        // construct vertices at intersections with the edges
        auto e_vert = [&ecoord](const int e, const int i) {
            const unsigned int l_coord[3]{ 1324855, 5299420, 16733440 };
            unsigned char flag = (l_coord[i] >> (2 * e)) & 3;
            if (flag == 3)
                return ecoord[e];
            else
                return (float)(flag);

        };

        // if there are three contours, then there is a tunnel and one
        // of the contours is not part of it.
        unsigned char _not_tunnel = 0xF;
        if (cnt_ == 3) {
            // loop over the contorus
            // triangulate the contour which is not part of
            // the tunnel
            const float uc_min = (ui[0] < ui[1]) ? ui[0] : ui[1];
            const float uc_max = (ui[0] < ui[1]) ? ui[1] : ui[0];
            for (int t = 0; t < (int)cnt_; t++) {
                if (get_cnt_size(t, c_) == 3) {
                    float umin = 2;
                    float umax = -2;
                    uint e0 = get_c(t, 0, c_);
                    uint e1 = get_c(t, 1, c_);
                    uint e2 = get_c(t, 2, c_);
                    const float u_e0 = e_vert(e0, 0);
                    const float u_e1 = e_vert(e1, 0);
                    const float u_e2 = e_vert(e2, 0);
                    umin = (u_e0 < umin) ? u_e0 : umin;
                    umin = (u_e1 < umin) ? u_e1 : umin;
                    umin = (u_e2 < umin) ? u_e2 : umin;
                    umax = (u_e0 > umax) ? u_e0 : umax;
                    umax = (u_e1 > umax) ? u_e1 : umax;
                    umax = (u_e2 > umax) ? u_e1 : umax;
                    if (uc_min > umax || uc_max < umin) {
                        // this contour is not part of the tunnel
                        _not_tunnel = t;
                        // save triangle in global memory
                        addHalfedges(he_, het_, v_gindex[e0], v_gindex[e1], v_gindex[e2]);
                        //const int a_ = atomicAdd(he_cnt, 3);
                        //addHalfedges(nr_he, he_e, he_table, he_ids, a_, v_gindex[e0], v_gindex[e1], v_gindex[e2]);
                        //t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[e0], v_gindex[e1], v_gindex[e2], 0);
                    }
                }
            }
        }

        // compute vertices of inner hexagon, save new vertices in list and compute and keep
        // global vertice index to build triangle connectivity later on.
        int tg_idx[6];
        for (int i = 0; i < 6; i++) {
            float4  hp;
            float4 hn;
            // local coordinates for trilinear interpolation
            const float u = hvt[i].x; const float v = hvt[i].y; const float w = hvt[i].z;
            hp.x = (1 - w)*((1 - v)*(p[0].x + u*(p[1].x - p[0].x)) + v*(p[2].x + u*(p[3].x - p[2].x))) + w*((1 - v)*(p[4].x + u*(p[5].x - p[4].x)) + v*(p[6].x + u*(p[7].x - p[6].x)));
            hp.y = (1 - w)*((1 - v)*(p[0].y + u*(p[1].y - p[0].y)) + v*(p[2].y + u*(p[3].y - p[2].y))) + w*((1 - v)*(p[4].y + u*(p[5].y - p[4].y)) + v*(p[6].y + u*(p[7].y - p[6].y)));
            hp.z = (1 - w)*((1 - v)*(p[0].z + u*(p[1].z - p[0].z)) + v*(p[2].z + u*(p[3].z - p[2].z))) + w*((1 - v)*(p[4].z + u*(p[5].z - p[4].z)) + v*(p[6].z + u*(p[7].z - p[6].z)));
            hn.x = (1 - w)*((1 - v)*(n[0].x + u*(n[1].x - n[0].x)) + v*(n[2].x + u*(n[3].x - n[2].x))) + w*((1 - v)*(n[4].x + u*(n[5].x - n[4].x)) + v*(n[6].x + u*(n[7].x - n[6].x)));
            hn.y = (1 - w)*((1 - v)*(n[0].y + u*(n[1].y - n[0].y)) + v*(n[2].y + u*(n[3].y - n[2].y))) + w*((1 - v)*(n[4].y + u*(n[5].y - n[4].y)) + v*(n[6].y + u*(n[7].y - n[6].y)));
            hn.z = (1 - w)*((1 - v)*(n[0].z + u*(n[1].z - n[0].z)) + v*(n[2].z + u*(n[3].z - n[2].z))) + w*((1 - v)*(n[4].z + u*(n[5].z - n[4].z)) + v*(n[6].z + u*(n[7].z - n[6].z)));
            // normalize normal
            const float factor = std::sqrt(hn.x * hn.x + hn.y * hn.y + hn.z * hn.z);
            hn.x = hn.x / factor;
            hn.y = hn.y / factor;
            hn.z = hn.z / factor;
            // the fourth coord.
            hp.w = 1.f;
            hn.w = 0.f;
            // this vertices are inner vertices
            tg_idx[i] = insert_vertex(int(9 * gl_index + v_count),ht_,v_,hp,hn);
            
            //int v_addr{ -1 };
            //if (insert_vertex_key(tg_idx[i], ht_, v_addr)) {
            //    //const int pos = atomicAdd(v_.t_size, 1);
            //    //v_.vertices[pos] = hp;
            //    //v_.normals[pos] = hn;
            //    //v_.v[pos].v = hp;
            //    //v_.v[pos].n = hn;
            //    //ht_.addr[v_addr] = pos;
            //    ht_.addr[v_addr] = v_.set(hp, hn);
            //}
            // update nr. of vertices
            v_count++;
        }


        // triangulate contours with inner hexagon
        unsigned char tcon_[12];
        for (int i = 0; i < (int)cnt_; i++) {
            if (_not_tunnel != i) { // contour belongs to tunnel
                const int cnt_sz = (int)get_cnt_size(i, c_);
                for (int r = 0; r < cnt_sz; r++) {
                    int index = -1;
                    double dist = 1000.;
                    uint ci = get_c(i, r, c_);
                    const double u_edge = e_vert(ci, 0);
                    const double v_edge = e_vert(ci, 1);
                    const double w_edge = e_vert(ci, 2);
                    for (int s = 0; s < 6; s++) {
                        const double uval = u_edge - hvt[s].x;
                        const double vval = v_edge - hvt[s].y;
                        const double wval = w_edge - hvt[s].z;
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
                        addHalfedges(he_, het_, v_gindex[tid1], v_gindex[tid2], tg_idx[cid1]);
                        //const int a_ = atomicAdd(he_cnt, 3);
                        //addHalfedges(nr_he, he_e, he_table, he_ids, a_, v_gindex[tid1], v_gindex[tid2], tg_idx[cid1]);
                        //t_.triangles[atomicAdd(t_.t_size, 1)] = make_int4(v_gindex[tid1], v_gindex[tid2], tg_idx[cid1], 0);
                    }
                    break;
                    case 1:
                    {
                        // measure diagonals
                        // triangulate along shortest diagonal
                        float u_edge = e_vert(tid1, 0);
                        float v_edge = e_vert(tid1, 1);
                        float w_edge = e_vert(tid1, 2);
                        const float l1 = (u_edge - hvt[cid2].x)*(u_edge - hvt[cid2].x) + (v_edge - hvt[cid2].y)*(v_edge - hvt[cid2].y) + (w_edge - hvt[cid2].z)*(w_edge - hvt[cid2].z);
                        u_edge = e_vert(tid2, 0);
                        v_edge = e_vert(tid2, 1);
                        w_edge = e_vert(tid2, 2);
                        const double l2 = (u_edge - hvt[cid1].x)*(u_edge - hvt[cid1].x) + (v_edge - hvt[cid1].y)*(v_edge - hvt[cid1].y) + (w_edge - hvt[cid1].z)*(w_edge - hvt[cid1].z);
                        if (l1 < l2) {
                            addHalfedges(he_, het_, v_gindex[tid1], v_gindex[tid2], tg_idx[cid2]);
                            addHalfedges(he_, het_, v_gindex[tid1], tg_idx[cid2], tg_idx[cid1]);
                        }
                        else {
                            addHalfedges(he_, het_, v_gindex[tid1], v_gindex[tid2], tg_idx[cid1]);
                            addHalfedges(he_, het_, v_gindex[tid2], tg_idx[cid2], tg_idx[cid1]);
                        }
                    }
                    break;
                    case 2:
                    {
                        const int cidm = midpointRingIntModulo(cid1, cid2);
                        addHalfedges(he_, het_, v_gindex[tid1], v_gindex[tid2], tg_idx[cidm]);
                        addHalfedges(he_, het_, v_gindex[tid1], tg_idx[cidm], tg_idx[cid1]);
                        addHalfedges(he_, het_, v_gindex[tid2], tg_idx[cid2], tg_idx[cidm]);
                    }
                    break;
                    } // switch
                } // for loop over the vertices of the contour
            } // if (_not_tunnel)
        } // for loop over contours
        if (cnt_ == 1) {
            // there is a single contour
            // triangulate and close inner hexagon
            addHalfedges(he_, het_, tg_idx[0], tg_idx[2], tg_idx[1]);
            addHalfedges(he_, het_, tg_idx[2], tg_idx[4], tg_idx[3]);
            addHalfedges(he_, het_, tg_idx[0], tg_idx[5], tg_idx[4]);
            addHalfedges(he_, het_, tg_idx[0], tg_idx[4], tg_idx[2]);
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
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)]);
                }
                break;
                case 4:
                {
                    
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)]);
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 2, c_)], v_gindex[get_c(i, 3, c_)]);
                }
                break;
                case 5:
                {
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)]);
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 2, c_)], v_gindex[get_c(i, 3, c_)]);
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 3, c_)], v_gindex[get_c(i, 4, c_)]);
                }
                break;
                case 6:
                {
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 3, c_)]);
                    addHalfedges(he_, het_, v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)], v_gindex[get_c(i, 3, c_)]);
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 3, c_)], v_gindex[get_c(i, 4, c_)]);
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 4, c_)], v_gindex[get_c(i, 5, c_)]);
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
            unsigned char fs[3][2]{ { (unsigned char)(q_sol & 1), (unsigned char)((q_sol >> 1) & 1) },{ (unsigned char)((q_sol >> 2) & 1), (unsigned char)((q_sol >> 3) & 1) },{ (unsigned char)((q_sol >> 4) & 1), (unsigned char)((q_sol >> 5) & 1) } };

            const unsigned char fc1 = fs[0][0] * fs[1][0] + fs[0][1] * fs[1][1];
            const unsigned char fc2 = fs[0][0] * fs[2][0] + fs[0][1] * fs[2][1];
            const unsigned char fc3 = fs[1][0] * fs[2][1] + fs[1][1] * fs[2][0];
            const unsigned char c_faces = fc1 + fc2 + fc3;
            float ucoord{};
            float vcoord{};
            float wcoord{};
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
            float4 ip;
            float4 in;
            ip.x = (1 - wcoord)*((1 - vcoord)*(p[0].x + ucoord*(p[1].x - p[0].x)) + vcoord*(p[2].x + ucoord*(p[3].x - p[2].x))) + wcoord*((1 - vcoord)*(p[4].x + ucoord*(p[5].x - p[4].x)) + vcoord*(p[6].x + ucoord*(p[7].x - p[6].x)));
            ip.y = (1 - wcoord)*((1 - vcoord)*(p[0].y + ucoord*(p[1].y - p[0].y)) + vcoord*(p[2].y + ucoord*(p[3].y - p[2].y))) + wcoord*((1 - vcoord)*(p[4].y + ucoord*(p[5].y - p[4].y)) + vcoord*(p[6].y + ucoord*(p[7].y - p[6].y)));
            ip.z = (1 - wcoord)*((1 - vcoord)*(p[0].z + ucoord*(p[1].z - p[0].z)) + vcoord*(p[2].z + ucoord*(p[3].z - p[2].z))) + wcoord*((1 - vcoord)*(p[4].z + ucoord*(p[5].z - p[4].z)) + vcoord*(p[6].z + ucoord*(p[7].z - p[6].z)));
            in.x = (1 - wcoord)*((1 - vcoord)*(n[0].x + ucoord*(n[1].x - n[0].x)) + vcoord*(n[2].x + ucoord*(n[3].x - n[2].x))) + wcoord*((1 - vcoord)*(n[4].x + ucoord*(n[5].x - n[4].x)) + vcoord*(n[6].x + ucoord*(n[7].x - n[6].x)));
            in.y = (1 - wcoord)*((1 - vcoord)*(n[0].y + ucoord*(n[1].y - n[0].y)) + vcoord*(n[2].y + ucoord*(n[3].y - n[2].y))) + wcoord*((1 - vcoord)*(n[4].y + ucoord*(n[5].y - n[4].y)) + vcoord*(n[6].y + ucoord*(n[7].y - n[6].y)));
            in.z = (1 - wcoord)*((1 - vcoord)*(n[0].z + ucoord*(n[1].z - n[0].z)) + vcoord*(n[2].z + ucoord*(n[3].z - n[2].z))) + wcoord*((1 - vcoord)*(n[4].z + ucoord*(n[5].z - n[4].z)) + vcoord*(n[6].z + ucoord*(n[7].z - n[6].z)));
            // normalize normal
            const float factor = std::sqrt(in.x * in.x + in.y * in.y + in.z * in.z);
            in.x = in.x / factor;
            in.y = in.y / factor;
            in.z = in.z / factor;
            // the fourth coordinate
            ip.w = 1.f;
            in.w = 0.f;
            // global index
            int gidx = int(9 * gl_index + v_count);
            // this point is only used if contours with more than three vertices
            // are present
            for (int i = 0; i < (int)cnt_; i++) {
                if (get_cnt_size(i, c_) > 3) {
                    gidx = insert_vertex(gidx, ht_, v_, ip, in);
                }
            }
            //bool pt_used{ false };

            // loop over the contorus
            for (int i = 0; i < (int)cnt_; i++) {
                const unsigned char cnt_sz = (unsigned char)get_cnt_size(i, c_);
                if (cnt_sz == 3) {
                    addHalfedges(he_, het_, v_gindex[get_c(i, 0, c_)], v_gindex[get_c(i, 1, c_)], v_gindex[get_c(i, 2, c_)]);
                }
                else {
                    //pt_used = true;
                    for (int t = 0; t < cnt_sz; t++) {
                        // add triangle to list
                        addHalfedges(he_, het_, v_gindex[get_c(i, t, c_)], v_gindex[get_c(i, (t + 1) % cnt_sz, c_)], gidx);
                    }
                }
            }
        } // else - there are saddle points
    }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      HOST CODE
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// read a volume data from file
void
p_mc::MarchingCubes::readDataFromFile(const std::string& i_file, std::array<int,3>& dim, std::array<float,3>& origin, std::array<float,3>& spacing, std::vector<float>& v_data)
{
	std::FILE* f{ nullptr };
	errno_t status = fopen_s(&f, i_file.c_str(), "rb");
	if (status != 0) {
		std::cerr << "ERROR: can't open file " << i_file.c_str() << std::endl;
		exit(1);
	}
	short x_size;
	short y_size;
	short z_size;
	std::fread(&x_size, sizeof(unsigned short), 1, f);
	std::fread(&y_size, sizeof(unsigned short), 1, f);
	std::fread(&z_size, sizeof(unsigned short), 1, f);

	float dx;
	float dy;
	float dz;
	std::fread(&dx, sizeof(float), 1, f);
	std::fread(&dy, sizeof(float), 1, f);
	std::fread(&dz, sizeof(float), 1, f);

	int v_size = x_size * y_size * z_size;
	ushort* v_buff = new ushort[v_size];
	std::fread(&v_buff[0], sizeof(unsigned short), v_size, f);

	// fill into host vector
	v_data.resize(v_size);
	for (int i = 0; i < v_size; i++) {
		v_data[i] = float(v_buff[i]);
	}
	std::fclose(f);
	delete[] v_buff;

	// set uniform grid data
    dim[0] = x_size;
    dim[1] = y_size;
    dim[2] = z_size;
    origin[0] = 0.f;
    origin[1] = 0.f;
    origin[2] = 0.f;
    spacing[0] = dx;
    spacing[1] = dy;
    spacing[2] = dz;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main function to process a volume data set
void
p_mc::MarchingCubes::mc_halfedge(const float i0,const std::string& i_file, int& nr_v, float** vertices, float** normals, int& nr_t, int4** h_hee, int** h_hev,int** h_hef)
{
	Halfedges he_;
    HalfedgeHashTable het_;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << " ... read data from file\n";
	std::vector<float> h_data;
    UGrid ugrid;
    std::array<int, 3> dims;
    std::array<float, 3> origin;
    std::array<float, 3> spacing;
	readDataFromFile(i_file, dims, origin, spacing, h_data);
    ugrid.size(dims[0], dims[1], dims[2]);
    ugrid.dx = spacing[0];
    ugrid.dy = spacing[1];
    ugrid.dz = spacing[2];
    ugrid.x0 = origin[0];
    ugrid.y0 = origin[1];
    ugrid.z0 = origin[2];

	// measure processing time
	CTimer ctimer1;
	CTimer ctimer2;

	// allocate 3D texture
	std::cout << " ... allocate 3D texture\n";
	const size_t x_size = (size_t)ugrid.idim;
	const size_t y_size = (size_t)ugrid.jdim;
	const size_t z_size = (size_t)ugrid.kdim;
	const size_t t_size = x_size * y_size * z_size;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// create texture buffer for 3D data
	// copy data to buffer
	cudaArray* d_data;
	cudaExtent extent = make_cudaExtent(x_size, y_size, z_size);
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
	cudaMalloc3DArray(&d_data, &desc, extent);
	//cudaCheckError();

	cudaMemcpy3DParms params{ 0 };
	params.srcPtr = make_cudaPitchedPtr(&(h_data[0]), x_size * sizeof(float), x_size, y_size);
	params.dstArray = d_data;
	params.extent = extent;
	params.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&params);
	//cudaCheckError();

	// create Texture object
	// Texture description
	cudaTextureDesc texDesc{};
    memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.readMode = cudaReadModeElementType;
	texDesc.filterMode = cudaFilterModePoint;
	// Texture resource description
	cudaResourceDesc resDesc{};
    memset(&resDesc, 0, sizeof(cudaResourceDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_data;
	// create Texture object
	cudaCreateTextureObject(&m_volume, &resDesc, &texDesc, nullptr);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// copy lookup tables
    //ctimer2.start();
	MC_lookup l_tables;
    initMC_lookup(l_tables, e_pattern, t_pattern, t_ambig);
    //ctimer2.stop();
    //ctimer2.print(std::string(" ... lookup tables"));

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // globa processing time
    std::cout << " ... compute isosurface\n";
    ctimer1.start();
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Atomic counter
	// Alloc memory for counting nr. of vertices
    //ctimer2.start();
    CellsIds cells;
    initCells(cells,(int)(t_size / 3));
    AmbiguousCells acells;
    initACells(acells,(int)(t_size / 4));
    int* d_vcount{ nullptr };
    int* d_acount{ nullptr };
    cudaMalloc(&d_vcount, sizeof(int));
    cudaMemset(d_vcount, 0, sizeof(int));
    cudaMalloc(&d_acount, sizeof(int));
    cudaMemset(d_acount, 0, sizeof(int));
	uint b_size = 512;
	uint g_size{ ((uint)t_size + b_size - 1) / b_size };
	mc_count << < g_size, b_size >> >(cells, acells, d_vcount, d_acount, i0, m_volume, ugrid, l_tables);
    //cudaCheckError();


	// count
    // each vertex is counted four times, except those at the boundaries
    // inner vertices are overestimated.
	//nr_v = (int)(warpReduce<int>(d_vcount, (int)t_size)) / 2;
	// read array size
	//int nr_cells{ 0 };
	//cudaMemcpy(&nr_cells, d_cpos, sizeof(int), cudaMemcpyDeviceToHost);
    //cudaCheckError();
    int nr1{ 0 };
    int nr2{ 0 };
    cudaMemcpy(&nr1, d_vcount, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(&nr2, d_acount, sizeof(int), cudaMemcpyDeviceToHost);
    nr_v = (nr1 + nr2) / 2;
    //ctimer2.stop();
    //ctimer2.print(std::string(" ... count vertices"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute triangles
    // 1. alloc memory for hash table
    // 2. alloc memory for vertices
    // 3. alloc memory for triangles
    // 4. alloc memory for computing cell ids of ambiguous cases
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. alloc and init hash table
    //ctimer2.start();
    VertexHashTable ht_;
    initVertexHashTable(ht_, nr_v); // appro. two times the number of vertices
    b_size = 512;
    g_size = (ht_.t_size + b_size - 1) / b_size;
    init_hash_table << < g_size, b_size >> >(ht_);
    //cudaCheckError();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. alloc and init vertices
    Vertices v_;
    initVertices(v_, nr_v);
    //cudaCheckError();
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. alloc and init triangles
    Triangles t_;
    initTriangles(t_, 2 * nr_v);
    //cudaCheckError();
    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute iso-surface
    const int nr_cells = size<CellsIds>(cells);
    const int nr_acells = size<AmbiguousCells>(acells);
	b_size = MC_BLOCKSIZE;
	g_size = (nr_cells + b_size - 1) / b_size;
    mc_slice << < g_size, b_size >> > (i0, m_volume, ugrid, l_tables, nr_cells, cells, ht_, v_, t_);
    //cudaCheckError();
    //ctimer2.stop();
    //ctimer2.print(std::string(" ... mc_slice()"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute ambiguous cases
    //ctimer2.start();
	b_size = AMB_BLOCKSIZE;
	g_size = (nr_acells + b_size - 1) / b_size;
	t_slice << < g_size, b_size >> >(i0, m_volume, ugrid, l_tables, nr_acells, acells, ht_, v_, t_);
	//cudaCheckError();
    //ctimer2.stop();
    //ctimer2.print(std::string(" ... ambiguous cases"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// don't need volume data any more
	//cudaFreeArray(d_data);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute shared vertex list for triangle mesh
    // indices of triangles have to be mapped to vertex index in vertex array
    // get number of vertices
    //ctimer2.start();
    nr_v = size<Vertices>(v_);
    // get number of triangles
    nr_t = size<Triangles>(t_);

	// map triangles indices
	b_size = 512;
	g_size = (3 * nr_t + b_size - 1) / b_size;
	map_triangles_fast <<< g_size, b_size >>>(nr_t,ht_, t_);
	//cudaCheckError();
    //ctimer2.stop();
    //ctimer2.print(std::string(" ... map triangles"));

    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// create halfedge data structure for triangle mesh
	// 1. there are three edges for each vertex
	//    an edge points to origin vertex
	//    an edge points to the face
	//    an edge points to next edge
	// 2. each vertex has points to a halfedge
	// 3. each triangle points to a halfedge
	// Data structure:
	// halfedge int4: 
	//    he.x = origin vertex
	//    he.y = face
	//    he.z = next
	//   he.w = tween
	// vertex int: point to one halfedge
	// face int: point to one halfedge of this face
	//ctimer2.start();
	const int nr_he = 3 * nr_t;
	initHalfedges(he_, nr_he, nr_v, nr_t);
	initHalfedgeHashTable(het_, 2 * nr_he); // the hash table is 1.5 the total nr. of halfedges
	//cudaCheckError();
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... allocate memory for halfedge data structure"));
	// for each triangle create three halfedges
	// compute unique halfedge twin ids and store in hash table
	//ctimer2.start();
	b_size = 256;
	g_size = (nr_he + b_size - 1) / b_size;
	create_halfedge << < g_size, b_size >> > (nr_he, t_, he_, het_);
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... create halfedges"));
	// connect each half edge with its twin
	// process each entrie in halfedge table
	//ctimer2.start();
	g_size = (het_.t_size + b_size - 1) / b_size;
	map_halfedge_twins_fast << < g_size, b_size >> > (he_, het_);
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... map halfedges"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute processing time
	ctimer1.stop();
    cudaDeviceSynchronize();
    ctimer1.print(std::string("tmc"));
    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// copy data back to host
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Shared vertex list
    float4* v_array = new float4[nr_v];
    float4* n_array = new float4[nr_v];
    int4* t_array = new int4[nr_t];

    cudaMemcpy(v_array, v_.vertices, nr_v * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(n_array, v_.normals,  nr_v * sizeof(float4), cudaMemcpyDeviceToHost);
    cudaMemcpy(t_array, t_.triangles, nr_t * sizeof(int4), cudaMemcpyDeviceToHost);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // halfedge data structure
    *h_hee = new int4[nr_he];
    *h_hev = new int[nr_v];
    *h_hef = new int[nr_t];
    cudaMemcpy(*h_hee, he_.he_e, nr_he * sizeof(int4), cudaMemcpyDeviceToHost);
    cudaMemcpy(*h_hev, he_.he_v, nr_v  * sizeof(int),  cudaMemcpyDeviceToHost);
    cudaMemcpy(*h_hef, he_.he_f, nr_t  * sizeof(int),  cudaMemcpyDeviceToHost);

    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of triangles " << nr_t << std::endl;
    std::cout << " ... total nr. of unambiguous cells " << nr1 << std::endl;
    std::cout << " ... total nr. of ambiguous cells " << nr2 << std::endl;

    *vertices = new float[3 * nr_v];
    *normals = new float[3 * nr_v];
    for (int id = 0; id < nr_v; id++) {
        // copy vertices
        (*vertices)[3 * id] = v_array[id].x;
        (*vertices)[3 * id + 1] = v_array[id].y;
        (*vertices)[3 * id + 2] = v_array[id].z;
        // copy normals
        (*normals)[3 * id] = -n_array[id].x;
        (*normals)[3 * id + 1] = -n_array[id].y;
        (*normals)[3 * id + 2] = -n_array[id].z;
    }

    std::cout << " ... done\n";

    // host memory
    // host memory
    delete[] v_array;
    delete[] n_array;
    delete[] t_array;
    //delete[] h_hee;
    //delete[] h_hev;
    //delete[] h_hef;
 

    // free common data
    // free memory
    cudaFreeArray(d_data);
    freeMC_lookup(l_tables);
    freeVertices(v_);
    freeTriangles(t_);
    freeHalfedges(he_);
    freeHalfedgeHashTable(het_);
    freeCells(cells);
    freeACells(acells);

    // arrays for vertex count
    cudaFree(d_acount);
    cudaFree(d_vcount);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      HALFEDGE RECONSTRUCTION
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main function to process a volume data set
void
p_mc::MarchingCubes::mc_sharedvertex(const float i0, const std::string& i_file, int& nr_v, float** vertices, float** normals, int& nr_t, int** triangles)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout << " ... read data from file\n";
	std::vector<float> h_data;
	UGrid ugrid;
	std::array<int, 3> dims;
	std::array<float, 3> origin;
	std::array<float, 3> spacing;
	readDataFromFile(i_file, dims, origin, spacing, h_data);
	ugrid.size(dims[0], dims[1], dims[2]);
	ugrid.dx = spacing[0];
	ugrid.dy = spacing[1];
	ugrid.dz = spacing[2];
	ugrid.x0 = origin[0];
	ugrid.y0 = origin[1];
	ugrid.z0 = origin[2];

	// measure processing time
	CTimer ctimer1;
	CTimer ctimer2;

	// allocate 3D texture
	std::cout << " ... allocate 3D texture\n";
	const size_t x_size = (size_t)ugrid.idim;
	const size_t y_size = (size_t)ugrid.jdim;
	const size_t z_size = (size_t)ugrid.kdim;
	const size_t t_size = x_size * y_size * z_size;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// create texture buffer for 3D data
	// copy data to buffer
	cudaArray* d_data;
	cudaExtent extent = make_cudaExtent(x_size, y_size, z_size);
	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
	cudaMalloc3DArray(&d_data, &desc, extent);
	//cudaCheckError();

	cudaMemcpy3DParms params{ 0 };
	params.srcPtr = make_cudaPitchedPtr(&(h_data[0]), x_size * sizeof(float), x_size, y_size);
	params.dstArray = d_data;
	params.extent = extent;
	params.kind = cudaMemcpyHostToDevice;
	cudaMemcpy3D(&params);
	//cudaCheckError();

	// create Texture object
	// Texture description
	cudaTextureDesc texDesc{};
	memset(&texDesc, 0, sizeof(cudaTextureDesc));
	texDesc.readMode = cudaReadModeElementType;
	texDesc.filterMode = cudaFilterModePoint;
	// Texture resource description
	cudaResourceDesc resDesc{};
	memset(&resDesc, 0, sizeof(cudaResourceDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_data;
	// create Texture object
	cudaCreateTextureObject(&m_volume, &resDesc, &texDesc, nullptr);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// copy lookup tables
	//ctimer2.start();
	MC_lookup l_tables;
	initMC_lookup(l_tables, e_pattern, t_pattern, t_ambig);
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... lookup tables"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// globa processing time
	std::cout << " ... compute isosurface\n";
	ctimer1.start();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Atomic counter
	// Alloc memory for counting nr. of vertices
	//ctimer2.start();
	CellsIds cells;
	initCells(cells, (int)(t_size / 3));
	AmbiguousCells acells;
	initACells(acells, (int)(t_size / 4));
	int* d_vcount{ nullptr };
	int* d_acount{ nullptr };
	cudaMalloc(&d_vcount, sizeof(int));
	cudaMemset(d_vcount, 0, sizeof(int));
	cudaMalloc(&d_acount, sizeof(int));
	cudaMemset(d_acount, 0, sizeof(int));
	uint b_size = 512;
	uint g_size{ ((uint)t_size + b_size - 1) / b_size };
	mc_count << < g_size, b_size >> >(cells, acells, d_vcount, d_acount, i0, m_volume, ugrid, l_tables);
	//cudaCheckError();


	// count
	// each vertex is counted four times, except those at the boundaries
	// inner vertices are overestimated.
	//nr_v = (int)(warpReduce<int>(d_vcount, (int)t_size)) / 2;
	// read array size
	//int nr_cells{ 0 };
	//cudaMemcpy(&nr_cells, d_cpos, sizeof(int), cudaMemcpyDeviceToHost);
	//cudaCheckError();
	int nr1{ 0 };
	int nr2{ 0 };
	cudaMemcpy(&nr1, d_vcount, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&nr2, d_acount, sizeof(int), cudaMemcpyDeviceToHost);
	nr_v = (nr1 + nr2) / 2;
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... count vertices"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute triangles
	// 1. alloc memory for hash table
	// 2. alloc memory for vertices
	// 3. alloc memory for triangles
	// 4. alloc memory for computing cell ids of ambiguous cases
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 1. alloc and init hash table
	//ctimer2.start();
	VertexHashTable ht_;
	initVertexHashTable(ht_, nr_v); // appro. two times the number of vertices
	b_size = 512;
	g_size = (ht_.t_size + b_size - 1) / b_size;
	init_hash_table << < g_size, b_size >> >(ht_);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 2. alloc and init vertices
	Vertices v_;
	initVertices(v_, nr_v);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 3. alloc and init triangles
	Triangles t_;
	initTriangles(t_, 2 * nr_v);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute iso-surface
	const int nr_cells = size<CellsIds>(cells);
	const int nr_acells = size<AmbiguousCells>(acells);
	b_size = MC_BLOCKSIZE;
	g_size = (nr_cells + b_size - 1) / b_size;
	mc_slice << < g_size, b_size >> > (i0, m_volume, ugrid, l_tables, nr_cells, cells, ht_, v_, t_);
	//cudaCheckError();
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... mc_slice()"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute ambiguous cases
	//ctimer2.start();
	b_size = AMB_BLOCKSIZE;
	g_size = (nr_acells + b_size - 1) / b_size;
	t_slice << < g_size, b_size >> >(i0, m_volume, ugrid, l_tables, nr_acells, acells, ht_, v_, t_);
	//cudaCheckError();
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... ambiguous cases"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// don't need volume data any more
	//cudaFreeArray(d_data);
	//cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute shared vertex list for triangle mesh
	// indices of triangles have to be mapped to vertex index in vertex array
	// get number of vertices
	//ctimer2.start();
	nr_v = size<Vertices>(v_);
	// get number of triangles
	nr_t = size<Triangles>(t_);

	// map triangles indices
	b_size = 512;
	g_size = (3 * nr_t + b_size - 1) / b_size;
	map_triangles_fast << < g_size, b_size >> >(nr_t, ht_, t_);
	//cudaCheckError();
	//ctimer2.stop();
	//ctimer2.print(std::string(" ... map triangles"));

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute processing time
	ctimer1.stop();
	cudaDeviceSynchronize();
	ctimer1.print(std::string("tmc"));

	
	// Shared vertex list
	float4* v_array = new float4[nr_v];
	float4* n_array = new float4[nr_v];
	//vAttr* vts = new vAttr[nr_v];
	int4* t_array = new int4[nr_t];

	cudaMemcpy(v_array, v_.vertices, nr_v * sizeof(float4), cudaMemcpyDeviceToHost);
	cudaMemcpy(n_array, v_.normals, nr_v * sizeof(float4), cudaMemcpyDeviceToHost);
	//cudaMemcpy(vts, v_.v, nr_v * sizeof(vAttr), cudaMemcpyDeviceToHost);
	cudaMemcpy(t_array, t_.triangles, nr_t * sizeof(int4), cudaMemcpyDeviceToHost);

	std::cout << " ... total nr. of vertices " << nr_v << std::endl;
	std::cout << " ... total nr. of triangles " << nr_t << std::endl;
	std::cout << " ... total nr. of unambiguous cells " << nr1 << std::endl;
	std::cout << " ... total nr. of ambiguous cells " << nr2 << std::endl;


	*vertices = new float[3 * nr_v];
	*normals = new float[3 * nr_v];
	*triangles = new int[3 * nr_t];
	for (int id = 0; id < nr_v; id++) {
		// copy vertices
		(*vertices)[3 * id] = v_array[id].x;
		(*vertices)[3 * id + 1] = v_array[id].y;
		(*vertices)[3 * id + 2] = v_array[id].z;
		(*vertices)[3 * id + 3] = 1.0f;
		//// copy normals
		(*normals)[3 * id] = -n_array[id].x;
		(*normals)[3 * id + 1] = -n_array[id].y;
		(*normals)[3 * id + 2] = -n_array[id].z;
		(*normals)[3 * id + 3] = 0.0f;
	}
	for (int id = 0; id < nr_t; id++) {
		(*triangles)[3 * id] = t_array[id].x;
		(*triangles)[3 * id + 1] = t_array[id].y;
		(*triangles)[3 * id + 2] = t_array[id].z;
	}
	std::cout << " ... done\n";
	// host memory
	delete[] v_array;
	delete[] n_array;
	//delete[] vts;
	delete[] t_array;

	// free common data
	// free memory
	cudaFreeArray(d_data);
	freeMC_lookup(l_tables);
	freeVertices(v_);
	freeTriangles(t_);
	freeCells(cells);
	freeACells(acells);

	// arrays for vertex count
	cudaFree(d_acount);
	cudaFree(d_vcount);
}
