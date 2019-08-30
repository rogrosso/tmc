#include "DualMarchingCubes.h"

// Lauch bounds
#define THREADS_PER_BLOCK 256

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
#define MC_BLOCKSIZE 128

// Mesh smoothing / remove elements with valence patterns 3333 and 3535
#define P_NEUT 0
#define P_3333 23
#define P_3535 29
#define P_NEIG 31

// MC ambiguous case
#define MC_AMBIGUOUS 105

// type aliases
// Introduce convenient aliases here
using uint = unsigned int;
using uchar = unsigned char;
using ushort = unsigned short;
using ullong = unsigned long long;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// error handling
#define cudaCheckError() { \
	cudaError_t e=cudaGetLastError(); \
	if(e!=cudaSuccess) { \
		std::ostringstream buf; \
		buf << "Cuda failure in file " << __FILE__ << ", line: " << __LINE__ << ", error: " << cudaGetErrorString(e) << "\n";  \
        std::cout << buf.str(); \
		exit(0); \
	} \
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
// Quality measure
struct QualityMeasure {
    /// buffer size
    int a_size{ 0 };
    /// element quality
    float* q_measure{ nullptr };
    std::shared_ptr<float> q_measure_;
    /// Constructors
    __host__ __device__ QualityMeasure() {  }
    /*__host__ __device__ QualityMeasure(const QualityMeasure& q)
    {
        this->a_size = q.a_size;
        this->q_measure = q.q_measure;
        this->q_measure_ = q.q_measure_;
    }*/
    __host__ QualityMeasure(const int sz) : a_size{ sz }
    {
        cudaMalloc(&q_measure, a_size * sizeof(float));
        cudaCheckError();
        q_measure_ = std::shared_ptr<float>(q_measure, cudaFree);
    }
    // Desctructor
    __host__ __device__ ~QualityMeasure()
    {
        a_size = 0;
        q_measure = nullptr;
        q_measure_.reset();
    }
    /// buffer size
    __host__ __device__ int buffSize() { return a_size; }
    /// set and get value
    __device__ float& operator[] (const int i) { return q_measure[i]; }
    __device__ const float& operator[] (const int i) const { return q_measure[i]; }
    /// set value by address
    __device__ void set(const int pos, const float val) { q_measure[pos] = val; }
    /// compute square of twice the area of a triangle
    __device__ float area(const float3 v0, const float3 v1, const float3 v2)
    {
        const float x1 = v1.x - v0.x;
        const float x2 = v2.x - v0.x;
        const float y1 = v1.y - v0.y;
        const float y2 = v2.y - v0.y;
        const float z1 = v1.z - v0.z;
        const float z2 = v2.z - v0.z;
        const float a1 = (y1 * z2 - z1 * y2) * (y1 * z2 - z1 * y2);
        const float a2 = (x1 * z2 - z1 * x2) * (x1 * z2 - z1 * x2);
        const float a3 = (x1 * y2 - y1 * x2) * (x1 * y2 - y1 * x2);
        return sqrtf(a1 + a2 + a3);
    }
    /// computes the mean ratio quality measure for triangles
    __device__ void mean_ratio(const int pos, const float3 v0, const float3  v1, const float3 v2)
    {
        const float l0 = (v1.x - v0.x) * (v1.x - v0.x) + (v1.y - v0.y) * (v1.y - v0.y) + (v1.z - v0.z) * (v1.z - v0.z);
        const float l1 = (v2.x - v1.x) * (v2.x - v1.x) + (v2.y - v1.y) * (v2.y - v1.y) + (v2.z - v1.z) * (v2.z - v1.z);
        const float l2 = (v0.x - v2.x) * (v0.x - v2.x) + (v0.y - v2.y) * (v0.y - v2.y) + (v0.z - v2.z) * (v0.z - v2.z);
        const float l = l0 + l1 + l2;
        const float a = area(v0, v1, v2);
        //return 2.0f * sqrtf(3.0f) * a / l;
        q_measure[pos] = 2.0f * sqrtf(3.0f) * a / l;
    }
    

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The hash table to construct quadrilaterals assigned to an edge in the uniform grid
struct QuadrilateralHashTable {
    /// buffer size
    int t_size{ 0 };
    /// buffer containing the keys
    int* keyBuff{ nullptr };
    std::shared_ptr<int> keyBuff_;
	/**
     * buffer containing vertex indices 
     * four consecutive indices build a quadrilateral
     */
	int* indexBuff{ nullptr };
    std::shared_ptr<int> indexBuff_;
	
    /// constructors
    __host__ __device__ QuadrilateralHashTable() {}
    /*__host__ __device__ QuadrilateralHashTable(const QuadrilateralHashTable& q) 
    {
        this->t_size = q.t_size;
        this->keyBuff = q.keyBuff;
        this->keyBuff_ = q.keyBuff_;
        this->indexBuff = q.indexBuff;
        this->indexBuff_ = q.indexBuff_;
    }*/
    __host__ QuadrilateralHashTable(const int sz) : t_size{sz}
    {
        cudaMalloc(&keyBuff, sz * sizeof(int));
        cudaMalloc(&indexBuff, 4 * sz * sizeof(int));
        cudaCheckError();
        keyBuff_ = std::shared_ptr<int>(keyBuff, cudaFree);
        indexBuff_ = std::shared_ptr<int>(indexBuff, cudaFree);
    }
    /// desctructor, a host and a device version
    __host__ __device__ ~QuadrilateralHashTable() 
    {
        t_size = 0;
        keyBuff = nullptr;
        keyBuff_.reset();
        indexBuff = nullptr;
        indexBuff_.reset();
    }
    /// total size of keys and vertex buffers
    __host__ __device__ int size() { return t_size;  }

    /// set default values
    __device__ void init(const int pos)
    {
        keyBuff[pos] = EMPTY_BUCKET_32;
        indexBuff[4 * pos] = -1;
        indexBuff[4 * pos + 1] = -1;
        indexBuff[4 * pos + 2] = -1;
        indexBuff[4 * pos + 3] = -1;
    }
    /// add a vertex index in the vertex buffer by key
    __device__ bool addVertex(const int key, const int pos, const int v)
    {
        const int start_address = hash_function(key);
        // open hashing
        int h = start_address;
        int e = 1;
        for (int loop = 0; loop < 128; loop++) {
            int old = atomicCAS(&keyBuff[h], EMPTY_BUCKET_32, key);
            if (old == EMPTY_BUCKET_32 || old == key) 
            {
                indexBuff[4 * h + pos] = v;
                return true;
            }
            else {
                // step with quadratic probing
                h = (h + e * e) % t_size;
                e = e + 1;
                if (h == start_address) {
                    printf("ERROR: can't find free bucket for %d\n", key);
                    return false;
                }
            }
        }
        return false;
    }
    /// hash function 
    __device__ int hash_function(const int key)
    {
        return ((3 * key) % 300000007) % t_size;
    }
    /// get vertex id
    __device__ int v0(const int pos) { return indexBuff[4 * pos]; }
    __device__ int v1(const int pos) { return indexBuff[4 * pos + 1]; }
    __device__ int v2(const int pos) { return indexBuff[4 * pos + 2]; }
    __device__ int v3(const int pos) { return indexBuff[4 * pos + 3]; }
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hash table to compute halfedge twin's neighbors 
struct HalfedgeHashTable {
    /// buffers size
    int t_size{ 0 };
    /// keys
    unsigned long long* key{ nullptr };
    std::shared_ptr<unsigned long long> key_{ nullptr };
    /// twins indices
    int2* twin{ nullptr };
    std::shared_ptr<int2> twin_{ nullptr };
    /// constructors
    __host__ __device__ HalfedgeHashTable() {}
    /*__host__ __device__ HalfedgeHashTable(const HalfedgeHashTable& h)
    {
        this->t_size = h.t_size;
        this->key = h.key;
        this->key_ = h.key_;
        this->twin = h.twin;
        this->twin_ = h.twin_;
    }*/
    __host__ HalfedgeHashTable(const int sz) : t_size{ sz }
    {
        cudaMalloc(&key, t_size * sizeof(unsigned long long));
        cudaMalloc(&twin, t_size * sizeof(int2));
        cudaCheckError();
        key_ = std::shared_ptr<unsigned long long>(key, cudaFree);
        twin_ = std::shared_ptr<int2>(twin, cudaFree);
    }
    /// destructor
    __host__ __device__ ~HalfedgeHashTable() 
    {
        t_size = 0;
        key = nullptr;
        key_.reset();
        twin = nullptr;
        twin_.reset();
    }
    /// get buffer size
    __host__ __device__ int size() { return t_size;  }
    /// change buffer size
    __host__ void resize(const int sz)
    {
        if (sz > t_size)
        {
            cudaMalloc(&key, sz * sizeof(unsigned long long));
            cudaMalloc(&twin, sz * sizeof(int2));
            cudaCheckError();
            key_.reset(key, cudaFree);
            twin_.reset(twin, cudaFree);
        }
        t_size = sz;
    }
	/// init buffers
    __device__ void init(const int pos)
    {
        key[pos] = EMPTY_BUCKET_64;
        twin[pos].x = -1;
        twin[pos].y = -1;
    }
    /// add twin edges to table
    __device__ bool addHalfedgeTwins(const int v0, const int v1, const int he)
    {
        unsigned long long k = setKey(v0, v1);
        int h{ hash1(k) };
        int e{ 1 };
        for (int i = 0; i < 128; i++)
        {
            unsigned long long old = atomicCAS(&key[h], EMPTY_BUCKET_64, k);
            if (old == EMPTY_BUCKET_64 || old == k)
            {
                if (v0 < v1) twin[h].x = he;
                else twin[h].y = he;
                return true;
            }
            else
            {
                // step quadratic
                h = (h + e * e) % t_size;
                e++;
            }
        }
        return false;
    }
    /// set 64bit key from two 32bit integers
    __device__ unsigned long long setKey(const int v0, const int v1)
    {
        if (v0 < v1)
            return (static_cast<unsigned long long>(v0) << 32) | (v1 & 0xffffffffL);
        else
            return (static_cast<unsigned long long>(v1) << 32) | (v0 & 0xffffffffL);
    }
    /// simple hash function for 64bit key
    __device__ int hash1(const unsigned long long k)
    {
        return static_cast<int>(k % t_size);
    }
    /// murmur like hash function
    __device__ int hash2(const int v0, const int v1)
    {
        unsigned long long h = setKey(v0, v1);
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53L;
        h ^= h >> 33;
        return static_cast<int>(h % t_size);
    }
    /// 64bit hash function
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
    /// access twins
    __device__ int t0(const int pos) { return twin[pos].x; }
    __device__ int t1(const int pos) { return twin[pos].y; }
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hash table to collect edges from shared vertex data
struct EdgeHashTable {
    int t_size{ 0 };
    /// keys
    unsigned long long* key{ nullptr };
    std::shared_ptr<unsigned long long> key_{ nullptr };
    /// twins indices
    int2* edge{ nullptr };
    std::shared_ptr<int2> edge_{ nullptr };
    /// constructors
    __host__ __device__ EdgeHashTable() {}
    __host__ __device__ EdgeHashTable(const EdgeHashTable& h)
    {
        this->t_size = h.t_size;
        this->key = h.key;
        this->key_ = h.key_;
        this->edge = h.edge;
        this->edge_ = h.edge_;
    }
    __host__ EdgeHashTable(const int sz) : t_size{ sz }
    {
        cudaMalloc(&key, t_size * sizeof(unsigned long long));
        cudaMalloc(&edge, t_size * sizeof(int2));
        cudaCheckError();
        key_ = std::shared_ptr<unsigned long long>(key, cudaFree);
        edge_ = std::shared_ptr<int2>(edge, cudaFree);
    }
    /// destructor
    __host__ __device__ ~EdgeHashTable()
    {
        t_size = 0;
        key = nullptr;
        key_.reset();
        edge = nullptr;
        edge_.reset();
    }
    /// get buffer size
    __host__ __device__ int size() { return t_size; }
    /// change buffer size
    __host__ void resize(const int sz)
    {
        if (sz > t_size)
        {
            cudaMalloc(&key, sz * sizeof(unsigned long long));
            cudaMalloc(&edge, sz * sizeof(int2));
            cudaCheckError();
            key_.reset(key, cudaFree);
            edge_.reset(edge, cudaFree);
        }
        t_size = sz;
    }
    /// init buffers
    __device__ void init(const int pos)
    {
        key[pos] = EMPTY_BUCKET_64;
        edge[pos].x = -1;
        edge[pos].y = -1;
    }
    /// add edge to table
    __device__ bool addEdge(const int v0, const int v1)
    {
        unsigned long long k = setKey(v0, v1);
        int h{ hash1(k) };
        int e{ 1 };
        for (int i = 0; i < 128; i++)
        {
            unsigned long long old = atomicCAS(&key[h], EMPTY_BUCKET_64, k);
            if (old == EMPTY_BUCKET_64)
            {
                edge[h].x = v0;
                edge[h].y = v1;
                return true;
            }
            else if (old == k)
            {
                return true;
            }
            else
            {
                // step quadratic
                h = (h + e * e) % t_size;
                e++;
            }
        }
        return false;
    }
    /// set 64bit key from two 32bit integers
    __device__ unsigned long long setKey(const int v0, const int v1)
    {
        if (v0 < v1)
            return (static_cast<unsigned long long>(v0) << 32) | (v1 & 0xffffffffL);
        else
            return (static_cast<unsigned long long>(v1) << 32) | (v0 & 0xffffffffL);
    }
    /// simple hash function for 64bit key
    __device__ int hash1(const unsigned long long k)
    {
        return static_cast<int>(k % t_size);
    }
    /// murmur like hash function
    __device__ int hash2(const int v0, const int v1)
    {
        unsigned long long h = setKey(v0, v1);
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccdL;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53L;
        h ^= h >> 33;
        return static_cast<int>(h % t_size);
    }
    /// 64bit hash function
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
    /// access twins
    __device__ int v0(const int pos) { return edge[pos].x; }
    __device__ int v1(const int pos) { return edge[pos].y; }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The hash table to obtain a unique vertex index, required by the standard MC algorithm
struct VertexHashTable {
    /// size of buffers
    int t_size{ 0 }; 
    /// key buffer
    int* key{ nullptr };
    std::shared_ptr<int> key_{ nullptr };
    /// index buffer
    int* index{ nullptr };
    std::shared_ptr<int> index_{ nullptr };
    /// constructiors
    __host__ __device__ VertexHashTable() {}
    /*__host__ __device__ VertexHashTable(const VertexHashTable& v)
    {
        this->t_size = v.t_size;
        this->key = v.key;
        this->key_ = v.key_;
        this->index = v.index;
        this->index_ = v.index_;
    }*/
    __host__ VertexHashTable(const int sz) : t_size{sz}
    {
        cudaMalloc(&key, sz * sizeof(int));
        cudaMalloc(&index, sz * sizeof(int));
        cudaCheckError();
        key_ = std::shared_ptr<int>(key, cudaFree); 
        index_ = std::shared_ptr<int>(index, cudaFree);
    }
    __host__ __device__ ~VertexHashTable()
    {
        t_size = 0;
        key = nullptr;
        key_.reset();
        index = nullptr;
        index_.reset();
    }
    /// get size of buffers
    __host__ __device__ int size() { return t_size; }
    /// set default values
    __device__ void init(const int pos)
    {
        key[pos] = EMPTY_BUCKET_32;
        index[pos] = -1; // invalid address
    }
    /// simple hash function
    __device__ int hash(const int k)
    {
        return ((3 * k) % 300000007) % t_size;
    }
    /// a bit mix hash function
    __device__ uint hash_function(uint key)
    {
        key = (key ^ 61) ^ (key >> 16);
        key = key + (key << 3);
        key = key ^ (key >> 4);
        key = key * 0x27d4eb2d;
        key = key ^ (key >> 15);
        return key;
    }
    __device__ bool addVertex(const int k, int* addr, const int pos)
    {
        const int start_address = hash(k);
        int h = start_address;
        int e = 1;
        for (int loop = 0; loop < 128; loop++) {
            int old = atomicCAS(&key[h], EMPTY_BUCKET_32, k);
            if (old == EMPTY_BUCKET_32) {
                addr[pos] = h;
                return true; // vertex has to be created
            }
            else if (old == k)
            {
                addr[pos] = h;
                return false; // vertex already exists, use address to map index later
            }
            else {
                // step with quadratic probing
                h = (h + e * e) % t_size;
                e = e + 1;
                if (h == start_address) {
                    addr[pos] = -1;
                    return false;
                }
            }
        }
        addr[pos] = -1;
        return false; // something went wrong
    }
    /// access vertex index
    __device__ int v(const int pos) { return index[pos]; }
    __device__ void set(const int pos, const int val) { index[pos] = val;  }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vertex array
struct Vertices {
    /// size of buffers
    int a_size;
    /// vertex buffer
	float3* vertices{ nullptr };
    std::shared_ptr<float3> vertices_{ nullptr };
    /// normal buffer
	float3* normals{ nullptr };
    std::shared_ptr<float3> normals_{ nullptr };
    /// erros for point back projection
	float3* error{ nullptr };
    std::shared_ptr<float3> error_{ nullptr };
    /// atomic counter
	int* t_size{ nullptr };
    std::shared_ptr<int> t_size_{ nullptr };
    /// total number of vertices
    int nr_v{ 0 };
    /// constructions
    __host__ __device__ Vertices() {}
    __host__ __device__ Vertices(const Vertices& v)
    {
        this->a_size = v.a_size;
        this->vertices = v.vertices;
        this->normals = v.normals;
        this->error = v.error;
        this->t_size = v.t_size;
        this->nr_v = v.nr_v;
        this->vertices_ = v.vertices_;
        this->normals_ = v.normals_;
        this->error_ = v.error_;
        this->t_size_ = v.t_size_;
    }
    __host__ Vertices(const int sz) : a_size{sz}, nr_v{0}
    {
        cudaMalloc(&vertices, sz * sizeof(float3));
        cudaMalloc(&normals, sz * sizeof(float3));
        cudaMalloc(&error, sz * sizeof(float3));
        cudaMalloc(&t_size, sizeof(int));
        cudaMemset(t_size, 0, sizeof(int));
        vertices_ = std::shared_ptr<float3>(vertices, cudaFree);
        normals_ = std::shared_ptr<float3>(normals, cudaFree);
        error_ = std::shared_ptr<float3>(error, cudaFree);
        t_size_ = std::shared_ptr<int>(t_size, cudaFree);
    }
    /// destructor
    __host__ __device__ ~Vertices() 
    {
        vertices_.reset();
        normals_.reset();
        error_.reset();
        t_size_.reset();
        vertices = nullptr;
        normals = nullptr;
        error = nullptr;
        t_size = nullptr;
        a_size = 0;
        nr_v = 0;
    }
    /// size of buffers: capacity in C++ vector
    __host__ int capacity() { return a_size;  }
    /// number of vertices
    __host__ int size() 
    { 
        cudaMemcpy(&nr_v, t_size, sizeof(int), cudaMemcpyDeviceToHost);
        return nr_v;
    }
    /// init atomic counter to reuse buffers
    __host__ void initAtomicCounter()
    {
        cudaMemset(t_size, 0, sizeof(int));
        nr_v = 0;
    }
    /// copy data 
    __host__ void copy(Vertices& v)
    {
        nr_v = v.size();
        if (nr_v > a_size) {
            cudaMalloc(&vertices, nr_v * sizeof(float3));
            cudaMalloc(&normals, nr_v * sizeof(float3));
            cudaMalloc(&error, nr_v * sizeof(float3));
            vertices_.reset(vertices);
            normals_.reset(normals);
            error_.reset(error);
            a_size = nr_v;
        }
        cudaMemcpy((void*)vertices, v.vertices, nr_v * sizeof(float3), cudaMemcpyDeviceToDevice);
        cudaMemcpy((void*)normals, v.normals, nr_v * sizeof(float3), cudaMemcpyDeviceToDevice);
        cudaMemcpy((void*)error, v.error, nr_v * sizeof(float3), cudaMemcpyDeviceToDevice);
        cudaMemcpy((void*)t_size, v.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
    }
    /// add a vertex to buffer
    __device__ int addVertex(const float3 v, const float3 n)
    {
        const int pos = atomicAdd(t_size, 1);
        vertices[pos] = v;
        normals[pos] = n;
        return pos;
    }
    /// add a vertex to buffer, consider the approximation error
    __device__ int addVertex(const float3 v, const float3 n, const float3 e)
    {
        const int pos = atomicAdd(t_size, 1);
        vertices[pos] = v;
        normals[pos] = n;
        error[pos] = e;
        return pos;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a triangle consist of three indices
struct Triangles {
    /// size of buffer
    int a_size{ 0 }; // size of buffer
    /// triangle index buffer
    int3* triangles{ nullptr };
    std::shared_ptr<int3> triangles_{ nullptr };
    /// atomic counter
    int* t_size{ nullptr };
    std::shared_ptr<int> t_size_{ nullptr };
    /// total number of triangles, contains same value a t_size
    int nr_t{ 0 };
    /// constructors
    __host__ __device__ Triangles() {}
   /* __host__ __device__ Triangles(const Triangles& t)
    {
        this->a_size = t.a_size;
        this->triangles = t.triangles;
        this->triangles_ = t.triangles_;
        this->t_size = t.t_size;
        this->t_size_ = t.t_size_;
        this->nr_t = t.nr_t;
    }*/
    __host__ Triangles(const int sz) : a_size{sz}, nr_t{0}
    {
        cudaMalloc(&triangles, sz * sizeof(int3));
        cudaMalloc(&t_size, sizeof(int));
        cudaCheckError();
        cudaMemset(t_size, 0, sizeof(int));
        cudaCheckError();
        triangles_ = std::shared_ptr<int3>(triangles, cudaFree);
        t_size_ = std::shared_ptr<int>(t_size, cudaFree);
    }
    /// desctrucotr
    __host__ __device__ ~Triangles()
    {
        a_size = 0;
        triangles = nullptr;
        triangles_.reset();
        t_size = nullptr;
        t_size_.reset();
        nr_t = 0;
    }
    /// get size of triangle index buffer
    __host__ __device__ int capacity() { return a_size; }
    /// total number of triangles
    __host__ int size()
    {
        cudaMemcpy(&nr_t, t_size, sizeof(int), cudaMemcpyDeviceToHost);
        return nr_t;
    }
    /// set default value to atomic counter
    __host__ void initAtomicCounter()
    {
        cudaMemset(t_size, 0, sizeof(int));
        nr_t = 0;
    }
    /// add a triangle to index buffer
    __device__ void addTriangle(const int pos, const int v0, const int v1, const int v2)
    {
        triangles[pos].x = v0;
        triangles[pos].y = v1;
        triangles[pos].z = v2;
    }
    /// add triangle to index buffer
    __device__ int addTriangle(const int v0, const int v1, const int v2)
    {
        const int pos = atomicAdd(t_size, 1);
        this->triangles[pos] = { v0, v1, v2 };
        return pos;
    }
    /// compute minimum angle of a triangle based on cosine rule
    __device__ float minAngle(const float3 v0, const float3 v1, const float3 v2)
    {
        const float a = sqrtf((v1.x - v0.x) * (v1.x - v0.x) + (v1.y - v0.y) * (v1.y - v0.y) + (v1.z - v0.z) * (v1.z - v0.z));
        const float b = sqrtf((v2.x - v1.x) * (v2.x - v1.x) + (v2.y - v1.y) * (v2.y - v1.y) + (v2.z - v1.z) * (v2.z - v1.z));
        const float c = sqrtf((v0.x - v2.x) * (v0.x - v2.x) + (v0.y - v2.y) * (v0.y - v2.y) + (v0.z - v2.z) * (v0.z - v2.z));
        const float A = acosf((b * b + c * c - a * a) / (2 * b * c));
        const float B = acosf((a * a + c * c - b * b) / (2 * a * c));
        const float C = acosf((b * b + a * a - c * c) / (2 * b * a));

        return fminf(fminf(A, B), C);
    }
    /// copy data
    __host__ void copy(Triangles& t)
    {
        nr_t = t.size();
        if (nr_t > a_size) {
            // needs to resize buffers
            cudaMalloc(&triangles, nr_t * sizeof(int3));
            cudaCheckError();
            triangles_.reset(triangles, cudaFree);
            a_size = nr_t;
        }
        cudaMemcpy((void*)triangles, t.triangles, nr_t * sizeof(float3), cudaMemcpyDeviceToDevice);
        cudaMemcpy((void*)t_size, t.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a quadrilateral consist of three indices
struct Quadrilaterals {
    /// size of buffer
    int a_size{ 0 };
    /// quadrilateral index buffer
    int4* quadrilaterals{ nullptr };
    std::shared_ptr<int4> quadrilaterals_{ nullptr };
    /// atomic counter
    int* t_size{ nullptr };
    std::shared_ptr<int> t_size_{ nullptr };
    /// total number of quadrilaterals
    int nr_q{ 0 };
    /// constructors
    __host__ __device__ Quadrilaterals() {}
    /*__host__ __device__ Quadrilaterals(const Quadrilaterals& q)
    {
        this->a_size = q.a_size;
        this->quadrilaterals = q.quadrilaterals;
        this->quadrilaterals_ = q.quadrilaterals_;
        this->t_size = q.t_size;
        this->t_size_ = q.t_size_;
        this->nr_q = q.nr_q;
    }*/
    __host__ Quadrilaterals(const int sz) : a_size{ sz }, nr_q{ 0 }
    {
        cudaMalloc(&quadrilaterals, sz * sizeof(int4));
        cudaMalloc(&t_size, sizeof(int));
        cudaCheckError();
        cudaMemset(t_size, 0, sizeof(int));
        cudaCheckError();
        quadrilaterals_ = std::shared_ptr<int4>(quadrilaterals, cudaFree);
        t_size_ = std::shared_ptr<int>(t_size, cudaFree);
    }
    /// destructor
    __host__ __device__ ~Quadrilaterals()
    {
        a_size = 0;
        quadrilaterals = nullptr;
        quadrilaterals_.reset();
        t_size = nullptr;
        t_size_.reset();
        nr_q = 0;
    }
    /// buffer size
    __host__ int capacity() { return a_size; }
    /// number of quadrilaterals
    __host__ int size() 
    {
        cudaMemcpy(&nr_q, t_size, sizeof(int), cudaMemcpyDeviceToHost);
        return nr_q;
    }
    /// set default value to atomic counter
    __host__ void initAtomicCounter()
    {
        cudaMemset(t_size, 0, sizeof(int));
        nr_q = 0;
    }
    /// add a quadrilateral to index buffer
    __device__ void addQuadrilateral(const int pos, const int v0, const int v1, const int v2, const int v3)
    {
        quadrilaterals[pos].x = v0;
        quadrilaterals[pos].y = v1;
        quadrilaterals[pos].z = v2;
        quadrilaterals[pos].w = v3;
    }
    /// add a quadrilateral to index buffer
    __device__ void addQuadrilateral(const int v0, const int v1, const int v2, const int v3)
    {
        const int pos = atomicAdd(t_size, 1);
        quadrilaterals[pos] = { v0,v1,v2,v3 };
    }
    /// copy data
    __host__ void copy(Quadrilaterals& q)
    {
        nr_q = q.size();
        if (nr_q > a_size) {
            // needs to resize buffers
            cudaMalloc(&quadrilaterals, nr_q * sizeof(int4));
            quadrilaterals_.reset(quadrilaterals, cudaFree);
            a_size = nr_q;
        }
        cudaMemcpy((void*)quadrilaterals, q.quadrilaterals, nr_q * sizeof(int4), cudaMemcpyDeviceToDevice);
        cudaMemcpy((void*)t_size, q.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// a quad edge consisting of start and end vertex, just for visualization purposes, directX can't render quads
struct Edges {
    /// size of buffer
	int a_size;
    /// index buffer
	int2* edges{ nullptr };
    std::shared_ptr<int2> edges_{ nullptr };
    /// atomic counter
	int* t_size{ nullptr };
    std::shared_ptr<int> t_size_{ nullptr };
    /// total number of edges
    int nr_e{ 0 };
    /// constructors
    __host__ __device__ Edges() {}
    /*__host__ __device__ Edges(const Edges& e) 
    {
        this->a_size = e.a_size;
        this->edges = e.edges;
        this->t_size = e.t_size;
        this->nr_e = e.nr_e;
    }*/
    __host__ Edges(const int sz) : a_size{sz}, nr_e{0}
    {
        cudaMalloc(&edges, sz * sizeof(int2));
        cudaMalloc(&t_size, sizeof(int));
        cudaCheckError();
        cudaMemset(t_size, 0, sizeof(int));
        edges_ = std::shared_ptr<int2>(edges, cudaFree);
        t_size_ = std::shared_ptr<int>(t_size, cudaFree);
    }
    /// destructor
    __host__ __device__ ~Edges()
    {
        a_size = 0;
        edges = nullptr;
        t_size = nullptr;
        edges_.reset();
        t_size_.reset();
    }
    /// get size of buffer
    __host__ int capacity() { return a_size; }
    /// get number of edges
    __host__ int size()
    {
        cudaMemcpy(&nr_e, t_size, sizeof(int), cudaMemcpyDeviceToHost);
        return nr_e;
    }
    /// set default value to atomic counter
    __host__ void initAtomicCounter()
    {
        cudaMemset(t_size, 0, sizeof(int));
        nr_e = 0;
    }
    /// add an edge
    __device__ int addEdge(const int v0, const int v1)
    {
        const int pos = atomicAdd(t_size, 1);
        edges[pos].x = v0;
        edges[pos].y = v1;
        return pos;
    }
    /// copy data
    __host__ void copy(Edges& e)
    {
        nr_e = e.size();
        cudaMemcpy((void*)edges, e.edges, nr_e * sizeof(int2), cudaMemcpyDeviceToDevice);
        cudaMemcpy((void*)t_size, e.t_size, sizeof(int), cudaMemcpyDeviceToDevice);
        cudaCheckError();
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Halfedge data structure
struct Halfedges {
    /// size of buffer
	int a_size{ 0 };
	/** halfedge int4: 
	//    he.x = origin vertex
	//    he.y = face
	//    he.z = next
	//    he.w = twin
    */
	int4* he_e{ nullptr };
    std::shared_ptr<int4> he_e_{ nullptr };
	/// total number of halfedges
    int nr_he{ 0 };
    /// constrcutors
    __host__ __device__ Halfedges() { }
    /*__host__ __device__ Halfedges(const Halfedges& he)
    {
        this->a_size = he.a_size;
        this->he_e = he.he_e;
        this->he_e_ = he.he_e_;
        this->nr_he = he.nr_he;
    }*/
    __host__ Halfedges(const int sz) : a_size{sz}, nr_he{0}
    {
        cudaMalloc(&he_e, sz * sizeof(int4));
        cudaCheckError();
        he_e_ = std::shared_ptr<int4>(he_e, cudaFree);
    }
    /// desctructor
    __host__ __device__ ~Halfedges() {
        a_size = 0;
        he_e = nullptr;
        he_e_.reset();
        nr_he = 0;
    }
    /// get buffer size
    __host__ int capacity() { return a_size;  }
    /// get total number of halfedges
    __host__ __device__ int size() { return nr_he; }
    /// set size, i.e. total number of halfedges
    __host__ bool setSize(const int sz)
    {
        if (sz > a_size)
        {
            return false;
        }
        else {
            nr_he = sz;
            return true;
        }
    }
    /// change buffer size
    __host__ void resize(const int sz)
    {
        if (sz > a_size)
        {
            a_size = sz;
            cudaMalloc(&he_e, a_size * sizeof(int4));
            cudaCheckError();
            he_e_.reset(he_e, cudaFree);
        }
        nr_he = sz;
    }
    /// add data to halfedge
    __device__ void addHalfedge(const int pos, const int v, const int f, const int he, const int t)
    {
        he_e[pos].x = v;
        he_e[pos].y = f;
        he_e[pos].z = he;
        he_e[pos].w = t;
    }
    /// add halfedge, set default value for twin
    /// twins are connected in a second kernel
    __device__ void addHalfedge(const int pos, const int v, const int f, const int he)
    {
        he_e[pos].x = v;
        he_e[pos].y = f;
        he_e[pos].z = he;
        he_e[pos].w = -1;
    }

};
struct HalfedgeVertices {
	/// size of buffer
    int a_size{ 0 };
    /// index buffer, contains index of 
    /// halfedge starting at vertex
	int* he_e{ nullptr };
    std::shared_ptr<int> he_e_{ nullptr };
    /// total number of vertices
    int nr_v;
    /// constructors
    __host__ __device__ HalfedgeVertices() {}
    /*__host__ __device__ HalfedgeVertices(const HalfedgeVertices& v)
    {
        this->a_size = v.a_size;
        this->he_e = v.he_e;
        this->he_e_ = v.he_e_;
        this->nr_v = v.nr_v;
    }*/
    __host__ HalfedgeVertices(const int sz) : a_size{sz}, nr_v{0}
    {
        cudaMalloc(&he_e, sz * sizeof(int));
        cudaCheckError();
        he_e_ = std::shared_ptr<int>(he_e, cudaFree);
    }
    /// destructor
    __host__ __device__ ~HalfedgeVertices()
    {
        a_size = 0;
        he_e = nullptr;
        he_e_.reset();
        nr_v = 0;
    }
    /// size of buffer
    __host__ int capacity() { return a_size; }
    /// total number of vertices
    __host__ int size() { return nr_v; }
    /// set size
    __host__ bool setSize(const int sz)
    {
        if (sz > a_size)
        {
            return false;
        }
        else {
            nr_v = sz;
            return true;
        }
    }
    /// change size of buffer
    __host__ void resize(const int sz)
    {
        if (sz > a_size)
        {
            a_size = sz;
            cudaMalloc(&he_e, sz * sizeof(int));
            he_e_.reset(he_e, cudaFree);
        }
        nr_v = sz;
    }
    /// add vertex
    __device__ void addVertex(const int pos, const int e)
    {
        he_e[pos] = e;
    }
};
struct HalfedgeFaces {
    /// size of buffer
	int a_size{ 0 };
    /// index buffer
	int* he_e{ nullptr };
    std::shared_ptr<int> he_e_{ nullptr };
    /// total number of faces
    int nr_f{ 0 };
    /// constructors
    __host__ __device__ HalfedgeFaces() {}
    /*__host__ __device__ HalfedgeFaces(const HalfedgeFaces& f)
    {
        this->a_size = f.a_size;
        this->he_e = f.he_e;
        this->he_e_ = f.he_e_;
        this->nr_f = f.nr_f;
    }*/
    __host__ HalfedgeFaces(const int sz) : a_size{ sz }, nr_f{ 0 }
    {
        cudaMalloc(&he_e, sz * sizeof(int));
        he_e_ = std::shared_ptr<int>(he_e, cudaFree);
    }
    /// destructor
    __host__ __device__ ~HalfedgeFaces()
    {
        a_size = 0;
        he_e = nullptr;
        he_e_.reset();
        nr_f = 0;
    }
    /// size of buffer
    __host__ int capacity() { return a_size; }
    /// total number of faces
    __host__ int size() { return nr_f; }
    /// change size of buffer, for convenience
    __host__ bool setSize(const int sz)
    {
        if (sz > a_size)
        {
            return false;
        }
        else {
            nr_f = sz;
            return true;
        }
    }
    /// change size of buffer
    __host__ void resize(const int sz)
    {
        if (sz > a_size)
        {
            a_size = sz;
            cudaMalloc(&he_e, a_size * sizeof(int));
            he_e_.reset(he_e, cudaFree);
        }
        nr_f = sz;
    }
    /// add vertex
    __device__ void addFace(const int pos, const int e)
    {
        he_e[pos] = e;
    }
};

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
    /// evaluate volume data at grid vertex
    __device__ float operator() (float* u, const int i, const int j, const int k) { return u[k * jdim * idim + j * idim + i]; }
    /// evaluate volume data at grid vertex
    __device__ float evaluate(float* u, const int i, const int j, const int k) { return u[k * jdim * idim + j * idim + i]; }
    /// trilinear interpolation of position
    __device__ float3 trilinear(const float3 p[8], const float u, const float v, const float w)
    {
        float3 po;
        po.x = (1 - w) * ((1 - v) * (p[0].x + u * (p[1].x - p[0].x)) + v * (p[2].x + u * (p[3].x - p[2].x))) + w * ((1 - v) * (p[4].x + u * (p[5].x - p[4].x)) + v * (p[6].x + u * (p[7].x - p[6].x)));
        po.y = (1 - w) * ((1 - v) * (p[0].y + u * (p[1].y - p[0].y)) + v * (p[2].y + u * (p[3].y - p[2].y))) + w * ((1 - v) * (p[4].y + u * (p[5].y - p[4].y)) + v * (p[6].y + u * (p[7].y - p[6].y)));
        po.z = (1 - w) * ((1 - v) * (p[0].z + u * (p[1].z - p[0].z)) + v * (p[2].z + u * (p[3].z - p[2].z))) + w * ((1 - v) * (p[4].z + u * (p[5].z - p[4].z)) + v * (p[6].z + u * (p[7].z - p[6].z)));
        return po;
    }

    /// trilinear interpolation of scalar
    __device__ float trilinear(const float p[8], const float u, const float v, const float w)
    {
        return (1 - w) * ((1 - v) * (p[0] + u * (p[1] - p[0])) + v * (p[2] + u * (p[3] - p[2]))) + w * ((1 - v) * (p[4] + u * (p[5] - p[4])) + v * (p[6] + u * (p[7] - p[6])));
    }
    /// compute the gradient of the scalar field by central difference
    __device__ void gradient(float3 n[8], float* u, const float s[8], const int i_index, const int j_index, const int k_index)
    {   
        int v0;
        float f = 2.f;
        auto u_index = [](const int dim, int i, float& f)
        {   
            if (i < 0) {
                f = 1.0f;
                return 0;
            }
            else if (i >= dim) {
                f = 1.0f;
                return dim-1;
            }
            else {
                f = 2.0f;
                return i;
            }
        };
        // 8 vertices
        // v0, x
        v0 = u_index(idim, i_index - 1, f);
        n[0].x = (s[1] - evaluate(u, v0, j_index, k_index)) / (f * dx);
        // v0, y
        v0 = u_index(jdim, j_index - 1, f);
        n[0].y = (s[2] - evaluate(u, i_index, v0, k_index)) / (f * dy);
        // v0, z
        v0 = u_index(kdim, k_index - 1, f);
        n[0].z = (s[4] - evaluate(u, i_index, j_index, v0)) / (f * dz);

        // v1, x
        v0 = u_index(idim, i_index + 2, f);
        n[1].x = (evaluate(u, v0, j_index, k_index) - s[0]) / (f * dx);
        // v1, y
        v0 = u_index(jdim, j_index - 1, f);
        n[1].y = (s[3] - evaluate(u, i_index + 1, v0, k_index)) / (f * dy);
        // v1, z
        v0 = u_index(kdim, k_index - 1, f);
        n[1].z = (s[5] - evaluate(u, i_index + 1, j_index, v0)) / (f * dz);

        // v2, x
        v0 = u_index(idim, i_index - 1, f);
        n[2].x = (s[3] - evaluate(u, v0, j_index + 1, k_index)) / (f * dx);
        // v2, y
        v0 = u_index(jdim, j_index + 2, f);
        n[2].y = (evaluate(u, i_index, v0, k_index) - s[0]) / (f * dy);
        // v2, z
        v0 = u_index(kdim, k_index - 1, f);
        n[2].z = (s[6] - evaluate(u, i_index, j_index + 1, v0)) / (f * dz);

        // v3, x
        v0 = u_index(idim, i_index + 2, f);
        n[3].x = (evaluate(u, v0, j_index + 1, k_index) - s[2]) / (f * dx);
        // v3, y
        v0 = u_index(jdim, j_index + 2, f);
        n[3].y = (evaluate(u, i_index + 1, v0, k_index) - s[1]) / (f * dy);
        // v3, z
        v0 = u_index(kdim, k_index - 1, f);
        n[3].z = (s[7] - evaluate(u, i_index + 1, j_index + 1, v0)) / (f * dz);

        // v4, x
        v0 = u_index(idim, i_index - 1, f);
        n[4].x = (s[5] - evaluate(u, v0, j_index, k_index + 1)) / (f * dx);
        // v4, y
        v0 = u_index(jdim, j_index - 1, f);
        n[4].y = (s[6] - evaluate(u, i_index, v0, k_index + 1)) / (f * dy);
        // v4, z
        v0 = u_index(kdim, k_index + 2, f);
        n[4].z = (evaluate(u, i_index, j_index, v0) - s[0]) / (f * dz);

        // v5, x
        v0 = u_index(idim, i_index + 2, f);
        n[5].x = (evaluate(u, v0, j_index, k_index + 1) - s[4]) / (f * dx);
        // v5, y
        v0 = u_index(jdim, j_index - 1, f);
        n[5].y = (s[7] - evaluate(u, i_index + 1, v0, k_index + 1)) / (f * dy);
        // v5, z
        v0 = u_index(kdim, k_index + 2, f);
        n[5].z = (evaluate(u, i_index + 1, j_index, v0) - s[1]) / (f * dz);

        // v6, x
        v0 = u_index(idim, i_index - 1, f);
        n[6].x = (s[7] - evaluate(u, v0, j_index + 1, k_index + 1)) / (f * dx);
        // v6, y
        v0 = u_index(jdim, j_index + 2, f);
        n[6].y = (evaluate(u, i_index, v0, k_index + 1) - s[4]) / (f * dy);
        // v6, z
        v0 = u_index(kdim, k_index + 2, f);
        n[6].z = (evaluate(u, i_index, j_index + 1, v0) - s[2]) / (f * dz);

        // v7, x
        v0 = u_index(idim, i_index + 2, f);
        n[7].x = (evaluate(u, v0, j_index + 1, k_index + 1) - s[6]) / (f * dx);
        // v7, y
        v0 = u_index(jdim, j_index + 2, f);
        n[7].y = (evaluate(u, i_index + 1, v0, k_index + 1) - s[5]) / (f * dy);
        // v7, z
        v0 = u_index(kdim, k_index + 2, f);
        n[7].z = (evaluate(u, i_index + 1, j_index + 1, v0) - s[3]) / (f * dz);
    }
    /// compute the vertices of a cell in the uniform grid by index
    __device__ void cell_vertices(float3 v[8], const int i, const int j, const int k)
    {
        v[0].x = x0 + i * dx;
        v[0].y = y0 + j * dy;
        v[0].z = z0 + k * dz;

        v[1].x = v[0].x + dx;
        v[1].y = v[0].y;
        v[1].z = v[0].z;

        v[2].x = v[0].x;
        v[2].y = v[0].y + dy;
        v[2].z = v[0].z;

        v[3].x = v[0].x + dx;
        v[3].y = v[0].y + dy;
        v[3].z = v[0].z;

        v[4].x = v[0].x;
        v[4].y = v[0].y;
        v[4].z = v[0].z + dz;

        v[5].x = v[0].x + dx;
        v[5].y = v[0].y;
        v[5].z = v[0].z + dz;

        v[6].x = v[0].x;
        v[6].y = v[0].y + dy;
        v[6].z = v[0].z + dz;

        v[7].x = v[0].x + dx;
        v[7].y = v[0].y + dy;
        v[7].z = v[0].z + dz;
    }
};
using UGrid = UniformGrid;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Standard Marching Cubes Lookup Tables
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct MarchingCubesLookupTables {
    ushort* e_pattern{ nullptr };
    std::shared_ptr<ushort> e_pattern_{ nullptr };
    char* t_pattern{ nullptr };
    std::shared_ptr<char> t_pattern_{ nullptr };
    /// constructors
    __host__ __device__ MarchingCubesLookupTables() {}
    __host__ MarchingCubesLookupTables(const std::array<ushort, 256>& e, const std::array<char, 4096>& t)
    {
        // alloc and init e_pattern
        cudaMalloc(&e_pattern, 256 * sizeof(ushort));
        cudaCheckError();
        cudaMemcpy(e_pattern, &e[0], 256 * sizeof(ushort), cudaMemcpyHostToDevice);
        cudaCheckError();
        // alloc and init t_pattern
        cudaMalloc(&t_pattern, 4096 * sizeof(char));
        cudaCheckError();
        cudaMemcpy(t_pattern, &t[0], 4096 * sizeof(char), cudaMemcpyHostToDevice);
        cudaCheckError();
        e_pattern_ = std::shared_ptr<ushort>(e_pattern, cudaFree);
        t_pattern_ = std::shared_ptr<char>(t_pattern, cudaFree);
    }
    /// desctructor
    __host__ __device__ ~MarchingCubesLookupTables()
    {
        e_pattern = nullptr;
        t_pattern = nullptr;
        e_pattern_.reset();
        t_pattern_.reset();
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Dual Marching Cubes
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct CellIntersection {
    /// mc polygons
    char* r_pattern{ nullptr };
    std::shared_ptr<char> r_pattern_;
    /// list of ambiguous cases
    uchar* t_ambig{ nullptr };
    std::shared_ptr<uchar> t_ambig_;
    /// constructor
    __host__ __device__ CellIntersection() {}
    __host__ __device__ CellIntersection(const CellIntersection& c)
    {
        this->r_pattern = c.r_pattern;
        this->t_ambig = c.t_ambig;
        this->r_pattern_ = c.r_pattern_;
        this->t_ambig_ = c.t_ambig_;
    }
    __host__ CellIntersection(const std::array<char, 4352>& p, const std::array<uchar,256>& a)
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
    __host__ __device__ ~CellIntersection()
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
    /// intersect cell, compute vertices, set edges
    __device__ void slice(const float i0, float* v_data, const int i_case, 
        const int i_index, const int j_index, const int k_index, 
        float u[8], UGrid& ugrid, 
        QuadrilateralHashTable& ht_, Vertices& v_)
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
        ugrid.gradient(n, v_data, u, i_index, j_index, k_index);
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
            float3 err{ 0,0,0 };
            float3 ni = make_float3(0, 0, 0);
            projectPoint(i0, i_index, j_index, k_index, ugrid, u, n, ei, ni, err);
            // save unique vertex
            const int v_addr = v_.addVertex(ei, ni, err);
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
    __device__ void projectPoint(const float i0, const int i, const int j, const int k, UGrid& ugrid, float f[8], float3 n[8], float3& p, float3& n_, float3& err)
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
        //err.x = fabsf(i0 - ii) / i0;
        err.x = fabsf(i0 - ii);
        const float delta = 0.1 * err.x;
        err.y = err.x;
        while (iter < max_iter && err.y > delta)
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
                err.y = err.x + 0.1;
                break;
            }
            ii = ugrid.trilinear(f, u, v, w);
            err.y = fabsf(i0 - ii);
            iter++;
        }

        if (err.y < err.x) {
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
        err.z = iter;
    }
 
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures for mesh smoothing
struct VertexMap {
    /// size of buffers
    int a_size{ 0 };
    /// total number of elements
    int nr_e{ 0 };
    /// vertex valence
    int* valence{ nullptr };
    std::shared_ptr<int> valence_{ nullptr };
    /// count elements sharing this vertex
    int* count{ nullptr };
    std::shared_ptr<int> count_{ nullptr };
    /// set vertex type
    int* type{ nullptr };
    std::shared_ptr<int> type_{ nullptr };
    /// vertex to which it will be mapped
	int* twin{ nullptr };
    std::shared_ptr<int> twin_{ nullptr };
    /// mapping addres by vertex removal
	int* map_addr{ nullptr };
    std::shared_ptr<int> map_addr_{ nullptr };
    /// constructors
    __host__ __device__ VertexMap(){ }
    /*__host__ __device__ VertexMap(const VertexMap& v)
    {
        this->a_size = v.a_size;
        this->nr_e = v.nr_e;
        this->valence = v.valence;
        this->valence_ = v.valence_;
        this->count = v.count;
        this->count_ = v.count_;
        this->type = v.type;
        this->type_ = v.type_;
        this->twin = v.twin;
        this->twin_ = v.twin_;
        this->map_addr = v.map_addr;
        this->map_addr_ = v.map_addr_;
    }*/
    __host__ VertexMap(const int size_) : a_size{ size_ }, nr_e{ size_ }
    {
        cudaMalloc(&valence, size_ * sizeof(int));
        cudaMalloc(&count, size_ * sizeof(int));
        cudaMalloc(&type, size_ * sizeof(int));
        cudaMalloc(&twin, size_ * sizeof(int));
        cudaMalloc(&map_addr, size_ * sizeof(int));
        cudaCheckError();
        valence_ = std::shared_ptr<int>(valence, cudaFree);
        count_ = std::shared_ptr<int>(count, cudaFree);
        type_ = std::shared_ptr<int>(type, cudaFree);
        map_addr_ = std::shared_ptr<int>(map_addr, cudaFree);
    }
    /// desctructor
    __host__ __device__ ~VertexMap()
    {
        a_size = 0;
        nr_e = 0;
        valence = nullptr;
        count = nullptr;
        type = nullptr;
        twin = nullptr;
        map_addr = nullptr;
        valence_.reset();
        count_.reset();
        type_.reset();
        map_addr_.reset();
    }
    /// size of buffers
    __host__ int capacity() { return a_size; }
    /// total number of elements
    __host__ __device__ int size() { return nr_e; }
    /// change number of elments, if necessary resize buffers
    __host__ void resize(const int size_)
    {
        if (size_ > a_size)
        {
            a_size = size_;
            cudaMalloc(&valence, size_ * sizeof(int));
            cudaMalloc(&count, size_ * sizeof(int));
            cudaMalloc(&type, size_ * sizeof(int));
            cudaMalloc(&twin, size_ * sizeof(int));
            cudaMalloc(&map_addr, size_ * sizeof(int));
            cudaCheckError();
            valence_.reset(valence, cudaFree);
            count_.reset(count, cudaFree);
            type_.reset(type, cudaFree);
            map_addr_.reset(map_addr, cudaFree);
        }
        nr_e = size_;
    }
    /// set default values
    __device__ void init(const int pos)
    {
        valence[pos] = 0;
        count[pos] = 0;
        type[pos] = P_NEUT;
        twin[pos] = -1;
        map_addr[pos] = -1;
    }
};

struct QuadrilateralMap {
    /// size of buffer
	int a_size{ 0 };
    /// number of elements
    int nr_q{ 0 };
    /// element type
	int* type{ nullptr };
    std::shared_ptr<int> type_{ nullptr };
    /// constructors
    __host__ __device__ QuadrilateralMap() {}
    /*__host__ __device__ QuadrilateralMap(const QuadrilateralMap& q)
    {
        this->a_size = q.a_size;
        this->nr_q = q.nr_q;
        this->type = q.type;
        this->type_ = q.type_;
    }*/
    __host__ QuadrilateralMap(const int size_) : a_size{size_}, nr_q{size_}
    {
        cudaMalloc(&type, size_ * sizeof(int));
        cudaCheckError();
        type_ = std::shared_ptr<int>(type, cudaFree);
    }
    /// destructor
    __host__ __device__ ~QuadrilateralMap()
    {
        a_size = 0;
        nr_q = 0;
        type = nullptr;
        type_.reset();
    }
    /// buffer size
    __host__ int capacity() { return a_size; }
    /// number of elements
    __host__ __device__ int size() { return nr_q; }
    /// change number of elements, if necessary resize buffers
    __host__ void resize(const int size_) 
    {
        if (size_ > a_size)
        {
            // resize buffer
            cudaMalloc(&type, size_ * sizeof(int));
            cudaCheckError();
            type_.reset(type, cudaFree);
            a_size = size_;
        }
        nr_q = size_;
    }
    __device__ void init(const int pos)
    {
        type[pos] = P_NEUT;
    }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CUDA timer
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
		std::ostringstream buf;
		buf << std::setprecision(7) << " ... time in ms: " << e_milliseconds << std::endl;
        std::cout << buf.str() << std::endl;
	}
	void __host__ print(std::string& m) {
		std::ostringstream buf;
		buf << std::setprecision(7) << " ... " << m << " time in ms: " << e_milliseconds << std::endl;
        std::cout << buf.str() << std::endl;
	}
	std::string getTime(const std::string& m) {
		std::ostringstream buf;
		buf << std::setprecision(7) << " ... " << m << " time in ms: " << e_milliseconds << std::endl;
		return buf.str();
	}

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      Compute Element Quality and Generate best triangle mesh out of a quadrilateral mesh
//      Use the MaxMin angle criterium
//      Compute triangle angles based on cosine rule
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__global__ void quadrilateral_to_triangle(Quadrilaterals q_, Vertices v_, Triangles t_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q)
        return;
    // get vertices
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    const float3 p0 = v_.vertices[v0];
    const float3 p1 = v_.vertices[v1];
    const float3 p2 = v_.vertices[v2];
    const float3 p3 = v_.vertices[v3];

    float a1_ = t_.minAngle(p0, p1, p2);
    float a2_ = t_.minAngle(p0, p2, p3);
    float b1_ = fminf(a1_, a2_);
    float b2_ = fmaxf(a1_, a2_);
    a1_ = t_.minAngle(p1, p3, p0);
    a2_ = t_.minAngle(p1, p2, p3);
    float c1_ = fminf(a1_, a2_);
    float c2_ = fmaxf(a1_, a2_);

    if (b1_ < c1_ || (b1_ == c1_ && b2_ <= c2_))
    {
        t_.addTriangle(2 * tid, v1, v3, v0);
        t_.addTriangle(2 * tid + 1, v1, v2, v3);
    }
    else 
    {
        t_.addTriangle(2 * tid, v0, v1, v2);
        t_.addTriangle(2 * tid + 1, v0, v2, v3);
    }
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute quality measure for a triangle mesh
__global__ void mean_ratio_measure(Triangles t_, Vertices v_, QualityMeasure q_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= t_.nr_t)
        return;
    // collect triangle vertices
    const int v0 = t_.triangles[tid].x;
    const int v1 = t_.triangles[tid].y;
    const int v2 = t_.triangles[tid].z;
    const float3 p0 = v_.vertices[v0];
    const float3 p1 = v_.vertices[v1];
    const float3 p2 = v_.vertices[v2];
    // compute quality measure
    q_.mean_ratio(tid, p0, p1, p2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      CUDA GLOBAL FUNCTIONS
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for quadrilaterals
__global__ void init_quadrilateral_hashtable(QuadrilateralHashTable ht_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= ht_.size())
		return;
    ht_.init(tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for vertices, this is a redundant data structure, just to differentiate
__global__ void init_vertex_hashtable(VertexHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size)
        return;
    ht_.init(tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for halfedges
__global__ void init_halfedge_hashtable(HalfedgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size)
        return;
    ht_.init(tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init hash table for edges
__global__ void init_edge_hashtable(EdgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= ht_.t_size)
        return;
    ht_.init(tid);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// init maps to smooth mesh by removing elements with valence 3-5-3-5
__global__ void init_VertexMap(VertexMap v_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= v_.size())
		return;
    v_.init(tid);
}

__global__ void init_QuadrilateralMap(QuadrilateralMap m_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= m_.a_size)
		return;
	m_.init(tid);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      DUAL MARCHING CUBES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hybrid version
//  - ambiguous cases are processed without lookup table
//  - unambiguous case are processed with the standard mc
__global__ void dual_mc(const float i0, float* v_data, const uint t_size, UGrid ugrid, CellIntersection c_, QuadrilateralHashTable ht_, Vertices v_)
{
    // use a 1d grid
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (t_size <= gl_index)
        return;

    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    if (i_index >= (ugrid.idim - 1) || j_index >= (ugrid.jdim - 1) || k_index >= (ugrid.kdim - 1))
    {
        return;
    }
    
    // scalar values at vertices
    float u[8];
    u[0] = ugrid(v_data, i_index, j_index, k_index);
    u[1] = ugrid(v_data, i_index + 1, j_index, k_index);
    u[2] = ugrid(v_data, i_index, j_index + 1, k_index);
    u[3] = ugrid(v_data, i_index + 1, j_index + 1, k_index);
    u[4] = ugrid(v_data, i_index, j_index, k_index + 1);
    u[5] = ugrid(v_data, i_index + 1, j_index, k_index + 1);
    u[6] = ugrid(v_data, i_index, j_index + 1, k_index + 1);
    u[7] = ugrid(v_data, i_index + 1, j_index + 1, k_index + 1);

    // 
    uchar i_case{ 0 };
    i_case = i_case + ((uint)(u[0] >= i0));
    i_case = i_case + ((uint)(u[1] >= i0)) * 2;
    i_case = i_case + ((uint)(u[2] >= i0)) * 4;
    i_case = i_case + ((uint)(u[3] >= i0)) * 8;
    i_case = i_case + ((uint)(u[4] >= i0)) * 16;
    i_case = i_case + ((uint)(u[5] >= i0)) * 32;
    i_case = i_case + ((uint)(u[6] >= i0)) * 64;
    i_case = i_case + ((uint)(u[7] >= i0)) * 128;

    if (i_case == 0 || i_case == 255)
        return;
    // intersect cell
    c_.slice(i0, v_data, i_case, i_index, j_index, k_index, u, ugrid, ht_, v_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      Standard MARCHING CUBES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// standard marching cubes
__global__ void standard_mc(const float i0, float* v_data, const uint t_size, UGrid ugrid, MarchingCubesLookupTables l_tables, VertexHashTable ht_, Vertices v_, Triangles t_)
{
    // use a 1d grid
    const int gl_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (gl_index >= t_size)
        return;

    const int i_index = ugrid.i_index(gl_index);
    const int j_index = ugrid.j_index(gl_index);
    const int k_index = ugrid.k_index(gl_index);
    if (i_index >= (ugrid.idim - 1) || j_index >= (ugrid.jdim - 1) || k_index >= (ugrid.kdim - 1))
    {
        return;
    }

    // scalar values at vertices
    float u[8];
    u[0] = ugrid(v_data, i_index, j_index, k_index);
    u[1] = ugrid(v_data, i_index + 1, j_index, k_index);
    u[2] = ugrid(v_data, i_index, j_index + 1, k_index);
    u[3] = ugrid(v_data, i_index + 1, j_index + 1, k_index);
    u[4] = ugrid(v_data, i_index, j_index, k_index + 1);
    u[5] = ugrid(v_data, i_index + 1, j_index, k_index + 1);
    u[6] = ugrid(v_data, i_index, j_index + 1, k_index + 1);
    u[7] = ugrid(v_data, i_index + 1, j_index + 1, k_index + 1);

    // compute case
    uchar i_case{ 0 };
    i_case = i_case + ((uint)(u[0] >= i0));
    i_case = i_case + ((uint)(u[1] >= i0)) * 2;
    i_case = i_case + ((uint)(u[2] >= i0)) * 4;
    i_case = i_case + ((uint)(u[3] >= i0)) * 8;
    i_case = i_case + ((uint)(u[4] >= i0)) * 16;
    i_case = i_case + ((uint)(u[5] >= i0)) * 32;
    i_case = i_case + ((uint)(u[6] >= i0)) * 64;
    i_case = i_case + ((uint)(u[7] >= i0)) * 128;

    if (i_case == 0 || i_case == 255)
        return;
    
    // compute cell intersection
    // compute unique edge global index
    auto e_glIndex = [](const int e, const int i_idx, const int j_idx, const int k_idx, UGrid ugrid)
    {
        const unsigned long long gei_pattern_ = 670526590282893600ull;
        const int i = i_idx + (int)((gei_pattern_ >> 5 * e) & 1); // global_edge_id[eg][0];
        const int j = j_idx + (int)((gei_pattern_ >> (5 * e + 1)) & 1); // global_edge_id[eg][1];
        const int k = k_idx + (int)((gei_pattern_ >> (5 * e + 2)) & 1); // global_edge_id[eg][2];
        const int offs = (int)((gei_pattern_ >> (5 * e + 3)) & 3);
        return (3 * ugrid.gl_index(i, j, k) + offs);
    };
    // table listing end vertices of an edge
    const unsigned char l_edges_[12]{ 16, 49, 50, 32, 84, 117, 118, 100, 64, 81, 115, 98 };
    const ushort e_ = l_tables.e_pattern[i_case];
    float3 n[8];
    ugrid.gradient(n, v_data, u, i_index, j_index, k_index);
    ushort flag{ 1 };
    int v_addr[12];
    for (int e = 0; e < 12; e++)
    {
        if (flag & e_)
        {
            const int e_id = e_glIndex(e, i_index, j_index, k_index, ugrid);
            // check if vertex was already generated
            const bool v_flag = ht_.addVertex(e_id, v_addr, e);
            if (v_flag) {
                // create vertex 
                const int v0 = (l_edges_[e] & 0xF);
                const int v1 = (l_edges_[e] >> 4) & 0xF;
                const float l = (i0 - u[v0]) / (u[v1] - u[v0]);
                const float3 vt0{ ugrid.x0 + (i_index + (v0 & 0x1)) * ugrid.dx, ugrid.y0 + (j_index + ((v0 & 0x2) >> 1)) * ugrid.dy, ugrid.z0 + (k_index + ((v0 & 0x4) >> 2)) * ugrid.dz };
                const float3 vt1{ ugrid.x0 + (i_index + (v1 & 0x1)) * ugrid.dx, ugrid.y0 + (j_index + ((v1 & 0x2) >> 1)) * ugrid.dy, ugrid.z0 + (k_index + ((v1 & 0x4) >> 2)) * ugrid.dz };
                float3 ei = make_float3(0, 0, 0);
                float3 ni = make_float3(0, 0, 0);
                ei.x = vt0.x + l * (vt1.x - vt0.x);
                ei.y = vt0.y + l * (vt1.y - vt0.y);
                ei.z = vt0.z + l * (vt1.z - vt0.z);
                ni.x = n[v0].x + l * (n[v1].x - n[v0].x);
                ni.y = n[v0].y + l * (n[v1].y - n[v0].y);
                ni.z = n[v0].z + l * (n[v1].z - n[v0].z);
                const float factor = sqrtf(ni.x * ni.x + ni.y * ni.y + ni.z * ni.z);
                ni.x /= factor;
                ni.y /= factor;
                ni.z /= factor;
                // create vertex
                const int addr_ = v_.addVertex(ei, ni);
                // map index
                ht_.set(v_addr[e], addr_);
            }
        }
        flag <<= 1;
    }
    // construct triangles
    char* t_pattern = &l_tables.t_pattern[ 16 * i_case ];
    for (int t = 0; t < 16; t += 3)
    {
        if (t_pattern[t] == -1)
            break;
        // add triangle to list
        const int v0 = v_addr[static_cast<int>(t_pattern[t])];
        const int v1 = v_addr[static_cast<int>(t_pattern[t+1])];
        const int v2 = v_addr[static_cast<int>(t_pattern[t+2])];
        t_.addTriangle(v0, v1, v2);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      SHARED VERTEX or INDEXED FACE DATA STRUCTURES
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id to position in vertex array
// Construct shared vertex list of triangles
__global__ void map_triangles(VertexHashTable ht_, Triangles t_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= t_.nr_t)
        return;
    const int v0 = t_.triangles[tid].x;
    const int v1 = t_.triangles[tid].y;
    const int v2 = t_.triangles[tid].z;

    t_.triangles[tid].x = ht_.v(v0);
    t_.triangles[tid].y = ht_.v(v1);
    t_.triangles[tid].z = ht_.v(v2);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id to position in vertex array
// Construct shared vertex list for quadrilaterals
__global__ void map_quadrilaterals(QuadrilateralHashTable ht_, Quadrilaterals q_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= ht_.size())
		return;
	// collect indices
    const int v0 = ht_.v0(tid);
    const int v1 = ht_.v1(tid);
    const int v2 = ht_.v2(tid);
    const int v3 = ht_.v3(tid);
    // check for boundary elements
	if (v0 == EMPTY_BUCKET_32 || v1 == EMPTY_BUCKET_32 || v2 == EMPTY_BUCKET_32 || v3 == EMPTY_BUCKET_32)
		return;
	// save quad
    q_.addQuadrilateral(v0, v1, v2, v3);
	/*const int q_addr = atomicAdd(q_.t_size, 1);
	q_.quadrilaterals[q_addr].x = v0;
	q_.quadrilaterals[q_addr].y = v1;
	q_.quadrilaterals[q_addr].z = v2;
	q_.quadrilaterals[q_addr].w = v3;*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Map vertex global id for edges, required for rendering purposes
__global__ void map_edges(EdgeHashTable ht_, Edges e_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= ht_.size())
		return;
    if (ht_.key[tid] == EMPTY_BUCKET_64)
        return;
	// add only if quadrilateral is complete
    e_.addEdge(ht_.edge[tid].x, ht_.edge[tid].y);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Save edges into hash table
__global__ void collect_edges(Quadrilaterals q_, EdgeHashTable ht_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= q_.nr_q)
		return;
    // collect vertices
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;

	// add edges to hash table
    ht_.addEdge(v0, v1);
    ht_.addEdge(v1, v2);
    ht_.addEdge(v2, v3);
    ht_.addEdge(v3, v0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute HalfEdge data structure
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// collect halfedges, faces and vertices
// Loop over quadrilaterals
__global__ void collect_halfedge_elements(Quadrilaterals q_, Halfedges he_, HalfedgeFaces he_f, HalfedgeVertices he_v, HalfedgeHashTable ht_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= q_.nr_q)
        return;
    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;
    // get an address and save the half edge
    // there are four halfedges for each quadrilateral
    const int he_addr = 4 * tid;
    //    he.x = origin vertex
    //    he.y = face
    //    he.z = next
    //   he.w = tween
    // 1. 
    he_.he_e[he_addr].x = v0;
    he_.he_e[he_addr].y = tid;
    he_.he_e[he_addr].z = he_addr + 1;
    he_.he_e[he_addr].w = -1;
    // 2.
    he_.he_e[he_addr + 1].x = v1;
    he_.he_e[he_addr + 1].y = tid;
    he_.he_e[he_addr + 1].z = he_addr + 2;
    he_.he_e[he_addr + 1].w = -1; // twins will be computed later
    // 3.
    he_.he_e[he_addr + 2].x = v2;
    he_.he_e[he_addr + 2].y = tid;
    he_.he_e[he_addr + 2].z = he_addr + 3;
    he_.he_e[he_addr + 2].w = -1; // twins will be computed later
    // 4.
    he_.he_e[he_addr + 3].x = v3;
    he_.he_e[he_addr + 3].y = tid;
    he_.he_e[he_addr + 3].z = he_addr;
    he_.he_e[he_addr + 3].w = -1; // twins will be computed later

    // collect faces
    he_f.he_e[tid] = he_addr;
    // collect vertices, don't care about race conditions
    // last writing get the index
    he_v.he_e[v0] = he_addr; // v0 is the origin of he_addr
    he_v.he_e[v1] = he_addr + 1; // v1 is the origin of he_addr+1
    he_v.he_e[v2] = he_addr + 2; // v2 is the origin of he_addr+2
    he_v.he_e[v3] = he_addr + 3; // v3 is the origin of he_addr+3

    // set hash tables to compute connectivity
    // 1.
    ht_.addHalfedgeTwins(v0, v1, he_addr);
    // 2.
    ht_.addHalfedgeTwins(v1, v2, he_addr+1);
    // 3. 
    ht_.addHalfedgeTwins(v2, v3, he_addr + 2);
    // 4. 
    ht_.addHalfedgeTwins(v3, v0, he_addr + 3);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute half edge twins
__global__ void connect_halfedge_twins(HalfedgeHashTable et_, Halfedges he_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= et_.size())
        return;
    
    const int he_0 = et_.twin[tid].x;
    const int he_1 = et_.twin[tid].y;
    if (he_0 != -1)  he_.he_e[he_0].w = he_1; 
    if (he_1 != -1)  he_.he_e[he_1].w = he_0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute vertex valence using edges
__global__ void vertex_valence(const int nr_e, Edges e_, VertexMap v_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= nr_e)
		return;
	const int v0 = e_.edges[tid].x;
	const int v1 = e_.edges[tid].y;
	atomicAdd(&v_.valence[v0], 1);
	atomicAdd(&v_.valence[v1], 1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute vertex valence using halfedges
__global__ void vertex_valence(Halfedges e_, VertexMap v_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= e_.size())
        return;
    const int v0 = e_.he_e[tid].x;
    atomicAdd(&v_.valence[v0], 1);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Mesh simplification
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// compute vertex count, i.e. neighbor elements are of type 3-X-3-Y sharing the same valence 3 vertex
__global__ void count_vertices_P3X3Y(const int nr_q, Quadrilaterals q_, VertexMap v_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= nr_q)
		return;
	const int v0 = q_.quadrilaterals[tid].x;
	const int v1 = q_.quadrilaterals[tid].y;
	const int v2 = q_.quadrilaterals[tid].z;
	const int v3 = q_.quadrilaterals[tid].w;

	const int valence0 = v_.valence[v0];
	const int valence1 = v_.valence[v1];
	const int valence2 = v_.valence[v2];
	const int valence3 = v_.valence[v3];


	if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5) {
		atomicAdd(&v_.count[v0], 1);
		atomicAdd(&v_.count[v2], 1);
	}
	else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3) {
		atomicAdd(&v_.count[v1], 1);
		atomicAdd(&v_.count[v3], 1);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// merge vertices if possible
// - if two elements with the pattern 3-X-3-Y are neighbors, they cannot be removed, this is indicated by count, count > 1 means more than
//   one element are sharering a vertex of valence 3
// - if the element can be removed, move the vertices with valence 3 to the midpoint of the line connecting the other two vertices
//   e.g. one can allways move the element v0 (valence(v0) = 3) or v1 (valence v1 = 3)
__global__ void merge_vertices_P3X3Y(const int nr_q, Quadrilaterals q_, QuadrilateralMap qm_, VertexMap m_, Vertices v_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_q)
        return;

    const int v0 = q_.quadrilaterals[tid].x;
    const int v1 = q_.quadrilaterals[tid].y;
    const int v2 = q_.quadrilaterals[tid].z;
    const int v3 = q_.quadrilaterals[tid].w;

    const int valence0 = m_.valence[v0];
    const int valence1 = m_.valence[v1];
    const int valence2 = m_.valence[v2];
    const int valence3 = m_.valence[v3];

    if (valence0 == 3 && valence1 >= 5 && valence2 == 3 && valence3 >= 5 && m_.count[v0] == 1 && m_.count[v2] == 1)
    {
        // compute new position and normals of v0
        v_.vertices[v0].x = v_.vertices[v1].x + 0.5 * (v_.vertices[v3].x - v_.vertices[v1].x);
        v_.vertices[v0].y = v_.vertices[v1].y + 0.5 * (v_.vertices[v3].y - v_.vertices[v1].y);
        v_.vertices[v0].z = v_.vertices[v1].z + 0.5 * (v_.vertices[v3].z - v_.vertices[v1].z);
        v_.normals[v0].x = v_.normals[v1].x + 0.5 * (v_.normals[v3].x - v_.normals[v1].x);
        v_.normals[v0].y = v_.normals[v1].y + 0.5 * (v_.normals[v3].y - v_.normals[v1].y);
        v_.normals[v0].z = v_.normals[v1].z + 0.5 * (v_.normals[v3].z - v_.normals[v1].z);

        float sz_ = v_.normals[v0].x * v_.normals[v0].x + v_.normals[v0].y * v_.normals[v0].y + v_.normals[v0].z * v_.normals[v0].z;
        sz_ = sqrtf(sz_);

        v_.normals[v0].x = v_.normals[v0].x / sz_;
        v_.normals[v0].y = v_.normals[v0].y / sz_;
        v_.normals[v0].z = v_.normals[v0].z / sz_;

        // mark v2 to be removed
        m_.type[v2] = P_3535;

        // set twins of v2 to be v0, to be able to remove element later
        m_.twin[v2] = v0;

        // element has to be removed
        qm_.type[tid] = P_3535;
    }
    else if (valence0 >= 5 && valence1 == 3 && valence2 >= 5 && valence3 == 3 && m_.count[v1] == 1 && m_.count[v3] == 1)
    {
        // compute new position and normal of v1
        v_.vertices[v1].x = v_.vertices[v0].x + 0.5 * (v_.vertices[v2].x - v_.vertices[v0].x);
        v_.vertices[v1].y = v_.vertices[v0].y + 0.5 * (v_.vertices[v2].y - v_.vertices[v0].y);
        v_.vertices[v1].z = v_.vertices[v0].z + 0.5 * (v_.vertices[v2].z - v_.vertices[v0].z);
        v_.normals[v1].x = v_.normals[v0].x + 0.5 * (v_.normals[v2].x - v_.normals[v0].x);
        v_.normals[v1].y = v_.normals[v0].y + 0.5 * (v_.normals[v2].y - v_.normals[v0].y);
        v_.normals[v1].z = v_.normals[v0].z + 0.5 * (v_.normals[v2].z - v_.normals[v0].z);

        float sz_ = v_.normals[v1].x * v_.normals[v1].x + v_.normals[v1].y * v_.normals[v1].y + v_.normals[v1].z * v_.normals[v1].z;
        sz_ = sqrtf(sz_);

        v_.normals[v1].x = v_.normals[v1].x / sz_;
        v_.normals[v1].y = v_.normals[v1].y / sz_;
        v_.normals[v1].z = v_.normals[v1].z / sz_;

        // mark v3 to be removed
        m_.type[v3] = P_3535;
        // set twins, remove v3, use addres of v1
        m_.twin[v3] = v1;

        // element has to be removed
        qm_.type[tid] = P_3535;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// new vertices 
__global__ void remove_vertices_P3X3Y(const int nr_v, Vertices v_, VertexMap m_, Vertices n_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_v)
        return;
    // 
    if (m_.type[tid] != P_3535) {
        // copy vertex to new list
        const int addr = atomicAdd(n_.t_size, 1);
        n_.vertices[addr] = v_.vertices[tid];
        n_.normals[addr] = v_.normals[tid];
        // keep address for mapping
        m_.map_addr[tid] = addr;
    }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// new elements
__global__ void remove_quadrilaterals_P3X3Y(const int nr_q, Quadrilaterals q_, QuadrilateralMap qm_, VertexMap m_, Quadrilaterals n_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_q)
        return;
    // check if quadrilateral has to be removed
    if (qm_.type[tid] == P_3535)
    {
        return; // element is removed from the list
    }
    else {
        // compute new quadrilateral
        const int v0 = q_.quadrilaterals[tid].x;
        const int v1 = q_.quadrilaterals[tid].y;
        const int v2 = q_.quadrilaterals[tid].z;
        const int v3 = q_.quadrilaterals[tid].w;
        int4 nq_;
        if (m_.type[v0] == P_3535) nq_.x = m_.map_addr[m_.twin[v0]];
        else  nq_.x = m_.map_addr[v0];
        if (m_.type[v1] == P_3535) nq_.y = m_.map_addr[m_.twin[v1]];
        else  nq_.y = m_.map_addr[v1];
        if (m_.type[v2] == P_3535) nq_.z = m_.map_addr[m_.twin[v2]];
        else  nq_.z = m_.map_addr[v2];
        if (m_.type[v3] == P_3535) nq_.w = m_.map_addr[m_.twin[v3]];
        else  nq_.w = m_.map_addr[v3];

        // search an address to store element
        const int addr = atomicAdd(n_.t_size, 1);
        n_.quadrilaterals[addr] = nq_;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// use halfedge data structure pattern is 3-3-3-3, mark element as type P_3333, and vertices as type P_3333
// if neighbor element is of the same type, the element can't be removed
// mark neighbor for remove: P_NEIG
__global__ void mark_elements_P3333(const int nr_q, HalfedgeFaces q_, Halfedges he_, VertexMap vm_, QuadrilateralMap qm_)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= nr_q)
        return;
    const int e0 = q_.he_e[tid];
    const int e1 = he_.he_e[e0].z;
    const int e2 = he_.he_e[e1].z;
    const int e3 = he_.he_e[e2].z;
    const int v0 = he_.he_e[e0].x;
    const int v1 = he_.he_e[e1].x;
    const int v2 = he_.he_e[e2].x;
    const int v3 = he_.he_e[e3].x;
    const int valence0 = vm_.valence[v0];
    const int valence1 = vm_.valence[v1];
    const int valence2 = vm_.valence[v2];
    const int valence3 = vm_.valence[v3];
    if (valence0 != 3 || valence1 != 3 || valence2 != 3 || valence3 != 3) {
        return;
    }
    // construct neighborhood
    // neighbor f0
    int twin = he_.he_e[e0].w;
    if (twin == -1) return;
    const int f0 = he_.he_e[twin].y;
    int next = he_.he_e[twin].z;
    next = he_.he_e[next].z;
    const int v4 = he_.he_e[next].x;
    next = he_.he_e[next].z;
    const int v5 = he_.he_e[next].x;
    // neighbor f1
    twin = he_.he_e[e1].w;
    if (twin == -1) return;
    const int f1 = he_.he_e[twin].y;
    // neighbor f2
    twin = he_.he_e[e2].w;
    if (twin == -1) return;
    const int f2 = he_.he_e[twin].y;
    next = he_.he_e[twin].z;
    next = he_.he_e[next].z;
    const int v6 = he_.he_e[next].x;
    next = he_.he_e[next].z;
    const int v7 = he_.he_e[next].x;
    // neighbor f3
    twin = he_.he_e[e3].w;
    if (twin == -1) return;
    const int f3 = he_.he_e[twin].y;

    // compute valence
    const int valence4 = vm_.valence[v4];
    const int valence5 = vm_.valence[v5];
    const int valence6 = vm_.valence[v6];
    const int valence7 = vm_.valence[v7];

    // check if neighor is of the sampe type
    bool flag0 = (valence4 == 3 && valence5 == 3);
    bool flag1 = (valence5 == 3 && valence6 == 3);
    bool flag2 = (valence6 == 3 && valence7 == 3);
    bool flag3 = (valence7 == 3 && valence4 == 3);
    if (flag0 || flag1 || flag2 || flag3) {
        // element can't be removed
        return;
    }
    // check for special case, where removing element would 
    // generate a non-manifold mesh
    if (valence4 == 4 && valence5 == 4 && valence6 == 4 && valence7 == 4) {
        return;
    }

    // mark element as type P_3333
    qm_.type[tid] = P_3333;
    // mark neighbors, they will be removed if possible
    qm_.type[f0] = P_NEIG;
    qm_.type[f1] = P_NEIG;
    qm_.type[f2] = P_NEIG;
    qm_.type[f3] = P_NEIG;
    // vertex type
    vm_.type[v0] = P_3333;
    vm_.type[v1] = P_3333;
    vm_.type[v2] = P_3333;
    vm_.type[v3] = P_3333;
    // vertex mappring
    vm_.twin[v0] = v4;
    vm_.twin[v1] = v5;
    vm_.twin[v2] = v6;
    vm_.twin[v3] = v7;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// remove vertices
__global__ void remove_vertices_P3333(const int nr_v, Vertices v_, VertexMap vm_, Vertices n_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= nr_v)
		return;
	if (vm_.type[tid] != P_3333) 
	{
		const int addr = atomicAdd(n_.t_size, 1);
		n_.vertices[addr] = v_.vertices[tid];
		n_.normals[addr] = v_.normals[tid];
		// keep address for mapping
		vm_.map_addr[tid] = addr;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// remove elements
__global__ void remove_quadrilaterals_P3333(const int nr_q, Quadrilaterals q_, VertexMap vm_, QuadrilateralMap qm_, Quadrilaterals n_)
{
	const int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= nr_q)
		return;
	if (qm_.type[tid] != P_NEIG)
	{
		// compute new quadrilateral
		const int v0 = q_.quadrilaterals[tid].x;
		const int v1 = q_.quadrilaterals[tid].y;
		const int v2 = q_.quadrilaterals[tid].z;
		const int v3 = q_.quadrilaterals[tid].w;
		int4 nq_;
        if (qm_.type[tid] == P_3333)
		{
			nq_.x = vm_.map_addr[vm_.twin[v0]];
			nq_.y = vm_.map_addr[vm_.twin[v1]];
			nq_.z = vm_.map_addr[vm_.twin[v2]];
			nq_.w = vm_.map_addr[vm_.twin[v3]];
		}
		else {
			nq_.x = vm_.map_addr[v0];
			nq_.y = vm_.map_addr[v1];
			nq_.z = vm_.map_addr[v2];
			nq_.w = vm_.map_addr[v3];
		}
		
		// search an address to store element
		const int addr = atomicAdd(n_.t_size, 1);
		n_.quadrilaterals[addr] = nq_;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      HOST CODE
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//      Dual Marching Cubes -- Host code
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main function to process a volume data set
void 
p_dmc::DualMarchingCubes::dualMC(const float i0, Mesh& mesh)
{
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// config
    std::ofstream o_file;
    bool errorFlag = false;
    bool valenceFlag = false;
    bool simplifyFlag = true;
    bool qualityFlag = false;
    bool objFlag = false;
   
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // set Uniform Grid helper structure
    UGrid ugrid;
    ugrid.size(dims[0], dims[1], dims[2]);
	ugrid.dx = spacing[0];
	ugrid.dy = spacing[1];
	ugrid.dz = spacing[2];
	ugrid.x0 = origin[0];
	ugrid.y0 = origin[1];
	ugrid.z0 = origin[2];
    const int t_size = dims[0] * dims[1] * dims[2];

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// measure processing time using CUDA methods
	CTimer ctimer;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// processing time
    std::cout << " ... compute iso-surface\n";

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute quadrilaterals
    // 0. alloc lookup tables
	// 1. alloc memory for hash table
	// 2. alloc memory for vertices
	// 3. alloc memory for quads
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 0. Cell intersection object
    CellIntersection c_(r_pattern, t_ambig);
    cudaCheckError();

    // 1. alloc and init hash table
	const int ht_size = static_cast<int>(20e6);
    QuadrilateralHashTable ht_(ht_size);
	//cudaCheckError();
	uint b_size = MC_BLOCKSIZE;
	uint g_size{ (static_cast<uint>(ht_.size()) + b_size - 1) / b_size };
	init_quadrilateral_hashtable << < g_size, b_size >> > (ht_);
    cudaDeviceSynchronize();
	cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 2. alloc and init vertices
    int nr_v{ 0 };
    Vertices v_(15e6);
    cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 3. alloc and init quadrilaterals
    int nr_q{ 0 };
    Quadrilaterals q_(15e6);
    cudaDeviceSynchronize();
	cudaCheckError();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 4. compute iso-surface 
    ctimer.start();
	b_size = MC_BLOCKSIZE;
	g_size = (t_size + b_size - 1) / b_size;
    dual_mc << < g_size, b_size >> > (i0, m_volume, t_size, ugrid, c_, ht_, v_);
	cudaDeviceSynchronize();
    cudaCheckError();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 5. compute shared vertex list for quadrilateral mesh
	// indices of quadrilateral vertices have to be mapped to global vertex index in vertex array
	// get number of vertices
	nr_v = v_.size();
    if (nr_v == 0)
    {
        std::cout << " ERROR: no vertices\n";
        return;
    }
	
	// map quadrilateral indices
	b_size = MC_BLOCKSIZE;
	g_size = (ht_.size() + b_size - 1) / b_size;
	map_quadrilaterals <<< g_size, b_size >>> (ht_, q_);
    cudaDeviceSynchronize();
    cudaCheckError();
	// get number of quadrilaterals
	nr_q = q_.size();
    
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// compute processing time
	ctimer.stop();
    ctimer.print(std::string(" ... Dual Marching Cubes"));
    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of quadrilaterals " << nr_q << std::endl;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate halfedge data structure
    int nr_e = 4 * nr_q; // each halfedge appears only once for each element
    HalfedgeHashTable et_;
    Halfedges he_;
    HalfedgeFaces he_f;
    HalfedgeVertices he_v;
    
    auto he_mesh = [&]( const int nr_v, const int nr_q, Quadrilaterals& q, Halfedges& he, HalfedgeFaces& f, HalfedgeVertices& v, HalfedgeHashTable& ht)
    {
        const int nr_e = 4 * nr_q;
        he.resize(nr_e);
        ht.resize(2 * nr_e);
        f.resize(nr_q);
        v.resize(nr_v);
        const int b_size = MC_BLOCKSIZE;
        int g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
        init_halfedge_hashtable << < g_size, b_size >> > (ht);
        cudaDeviceSynchronize();
        cudaCheckError();
        // collect all halfedge elements
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        collect_halfedge_elements << < g_size, b_size >> > (q, he, he_f, he_v, ht);
        cudaDeviceSynchronize();
        cudaCheckError();
        // connect halfedges
        g_size = (static_cast<uint>(ht.t_size) + b_size - 1) / b_size;
        connect_halfedge_twins << < g_size, b_size >> > (ht, he);
        cudaDeviceSynchronize();
        cudaCheckError();
    };

    // compute halfedge mesh data structure
    ctimer.start();
    he_mesh(nr_v, nr_q, q_, he_, he_f, he_v, et_);
    cudaDeviceSynchronize();
    cudaCheckError();
    ctimer.stop();

    ctimer.print(std::string("... halfedge data structure"));

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Remove elements with valence pattern 3-x-3-y and 3-3-3-3
    if (simplifyFlag)
    {
        Vertices nv_(nr_v);
        Quadrilaterals nq_(nr_q);
        VertexMap vm_(nr_v);
        QuadrilateralMap qm_(nr_q);
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Remove elements with valence pattern 3-5-3-5
        // Do not remove elements which have neighbors with the same valence pattern, 
        // sharing the vertices with valence 3
        ctimer.start();
        b_size = MC_BLOCKSIZE;
        g_size = (static_cast<uint>(nr_v) + b_size - 1) / b_size;
        init_VertexMap << < g_size, b_size >> > (vm_);
        cudaCheckError();
        //
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        init_QuadrilateralMap << < g_size, b_size >> > (qm_);
        cudaDeviceSynchronize();
        cudaCheckError();
        //
        g_size = (static_cast<uint>(nr_e) + b_size - 1) / b_size;
        vertex_valence << < g_size, b_size >> > (he_, vm_);
        cudaDeviceSynchronize();
        cudaCheckError();
        //
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        count_vertices_P3X3Y << < g_size, b_size >> > (nr_q, q_, vm_);
        cudaDeviceSynchronize();
        cudaCheckError();
        //
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        merge_vertices_P3X3Y << < g_size, b_size >> > (nr_q, q_, qm_, vm_, v_);
        cudaDeviceSynchronize();
        cudaCheckError();
        //
        g_size = (static_cast<uint>(nr_v) + b_size - 1) / b_size;
        remove_vertices_P3X3Y << < g_size, b_size >> > (nr_v, v_, vm_, nv_);
        cudaDeviceSynchronize();
        cudaCheckError();
        //
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        remove_quadrilaterals_P3X3Y << < g_size, b_size >> > (nr_q, q_, qm_, vm_, nq_);
        cudaDeviceSynchronize();
        cudaCheckError();

        //// ======= Generate data structures for rendering ======= //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// For simplicity copy smoothed device data from device to device onto initial vertex
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        nr_v = nv_.size();
        nr_q = nq_.size();
        v_.copy(nv_);
        q_.copy(nq_);
        nv_.initAtomicCounter();
        nq_.initAtomicCounter();
        cudaCheckError();
        ctimer.stop();
        ctimer.print(std::string("... remove vertex valence pattern 3-5-3-5"));
        std::cout << " ... total nr. of vertices " << nr_v << std::endl;
        std::cout << " ... total nr. of quadrilaterals " << nr_q << std::endl;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Generate halfedges
        ctimer.start();
        he_mesh(nr_v, nr_q, q_, he_, he_f, he_v, et_);
        ctimer.stop();
        // print time
        ctimer.print(std::string("... halfedge structure after simplification of  pattern 3-5-3-5"));

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Remove elements with valence pattern 3-3-3-3
        // Do not remove elements which have neighbors with the same valence pattern, sharing the vertices with
        // free memory
        vm_.resize(nr_v);
        qm_.resize(nr_q);

        // initialize
        b_size = MC_BLOCKSIZE;
        g_size = (static_cast<uint>(nr_v) + b_size - 1) / b_size;
        init_VertexMap << < g_size, b_size >> > (vm_);
        cudaCheckError();

        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        init_QuadrilateralMap << < g_size, b_size >> > (qm_);
        cudaDeviceSynchronize();
        cudaCheckError();

        // compute vertex valence
        g_size = (static_cast<uint>(nr_e) + b_size - 1) / b_size;
        vertex_valence << < g_size, b_size >> > (he_, vm_);
        cudaDeviceSynchronize();
        cudaCheckError();

        // new code here
        ctimer.start();
        //markElement_HE_P3333(const int nr_q, HalfEdgeFaces q_, HalfEdges he_, VertexMap vm_, QuadrilateralMap qm_)
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        mark_elements_P3333<<< g_size, b_size >> > (nr_q, he_f, he_, vm_, qm_);
        cudaDeviceSynchronize();
        cudaCheckError();

        // remove vertices
        g_size = (static_cast<uint>(nr_v) + b_size - 1) / b_size;
        remove_vertices_P3333 << < g_size, b_size >> > (nr_v, v_, vm_, nv_);
        cudaDeviceSynchronize();
        cudaCheckError();

        // remove elements
        g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
        remove_quadrilaterals_P3333 << < g_size, b_size >> > (nr_q, q_, vm_, qm_, nq_);
        cudaDeviceSynchronize();
        cudaCheckError();

        // copy elements
        nr_v = nv_.size();
        nr_q = nq_.size();
        v_.copy(nv_);
        q_.copy(nq_);
        cudaDeviceSynchronize();
        cudaCheckError();
        ctimer.stop();;
        ctimer.print(std::string(" ... remove vertex valence pattern 3-3-3-3"));
        std::cout << " ... total nr. of vertices " << nr_v << std::endl;
        std::cout << " ... total nr. of quadrilaterals " << nr_q << std::endl;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// Generate halfedges
        ctimer.start();
        he_mesh(nr_v, nr_q, q_, he_, he_f, he_v, et_);
        cudaDeviceSynchronize();
        cudaCheckError();
        ctimer.stop();
        ctimer.print(std::string(" halfedge data structure after simplification of pattern 3-3-3-3"));
    }// flag_simplify

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate a triangle mesh by optimal subdivision of quadrilaterals into trinalges
    const int nr_t = 2 * nr_q;
    Triangles t_(nr_t);
    g_size = (static_cast<uint>(t_.a_size) + b_size - 1) / b_size;
    quadrilateral_to_triangle<<< g_size, b_size>>>(q_, v_, t_);
    cudaDeviceSynchronize();
    cudaCheckError();
    t_.nr_t = nr_t;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate edges for rendering purpose
    EdgeHashTable eht_(5 * nr_q);
    Edges qe_(5 * nr_q);
    g_size = (static_cast<uint>(eht_.size()) + b_size - 1) / b_size;
    init_edge_hashtable << < g_size, b_size >> > (eht_);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(nr_q) + b_size - 1) / b_size;
    collect_edges << < g_size, b_size >> > (q_, eht_);
    cudaDeviceSynchronize();
    cudaCheckError();
    g_size = (static_cast<uint>(eht_.size()) + b_size - 1) / b_size;
    map_edges << < g_size, b_size >> > (eht_, qe_);
    cudaDeviceSynchronize();
    cudaCheckError();
    int nr_qe = qe_.size();
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute mesh quality
    if (qualityFlag)
    {
        QualityMeasure m_(nr_t);
        g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
        mean_ratio_measure << < g_size, b_size >> > (t_, v_, m_);
        cudaDeviceSynchronize();
        cudaCheckError();
        // save quality measure data
        float* qm_data = new float[nr_t];
        cudaMemcpy(qm_data, m_.q_measure, nr_t * sizeof(float), cudaMemcpyDeviceToHost);
        // write to file
        o_file.open("./data/models/quality.bin", std::ofstream::binary);
        o_file.write((char*)& nr_t, sizeof(int));
        o_file.write((char*)qm_data, nr_t * sizeof(float));
        o_file.close();
        delete[] qm_data;
    }

    // Resume
    std::cout << " ... Total number of elements generated " << std::endl;
    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of quadrilaterals " << nr_q << std::endl;
    std::cout << " ... total nr. of triangles " << nr_t << std::endl;
    std::cout << " ... total nr. of edges " << nr_qe << std::endl;

	// Copy data from device
	// copy first elements
    mesh.resizeQuadrilaterals(nr_q);
    mesh.resizeTriangles(nr_t);
    mesh.resizeLines(nr_qe);
    cudaMemcpy(mesh.quadrilateralsData(), q_.quadrilaterals, nr_q * sizeof(int4), cudaMemcpyDeviceToHost);
    cudaMemcpy(mesh.trianglesData(), t_.triangles, nr_t * sizeof(int3), cudaMemcpyDeviceToHost); 
    cudaMemcpy(mesh.linesData(), qe_.edges, nr_qe * sizeof(int2), cudaMemcpyDeviceToHost); 

  
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// copy vertices to mesh
    // Note: vertices has to be copied, they are double in mesh and float in CUDA
    mesh.resizeVertices(nr_v);
    float3* v_array = new float3[nr_v];
    float3* n_array = new float3[nr_v];
    cudaMemcpy(v_array, v_.vertices, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    cudaMemcpy(n_array, v_.normals, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
	for (int id = 0; id < nr_v; id++) {
		// copy vertices
		Mesh::Vertex v{ v_array[id].x,v_array[id].y, v_array[id].z };
		Mesh::Normal n{ -n_array[id].x, -n_array[id].y, -n_array[id].z };
		mesh.addVertex(id, v, n);
	}
    delete[] v_array;
    delete[] n_array;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// done!
	std::cout << " ... done\n";

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// standard marching cubes
void
p_dmc::DualMarchingCubes::standardMC(const float i0, Mesh& mesh)
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //std::vector<float> h_data;
    UGrid ugrid;
    //std::array<int, 3> dims;
    //std::array<float, 3> origin;
    //std::array<float, 3> spacing;
    //readDataFromFile(i_file, dims, origin, spacing, h_data);
    ugrid.size(dims[0], dims[1], dims[2]);
    ugrid.dx = spacing[0];
    ugrid.dy = spacing[1];
    ugrid.dz = spacing[2];
    ugrid.x0 = origin[0];
    ugrid.y0 = origin[1];
    ugrid.z0 = origin[2];
    const int t_size = dims[0] * dims[1] * dims[2];
    // measure processing time
    CTimer ctimer1;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // alloc and init data structures
    VertexHashTable ht_(25000000);
    Vertices v_(10000000);
    Triangles t_(20000000);
    MarchingCubesLookupTables l_tables(e_pattern, t_pattern);

    uint b_size = MC_BLOCKSIZE;
    uint g_size = (static_cast<uint>(ht_.t_size) + b_size - 1) / b_size;
    init_vertex_hashtable << < g_size, b_size >> > (ht_);
    cudaDeviceSynchronize();
    cudaCheckError();

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // processing time
    std::cout << " ... compute isosurface\n" << std::endl;;
    ctimer1.start();

    // compute MC
    g_size = (t_size + b_size - 1) / b_size;
    standard_mc << < g_size, b_size >> > (i0, m_volume, t_size, ugrid, l_tables, ht_, v_, t_);
    cudaDeviceSynchronize();
    cudaCheckError();

    // collect tirangles
    int nr_v = size<Vertices>(v_);
    int nr_t = size<Triangles>(t_);
    t_.nr_t = nr_t;
    g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
    map_triangles << < g_size, b_size >> > (ht_, t_);
    cudaDeviceSynchronize();
    cudaCheckError();

    ctimer1.stop();
    ctimer1.print(std::string("Standard Marching Cubes"));
    std::cout << " ... total nr. of vertices " << nr_v << std::endl;
    std::cout << " ... total nr. of triangles " << nr_t << std::endl;

    // compute elment quality
    float* qm_data{ nullptr };
    
    QualityMeasure m_(nr_t);
    //allocQualityMeasure(m_, nr_t);
    g_size = (static_cast<uint>(nr_t) + b_size - 1) / b_size;
    mean_ratio_measure << < g_size, b_size >> > (t_, v_, m_);
    // save quality measure data
    qm_data = new float[nr_t];
    cudaMemcpy(qm_data, m_.q_measure, nr_t * sizeof(float), cudaMemcpyDeviceToHost);
    
    // copy data into mesh
    // Copy data from device
    float3* v_array = new float3[nr_v];
    float3* n_array = new float3[nr_v];
    int3* t_array = new int3[nr_t];

    cudaMemcpy(v_array, v_.vertices, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    cudaMemcpy(n_array, v_.normals, nr_v * sizeof(float3), cudaMemcpyDeviceToHost);
    cudaMemcpy(t_array, t_.triangles, nr_t * sizeof(int3), cudaMemcpyDeviceToHost);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate mesh structure
    mesh.resizeVertices(nr_v);
    for (int id = 0; id < nr_v; id++) {
        // copy vertices
        Mesh::Vertex v{ v_array[id].x,v_array[id].y, v_array[id].z };
        Mesh::Normal n{ -n_array[id].x, -n_array[id].y, -n_array[id].z };
        mesh.addVertex(id, v, n);
    }
    // generate two triangles for each quad, directX can render only triangles
    mesh.resizeTriangles(nr_t);
    for (int id = 0; id < nr_t; id++) {
        const int v0 = t_array[id].x;
        const int v1 = t_array[id].y;
        const int v2 = t_array[id].z;
        Mesh::Triangle t1{ v0, v1, v2 };
        mesh.addTriangle(id, t1);
    }
}


void p_dmc::DualMarchingCubes::generateData(std::array<int, 3> & dim)
{
    
    const int nx = dim[0];
    const int ny = dim[1];
    const int nz = dim[2];
    const float dx = 2.0f / (nx - 1.0f);
    const float dy = 2.0f / (ny - 1.0f);
    const float dz = 2.0f / (nz - 1.0f);
    m_dx = static_cast<float>(dx);
    m_dy = static_cast<float>(dy);
    m_dz = static_cast<float>(dz);
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    m_bbox[0] = Point{ -1, -1, -1 };
    m_bbox[1] = Point{ 1, -1, -1 };
    m_bbox[2] = Point{ -1, 1, -1 };
    m_bbox[3] = Point{ 1, 1, -1 };
    m_bbox[4] = Point{ -1, -1, 1 };
    m_bbox[5] = Point{ 1, -1, 1 };
    m_bbox[6] = Point{ -1, 1, 1 };
    m_bbox[7] = Point{ 1, 1, 1 };


    auto f_ = [](const float x_, const float y_, const float z_)
    {
        double alpha = 1.0;
        double x = (x_ + 1.0) / 2.0;
        double y = (y_ + 1.0) / 2.0;
        double z = (z_ + 1.0) / 2.0;
        x = alpha * (4.0 * x - 2.0);
        y = alpha * (4.0 * y - 2.0);
        z = alpha * (4.0 * z - 2.0);
        double v = 2 * y * (y * y - 3 * x * x) * (1 - z * z) + (x * x + y * y) * (x * x + y * y) - (9 * z * z - 1) * (1 - z * z);
        return static_cast<float>(v);
    };

    float xmin = -1.0f;
    float ymin = -1.0f;
    float zmin = -1.0f;
    size_t size_ = static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
    float* h_data = new float[size_];
    float x = xmin;
    for (int i = 0; i < nx; i++)
    {
        float y = ymin;
        for (int j = 0; j < ny; j++)
        {
            float z = zmin;
            for (int k = 0; k < nz; k++) {
                const int gl_index = k * nx * ny + j * ny + i;
                h_data[gl_index] = f_(x, y, z);
                z += dz;
            }
            y += dy;
        }
        x += dx;
    }
    // copy data to device
    cudaMalloc(&m_volume, size_ * sizeof(float));
    cudaMemcpy((void*)m_volume, (void*)h_data, size_ * sizeof(float), cudaMemcpyHostToDevice);
    cudaCheckError();
    
    // set uniform grid data
    dims[0] = nx;
    dims[1] = ny;
    dims[2] = nz;
    origin[0] = m_bbox[0][0];
    origin[1] = m_bbox[0][1];
    origin[2] = m_bbox[0][2];
    spacing[0] = dx;
    spacing[1] = dy;
    spacing[2] = dz;
}