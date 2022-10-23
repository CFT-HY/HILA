#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_

/// Lattice related data that needs to be communicated
/// to kernels
struct backend_lattice_struct {
    /// Storage for the neighbour indexes. Stored on device
    unsigned *d_neighb[NDIRS];

#ifdef SPECIAL_BOUNDARY_CONDITIONS
    /// Neighbour indexes with special boundaries. Stored on device
    unsigned *d_neighb_special[NDIRS];
#endif

    /// The full number of elements in a field, including haloes.
    /// This is necessary for structure-of-arrays -storage
    unsigned field_alloc_size;
    /// beginning and end of this loop (using lattice to communicate,
    /// which may not be the clearest choice.)
    int loop_begin, loop_end;

#ifdef EVEN_SITES_FIRST
    /// Finally a pointer to the list of coordinates, stored on device
    CoordinateVector *d_coordinates;

//#if defined(__CUDACC__) || defined(__HIPCC__)
#if defined(CUDA) || defined(HIP)

    /// get the coordinates at a given site
    __host__ __device__ const CoordinateVector &coordinates(unsigned idx) const {
        return d_coordinates[idx];
    }
    __host__ __device__ int coordinate(unsigned idx, Direction dir) const {
        return d_coordinates[idx][dir];
    }

#endif

#else  
    // Now not EVEN_SITES_FIRST

    // these defined in hila_gpu.cpp
    __device__ const CoordinateVector coordinates(unsigned idx) const;
    __device__ int coordinate(unsigned idx, Direction dir) const;

#endif

    /// setup the backend lattice data
    void setup(const lattice_struct &lattice);

};

//#if defined(__CUDACC__) || defined(__HIPCC__)
#if defined(CUDA) || defined(HIP)

// define also loop_lattice_size() and _volume() methods for cuda

__host__ __device__ int loop_lattice_size(Direction d);
__host__ __device__ CoordinateVector loop_lattice_size(void);
__host__ __device__ int64_t loop_lattice_volume(void);

#endif

#endif