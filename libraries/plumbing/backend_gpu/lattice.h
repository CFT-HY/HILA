#ifndef HILA_BACKEND_LATTICE_GPU_H_
#define HILA_BACKEND_LATTICE_GPU_H_

#include "coordinates.h"
#include "globals.h"

/// Define some global variables

#ifdef IN_HILA_GPU
#define HILA_GPU_EXTERN /* nothing */
#else
#define HILA_GPU_EXTERN extern
#endif


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

    // /// get the coordinates at a given site
    // __host__ __device__ const CoordinateVector &coordinates(unsigned idx) const {
    //     return _dev_coordinates[idx];
    // }
    // __host__ __device__ int coordinate(unsigned idx, Direction dir) const {
    //     return _dev_coordinates[idx][dir];
    // }

#else
    // Now not EVEN_SITES_FIRST

    // these defined in hila_gpu.cpp
    __device__ const CoordinateVector coordinates(unsigned idx) const;
    __device__ int coordinate(unsigned idx, Direction dir) const;

#endif

    /// setup the backend lattice data
    void setup(lattice_struct &lattice);

    void set_device_globals(const lattice_struct &lattice);
};

////////////////////////////////////////////////////////////////////////////
// Define here some globals and inline functions

HILA_GPU_EXTERN __constant__ unsigned _dev_field_alloc_size;
#ifdef EVEN_SITES_FIRST
HILA_GPU_EXTERN __constant__ CoordinateVector *_dev_coordinates;
#endif

template <typename T>
inline __device__ const CoordinateVector &get_dev_coordinates(T idx) {
     return _dev_coordinates[idx];
}

template <typename T, typename D>
inline __device__ int get_dev_coordinate(T idx, D dir) {
     return _dev_coordinates[idx][dir];
}


// make this template to avoid code generation if not needed
#ifdef EVEN_SITES_FIRST
template <typename T>
inline __device__ Parity get_dev_site_parity(T idx) {
    auto cv = get_dev_coordinates(idx);
    int s = 0;
    foralldir(d) s += cv.e(d);
    if (s % 2 == 0)
        return Parity::even;
    else
        return Parity::odd;
}
#endif

// Save "constants" lattice size and volume here
HILA_GPU_EXTERN hila::global<int64_t> _d_volume;
HILA_GPU_EXTERN hila::global<CoordinateVector> _d_size;

#ifndef EVEN_SITES_FIRST
HILA_GPU_EXTERN hila::global<CoordinateVector> _d_nodesize;
HILA_GPU_EXTERN hila::global<CoordinateVector> _d_nodemin;
HILA_GPU_EXTERN hila::global<CoordinateVector> _d_nodefactor;
#endif


// Then, define global functions loop_lattice_size() and _volume()
#pragma hila loop_function
inline int loop_lattice_size(Direction dir) {
    return _d_size()[dir];
}

#pragma hila loop_function
inline CoordinateVector loop_lattice_size(void) {
    return _d_size();
}

#pragma hila loop_function
inline int64_t loop_lattice_volume(void) {
    return _d_volume();
}

#ifndef EVEN_SITES_FIRST

#pragma hila loop_function
inline const CoordinateVector backend_lattice_struct::coordinates(unsigned idx) const {
    CoordinateVector c;
    unsigned vdiv, ndiv;

    vdiv = idx;
    for (int d = 0; d < NDIM - 1; ++d) {
        ndiv = vdiv / _d_nodesize()[d];
        c[d] = vdiv - ndiv * _d_nodesize()[d] + _d_nodemin()[d];
        vdiv = ndiv;
    }
    c[NDIM - 1] = vdiv + _d_nodemin()[NDIM - 1];

    return c;
}

#pragma hila loop_function
inline int backend_lattice_struct::coordinate(unsigned idx, Direction dir) const {
    return (idx / _d_nodefactor()[dir]) % _d_nodesize()[dir] + _d_nodemin()[dir];
}

#endif // not EVEN_SITES_FIRST


#endif