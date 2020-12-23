#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_


/// Lattice related data that needs to be communicated
/// to kernels
struct backend_lattice_struct {
  unsigned * d_neighb[NDIRS];
  unsigned * d_neighb_special[NDIRS];
  unsigned field_alloc_size;
  int loop_begin, loop_end;
  coordinate_vector * d_coordinates;

  void setup(lattice_struct *lattice);

  #ifdef __CUDACC__
  __host__ __device__
  coordinate_vector coordinates( unsigned idx ){
    return d_coordinates[idx];
  }
  #endif
};


#endif