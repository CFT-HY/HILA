#ifndef _BACKEND_LATTICE_H_
#define _BACKEND_LATTICE_H_


/* Lattice related data that needs to be communicated
 * to kernels
 */
struct backend_lattice_struct {
  unsigned * d_neighb[NDIRS];
  unsigned field_alloc_size;
  int loop_begin, loop_end;
  location * d_coordinates;

  void setup(lattice_struct lattice);

  #pragma transformer loop_function
  location coordinates( unsigned idx ){
    return d_coordinates[idx];
  }
};


#endif