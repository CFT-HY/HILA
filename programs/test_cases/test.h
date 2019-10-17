#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>

#include "../plumbing/field.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/cmplx.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;

const int nd[4] = { 10, 10, 10, 4 };


inline void checkLatticeSetup(){
	for (int dir = 0; dir < NDIRS; dir++){
        //check that neighbor arrays are allocated
		assert(lattice->neighb[dir]!=nullptr);
        #ifdef CUDA
        assert(lattice->device_info.d_neighb[dir]!=nullptr)
        #endif
	}
    for(int dir = 0; dir < NDIM; dir++){
    	assert(lattice->size(dir)==nd[dir]);
    }
}

inline void test_setup(){
    #if NDIM==1
    lattice->setup( nd[0] );
    #elif NDIM==2
    lattice->setup( nd[0], nd[1] );
    #elif NDIM==3
    lattice->setup( nd[0], nd[1], nd[2] );
    #elif NDIM==4
    lattice->setup( nd[0], nd[1], nd[2], nd[3] );
    #endif
    checkLatticeSetup();
}
