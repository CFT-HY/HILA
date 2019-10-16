#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>

/////////////////////
/// test_case 1
/// 2D field of matrices
/// Coverage:
/// - reduction + onsites env.
/// - EVEN + ODD accessors / neighbor access
/// - field with matrix elements
/////////////////////

#define NDIM 2

#include "../plumbing/field.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/cmplx.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;

const int nx = 10, ny = 10;

void checkLatticeSetup(){
	for (int dir = 0; dir < NDIRS; dir++){
        //check that neighbor arrays are allocated
		assert(lattice->neighb[dir]!=nullptr);
        #ifdef CUDA
        assert(lattice->device_info.d_neighb[dir]!=nullptr)
        #endif
	}
	assert(lattice->size(0)==nx);
	assert(lattice->size(1)==ny);
}

int main(){
    int sum = 0;
    matrix<2,2,double> a;
    a.c[0][0] = 0;
    a.c[0][1] = -1;
    a.c[1][0] = 1;
    a.c[1][1] = 0;

    lattice->setup( nx, ny );
    checkLatticeSetup();

    field<matrix<2,2,double> > matrices;

    assert(matrices.fs==nullptr); //check that fieldstruct allocated only after assignment
    matrices[EVEN] = a; //90 degree rotations
    matrices[ODD] = a.conjugate(); //-90 degree rotations
    assert(matrices.fs!=nullptr);

    matrices[EVEN]*=matrices[X + XUP]; //left & right rotations cancel out 
    matrices[ODD]*=matrices[X].conjugate();

    onsites(ALL){
        sum += (int) matrices[X].trace(); //reduction operation
    }

    assert(sum==nx*ny*2);
    return 0;
}
