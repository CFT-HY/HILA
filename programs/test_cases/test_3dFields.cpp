#include <iostream>
#include <string>
#include <math.h>

/////////////////////
/// test_case 2
/// manipulation of 3d fields
/// Coverage:
/// - foralldir & onsites env
/// - operations between multiple fields
/////////////////////

#define NDIM 3

#include "../plumbing/field.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/cmplx.h"

// Direct output to stdout
std::ostream &hila::output = std::cout;

// Define the lattice global variable
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;

const int nx = 10, ny = 10, nz = 10;

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
	assert(lattice->size(2)==nz);
}

int main(){
    cmplx<double> sum = 0;
    lattice->setup(nx, ny, nz);
    checkLatticeSetup();
    field<cmplx<double>> s1, s2, s3;

    s1[ALL] = 0.0;
    s2[EVEN] = 1.0;
    s3[ODD] = 1.0;

    s1 = s2 + s3; //now all sites in s1 should be set to 1
    onsites(ALL){
        sum+=s1[X];
    }
    assert(sum.re==(double)nx*ny*nz); 

    onsites(ODD){
        foralldir(d){
            s2[X]-=s1[X+(direction)d];
        }
    }
    foralldir(d){
        onsites(ALL){

        }
    }

}
