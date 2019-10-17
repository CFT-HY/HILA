#include <iostream>
#include <string>
#include <math.h>
#include <assert.h>

/////////////////////
/// Benchmark 1
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

const int nx = 32, ny = 32;

int main(){
    lattice->setup( nx, ny );

    field<matrix<2,2,double> > matrix1;
    field<matrix<2,2,double> > matrix2;
    field<matrix<2,2,double> > matrix3;

    matrix1[ALL] = 1; 
    matrix2[ALL] = 1; 

    for( int i=0; i<100; i++){
        matrix3[ALL] = matrix1[X]*matrix2[X];
    }

    return 0;
}
