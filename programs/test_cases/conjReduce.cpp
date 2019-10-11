#include <iostream>
#include <string>
#include <math.h>

/////////////////////
/// test_case 1
/// neighbor access on 2d array
/// return value 1 = Fail, 0 = Success
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

int main(){
    int nx = 8, ny = 8;
    int sum = 0;
    matrix<2,2,double> a;
    a.c[0][0] = 0; 
    a.c[0][1] = -1;
    a.c[1][0] = 1;
    a.c[1][1] = 0;

    lattice->setup( nx, ny );
    field<matrix<2,2,double> > matrices;

    matrices[EVEN] = a; //90 degree rotations
    matrices[ODD] = a.conjugate(); //-90 degree rotations

    matrices[EVEN]*=matrices[X + XUP]; //should turn to identity 
    matrices[ODD]*=matrices[X].conjugate(); //should also give identity
    
    onsites(ALL){
        sum += (int) matrices[X].trace(); //reduction operation
    }

    if(sum!=nx*ny*2){
        exit(1);
    } 

    return 0;
}