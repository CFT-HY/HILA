#include "test.h"
/////////////////////
/// test_case 1
/// 2D field of matrices
/// Coverage:
/// - reduction + onsites env.
/// - EVEN + ODD accessors / neighbor access
/// - field with matrix elements
/////////////////////

int main(){
    int sum = 0;
    matrix<2,2,double> a;

    test_setup();

    a.c[0][0] = 0;
    a.c[0][1] = -1;
    a.c[1][0] = 1;
    a.c[1][1] = 0;

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

    assert(sum==lattice->volume()*2);
    return 0;
}
