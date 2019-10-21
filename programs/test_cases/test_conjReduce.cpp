#include "test.h"
/////////////////////
/// test_case 1
/// Coverage:
/// - reduction + onsites env.
/// - EVEN + ODD accessors / neighbor access
/// - field with matrix elements
/////////////////////

int main(){
    int sum = 0;
    matrix<2,2,double> a;
    field<matrix<2,2,double> > matrices;

    field<int> coordinate, nb_coordinate1, nb_coordinate2;
    test_setup();

    // Test that neighbours are fetched correctly
    foralldir(d){
        direction dir = (direction)d;
        onsites(ALL){
            location l = coordinates(X);
            coordinate[X] = l[d];
            nb_coordinate1[X] = (l[d] + 1) % nd[d];
        }

        nb_coordinate2[ALL] = coordinate[X+dir];

        onsites(ALL){
            sum += nb_coordinate1[X]-nb_coordinate2[X];
        }
        assert(sum==0); // Value fetched from neighbour is correct
    }

    a.c[0][0] = 0;
    a.c[0][1] = -1;
    a.c[1][0] = 1;
    a.c[1][1] = 0;


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
