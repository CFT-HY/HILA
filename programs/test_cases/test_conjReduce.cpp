#include "test.h"
/////////////////////
/// test_case 1
/// 2D field of matrices
/// Coverage:
/// - reduction + onsites env.
/// - EVEN + ODD accessors / neighbor access
/// - field with matrix elements
/// - NOTE: Assumes periodic boundary conditions
/////////////////////

//TODO: rng definition that works for MPI, GPU and CPU

int main(){
    seed_mersenne(4l);
    double sum = 0;
    test_setup();

    std::cout << lattice->volume() << "\n";

    field<matrix<2,2,double> > matrices;

    assert(matrices.fs==nullptr); //check that fieldstruct allocated only after assignment
    onsites(EVEN){
        matrix<2,2,double> a;
        double theta = 2.0*M_PI*mersenne(); //make a random rotation matrix at each even site
        a.c[0][0] = cos(theta);
        a.c[0][1] = -sin(theta);
        a.c[1][0] = sin(theta);
        a.c[1][1] = cos(theta);
        matrices[X] = a;
    }
    assert(matrices.fs!=nullptr);

    matrices[ODD] = matrices[X + XDOWN].conjugate(); //odd matrices are the conjugate of the previous even site in the X direction

    onsites(EVEN){
        matrices[X]*=matrices[X + XUP]; //multiply the even sites with the matrices stored in the neighboring odd site in X direction
    }
    onsites(ODD){
        matrices[X]*=matrices[X].conjugate(); //multiply odd sites by their conjugate
    }

    //now all sites should be unit matrices

    onsites(ALL){
        sum += matrices[X].trace();
    }

    assert(((int) round(sum)) == lattice->volume()*2);

    return 0;
}
