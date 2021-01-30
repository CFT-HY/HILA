#include "test.h"
#include<string>
#include "plumbing/param_input.h"

/////////////////////
/// test_case 1
/// Coverage:
/// - reduction + onsites env.
/// - EVEN + ODD accessors / neighbor access
/// - Field with matrix elements
/// - NOTE: Assumes periodic boundary conditions
/////////////////////

//TODO: rng definition that works for MPI, GPU and CPU

int main(int argc, char **argv){
    seed_mersenne(4l);
    double dsum = 0;
    int isum = 0;
    test_setup(argc, argv);

    // Test input class - should be possible to call before or after setup
    input a("params.txt");
    int nx = a.get("nx");
    std::string out = a.get("output");

    Field<Matrix<2,2,double> > matrices;
    Field<int> coordinate, nb_coordinate1, nb_coordinate2;
    Field<SU<4, double>> sufield;
    
    //check that fieldstruct allocated only after assignment
    assert(matrices.fs==nullptr);

    // Test that neighbours are fetched correctly
    foralldir(dir){
        onsites(ALL){
          element<CoordinateVector> l = X.coordinates();
          coordinate[X] = l[dir];
          nb_coordinate1[X] = (l[dir] + 1) % nd[dir];
        }

        nb_coordinate2[ALL] = coordinate[X+dir];

        onsites(ALL){
            element<int> diff = nb_coordinate1[X]-nb_coordinate2[X];
            isum += diff*diff;
        }
        assert(isum==0 && "Value fetched from neighbour is correct");
    }

    // If MPI is defined, check that gathers are counted correctly
    #ifdef MPI
    coordinate[ALL] = 1;
    lattices[0]->n_gather_done = 0;
    lattices[0]->n_gather_avoided = 0;

    nb_coordinate2[ALL] = coordinate[X-e_x];
    assert(lattices[0]->n_gather_done==1);
    assert(lattices[0]->n_gather_avoided==0);
    nb_coordinate2[ALL] = coordinate[X-e_x];
    assert(lattices[0]->n_gather_done==1);
    assert(lattices[0]->n_gather_avoided==1);
    #endif
    
    

    // Test matrix multiplication and neighbour fetches
    // Calculates M(X) * M.congugate(X+dir)
    onsites(ALL){
        if(disable_avx[X]==0){};
        element<Matrix<2,2,double>> a;
        element<double> theta = 2.0*M_PI*hila_random(); //make a random rotation matrix at each even site
        a.c[0][0] =  cos(theta);
        a.c[0][1] = -sin(theta);
        a.c[1][0] =  sin(theta);
        a.c[1][1] =  cos(theta);
        matrices[X] = a;
    }
    assert(matrices.fs!=nullptr);

    matrices[ODD] = matrices[X - e_x].conjugate(); //odd matrices are the conjugate of the previous even site in the X direction

    onsites(EVEN){
        matrices[X]*=matrices[X + e_x]; //multiply the even sites with the matrices stored in the neighboring odd site in X direction
    }
    onsites(ODD){
        matrices[X]*=matrices[X].conjugate(); //multiply odd sites by their conjugate
    }

    //now all sites should be unit matrices

    dsum=0;
    onsites(ALL){
        dsum += matrices[X].trace();
    }

    assert(((int) dsum) == lattice->volume()*2 && "Matrix conjugate multiplication");
    
    hila::finishrun();
}
