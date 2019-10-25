#ifndef __DIRAC_H__
#define __DIRAC_H__

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../plumbing/field.h"

template<typename mtype, typename vtype>
void dirac_stagggered(
    const mtype &gauge,
    const double mass,
    vtype &v_in,
    const vtype &v_out)
{
    static field<double> eta[NDIM]; // The staggered phase
    static bool initialized = false;

    // Initialize the staggered eta field
    if(!initialized){
        foralldir(d){
            onsites(ALL){
                location l = coordinates(X);
                int sumcoord = 0;
                for(int d2=0;d2<d;d2++){
                    sumcoord += l[d];
                }
                if( sumcoord %2 ){
                    eta[d][X] = 1;
                } else {
                    eta[d][X] =-1;
                }
            }
        }
        initialized = true;
    }
    

    // Apply the mass diagonally
    v_out[ALL] = mass * v_in[X];

    foralldir(dir){
        direction odir = opp_dir( (direction)dir );
        direction odir2 = opp_dir( (direction)dir );
        // Positive directions: get the vector and multiply by matrix stored here
        v_out[ALL] += 0.5*eta[dir][X]*v_in[X+dir]*gauge[X];
        // Negative directions: get both form neighbour
        v_out[ALL] -= 0.5*eta[dir][X]*v_in[X+odir]*gauge[X+odir2].conjugate();
    }
}

#endif