#include "test.h"

/////////////////////
/// test_case 2
/// manipulation of 3d fields
/// Coverage:
/// - foralldir & onsites env
/// - operations between multiple fields
/////////////////////

int main(){
    cmplx<double> sum = 0;
    field<cmplx<double>> s1, s2, s3;
    
    test_setup();

    s1[ALL] = 0.0;
    s2[EVEN] = 1.0;
    s3[ODD] = 1.0;

    s1[ALL] = s2[X] + s3[X];
    //s1 = s2 + s3; //now all sites in s1 should be set to 1
    onsites(ALL){
        sum+=s1[X];
    }
    assert(sum.re==(double)lattice->volume()); 

    foralldir(d){
        direction dir = (direction)d;
        onsites(ODD){
            s2[X]-=s1[X+dir];
        }
    }
    //foralldir(d){
    //    onsites(ALL){
//
    //    }
    //}

}
