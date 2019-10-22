#include "test.h"

/////////////////////
/// test_case 2
/// manipulation of fields in all dirs
/// Coverage:
/// - directions, onsites and foralldir environments
/// - operations between fields 
/// - TODO: implement test for foralldir env inside onsites
/////////////////////

int main(){

    //check that you can increment a direction correctly
    #if NDIM > 1
    direction d = XUP;
    direction d2 = (direction) (NDIRS - 1);
    d++; 
    d2++;
    assert(d==YUP);
    assert(XUP==0);
    assert(d2==XUP);
    #endif

    cmplx<double> sum = 0;
    field<cmplx<double>> s1, s2, s3;

    test_setup();

    s1[ALL] = 0.0;
    s2[EVEN] = 1.0;
    s3[ODD] = 1.0;

    s1[ALL] = s2[X] + s3[X];

    onsites(ALL){
        sum+=s1[X];
    }
    assert(sum.re==(double)lattice->volume());

    foralldir(d){
        onsites(ODD){
      	    s2[X] += s1[X + d];
	    }
	    onsites(EVEN){
		    s3[X] += s1[X + d];
	    }
    }

    s1[ALL] = s2[X] + s3[X];

    onsites(ALL){
	cmplx<double> val = s1[X];
	    assert(val.re==NDIM+1);
    }
}
