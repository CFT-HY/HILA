#include "test.h"

/////////////////////
/// test_case 2
/// Coverage:
/// - directions, onsites and foralldir environments
/// - operations between fields 
/// - TODO: implement test for foralldir env inside onsites
/////////////////////

int main(){

    //check that you can increment a direction correctly
    direction d = XUP;
    direction d2 = (direction) (NDIRS - 1);
    #if NDIM > 1
    d++; 
    d2++;
    assert(d==YUP);
    assert(XUP==0);
    assert(d2==XUP);
    #endif

    double sum = 0;
    field<cmplx<double>> s1, s2, s3;
    field<cmplx<double>> s4[3];

    test_setup();

    // Test field assingment
    s1 = 0.0;
    s2 = 1.0;
    s3 = 1.0;

    // Test sum and move constructor
    s1 = s2 + s3;

    // Test field assignment to field element (reduction)
    onsites(ALL){
        sum+=s1[X].re;
    }
    assert(sum==2*(double)lattice->volume());
    s1=0; s2=0; s3=0;
    sum = 0;

    // Test field-parity expressions
    s1[ALL] = 0.0;
    s2[EVEN] = 1.0;
    s3[ODD] = 1.0;

    s1[ALL] = s2[X] + s3[X];

    onsites(ALL){
        sum+=s1[X].re;
    }
    assert(sum==(double)lattice->volume());

    foralldir(d){
        onsites(ODD){
      	    s2[X] += s1[X + d];
	    }
	    onsites(EVEN){
		    s3[X] += s1[X + d];
	    }
    }

    s1[ALL] = s2[X] + s3[X];

    sum = 0;
    onsites(ALL){
        double diff = s1[X].re - (NDIM+1);
	    sum += diff*diff;
    }
	assert(sum==0);


    // Test foralldir loop in a onsites environment
    onsites(ALL){
        double csum = 0;
        foralldir(d){
            csum += coordinates(X)[d];
        }
        s1[X] = csum;
        s2[X] = 0;
        s3[X] = 0;
    }

    foralldir(d){
        onsites(ALL){
            s2[X] += s1[X + d];
        }
    }

    onsites(ALL){
        foralldir(d){
            s3[X] += s1[X + d];
        }
    }

    sum = 0;
    onsites(ALL){
        double diff = (s2[X]-s3[X]).re;
        sum += diff*diff;
    }
    assert(sum == 0);


    // Test referring to an array of fields

    s4[0] = s1;
    s4[1] = s1;

    s4[2][ALL] = s4[0][X] - s4[1][X];

    sum = 0;
    onsites(ALL){
        sum += (s4[2][X]*s4[2][X]).re;
    }
    assert(sum == 0);
}
