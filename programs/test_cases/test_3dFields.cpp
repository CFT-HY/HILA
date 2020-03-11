#include "test.h"

/////////////////////
/// test_case 2
/// Coverage:
/// - directions, onsites and foralldir environments
/// - operations between fields 
/// - foralldir env inside onsites
/// - referring to an array of fields in a loop
/// - calling a function with const parameters
///   - requiring communication of a const field
/// - calling a function from a loop
/////////////////////

template<typename A, typename B, typename C>
void sum_test_function(A &a, const B &b, const C &c){
    onsites(ALL){
        a[X] = b[X] + c[X+XUP];
    }
}

template<typename T>
T test_template_function(T a){
  return 2*a;
}

element<cmplx<double>> test_nontemplate_function(element<cmplx<double>> a){
  element<cmplx<double>> b = a;
  return 2*a;
}


int main(int argc, char **argv){

    //check that you can increment a direction correctly
    direction d = XUP;
    direction d2 = (direction) (NDIRS - 2);
    #if NDIM > 1
    d=next_direction(d); 
    d2=next_direction(d2);
    assert(d==YUP);
    assert(XUP==0);
    assert(d2==XDOWN);
    #endif

    double sum = 0;
    field<cmplx<double>> s1, s2, s3;
    field<cmplx<double>> s4[3];

    test_setup(argc, argv);

    // Test field assingment
    s1 = 0.0;
    s2 = 1.0;
    s3 = 1.0;

    // Test sum and move constructor
    s1 = s2 + s3;

    onsites(ALL){
        sum+=s1[X].re;
    }
    assert(sum==2*(double)lattice->volume() && "onsites reduction");
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
    assert(sum==(double)lattice->volume() && "test setting field with parity");

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
        element<double> diff = s1[X].re - (NDIM+1);
	    sum += diff*diff;
    }
	  assert(sum==0 && "test neighbour fetch");
    

    // Test starting communication manually

    s1[EVEN] = 1.0;
    s2[EVEN] = 1.0;
    s2[ODD] = -s1[X+XUP];
    s2.start_move(XUP,ODD);

    sum = 0;
    onsites(ALL){
	    sum += s2[X].re;
    }
	  assert(sum==0);

    // Test referring to an array of fields

    s4[0] = s1;
    s4[1] = s1;

    s4[2][ALL] = s4[0][X] - s4[1][X];

    sum = 0;
    onsites(ALL){
        sum += (s4[2][X]*s4[2][X]).re;
    }
    assert(sum == 0);

    //Test function call outside loop
    s1[ALL] = 0.0;
    s2[ALL] = 1.0;
    sum_test_function( s3, s1, s2 ); //s3 = s1 + s2
    onsites(ALL){
        element<double> diff = s3[X].re - 1.0;
        sum += diff*diff;
    }
    assert(sum == 0);

    //Test function calls in loop
    s1[ALL] = 1.0;
    s2[ALL] = 1.0;
    onsites(ALL){
      s1[X] = test_template_function(s1[X]);
      s2[X] = test_nontemplate_function(s2[X]);
    }
    onsites(ALL){
        element<double> diff1 = s1[X].re - 2.0;
        element<double> diff2 = s2[X].re - 2.0;
        sum += diff1*diff1 + diff2*diff2;
    }
    assert(sum == 0);


    // Test array reduction
    field<double> dfield;
    dfield[ALL] = 1;
    double *arraysum = (double*) malloc(nd[TUP]*sizeof(double));
    onsites(ALL){
      element<coordinate_vector> l = coordinates(X);
      element<int> t = l[TUP];
      arraysum[t] += dfield[X];
    }
    for(int t=0; t<nd[TUP]; t++){
      assert(arraysum[t] == nd[XUP]*nd[YUP]*nd[ZUP]);
    }

    finishrun();
}
