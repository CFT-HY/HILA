#include "hila.h"
#include "catch.hpp"

using MyType = float;
class ComplexTest {

public:

    double const dummy_re = hila::random();
    double const dummy_im = hila::random();
    Complex<MyType> dummy_complex;

    void fill_dummy_complex(MyType re, MyType im) {
        dummy_complex.re = re;
        dummy_complex.im = im;
    }
};

TEST_CASE_METHOD(ComplexTest, "Complex constructor", "[Complex]") {
    SECTION("Constructor from Complex") {
        fill_dummy_complex(dummy_re,dummy_im);
        Complex<MyType> temporary_complex(dummy_complex);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Constructor from real") {
        fill_dummy_complex(dummy_re,0);
        Complex<MyType> temporary_complex(dummy_re);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Construction from zero") {
        fill_dummy_complex(0,0);
        Complex<MyType> temporary_complex(0);
        REQUIRE(temporary_complex == dummy_complex);  
    }
    SECTION("Constructor from real and imaginary") {
        fill_dummy_complex(dummy_re,dummy_im);
        Complex<MyType> temporary_complex(dummy_re,dummy_im);
        REQUIRE(temporary_complex == dummy_complex);
    }
}

TEST_CASE_METHOD(ComplexTest, "Complex assignment", "[Complex]") {
    SECTION("Assignment from Complex") {
        fill_dummy_complex(dummy_re,dummy_im);
        Complex<MyType> temporary_complex;
        temporary_complex = dummy_complex;
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Assignment from real") {
        dummy_complex = dummy_re;
        REQUIRE(dummy_complex.re == Approx(dummy_re));
        REQUIRE(dummy_complex.im == 0);
    }
    SECTION("Assignment from zero") {
        dummy_complex = 0;
        REQUIRE(dummy_complex.re == 0);
        REQUIRE(dummy_complex.im == 0); 
    }
    SECTION("Assignment from real and imaginary") {
        fill_dummy_complex(dummy_re,dummy_im);
        REQUIRE(dummy_complex.re == Approx(dummy_re));
        REQUIRE(dummy_complex.im == Approx(dummy_im));
    }
    SECTION("Assignment from Complex<A>") {
        fill_dummy_complex(dummy_re,dummy_im);
        Complex<double> temporary_complex;
        temporary_complex = dummy_complex;
        REQUIRE(temporary_complex.re == Approx(dummy_re));
        REQUIRE(temporary_complex.im == Approx(dummy_im));    
    }
    SECTION("Assignment from complex notation") {
        fill_dummy_complex(1,1);
        Complex<double> temporary_complex;
        temporary_complex = 1 + 1_i;
        REQUIRE(temporary_complex == dummy_complex);
    }
}

TEST_CASE_METHOD(ComplexTest, "Complex mathematical methods", "[Complex]") {
    SECTION("Complex arithmetic") {
        fill_dummy_complex(2,2);
        REQUIRE(((1 + 1_i) += (1 + 1_i)) == dummy_complex);
        REQUIRE(((3 + 3_i) -= (1 + 1_i)) == dummy_complex);
        REQUIRE(((1 + 0_i) *= (2 + 2_i)) == dummy_complex);
        REQUIRE(((0 + 8_i) /= (2 + 2_i)) == dummy_complex);
        REQUIRE(((1 + 1_i) +  (1 + 1_i)) == dummy_complex);
        REQUIRE(((3 + 3_i) -  (1 + 1_i)) == dummy_complex);
        REQUIRE(((1 + 0_i) *  (2 + 2_i)) == dummy_complex);
        REQUIRE(((0 + 8_i) /  (2 + 2_i)) == dummy_complex);
    }
    SECTION("Constant arithmetic") {
        fill_dummy_complex(2,0);
        REQUIRE(((1 + 0_i) += 1) == dummy_complex);
        REQUIRE(((3 + 0_i) -= 1) == dummy_complex);
        REQUIRE(((1 + 0_i) *= 2) == dummy_complex);
        REQUIRE(((4 + 0_i) /= 2) == dummy_complex);
        REQUIRE(((1 + 0_i) +  1) == dummy_complex);
        REQUIRE(((3 + 0_i) -  1) == dummy_complex);
        REQUIRE(((1 + 0_i) *  2) == dummy_complex);
        REQUIRE(((4 + 0_i) /  2) == dummy_complex);
    }
    SECTION("Complex overloaded functions and class methods") {
        INFO("Testing functions squarenorm, abs, arg, conj, dagger");
        fill_dummy_complex(1,1);
        REQUIRE(dummy_complex.squarenorm() == 2);
        REQUIRE(squarenorm(dummy_complex) == 2);
        REQUIRE(dummy_complex.abs() == Approx(sqrt(2))); 
        REQUIRE(abs(dummy_complex) == Approx(sqrt(2)));
        REQUIRE(dummy_complex.arg() == Approx(acos(-1)/4));
        REQUIRE(arg(dummy_complex) == Approx(acos(-1)/4));
        REQUIRE(dummy_complex.conj().im == -1);
        REQUIRE(conj(dummy_complex).im == -1);
        REQUIRE(dummy_complex.dagger().im == -1);
        REQUIRE(dagger(dummy_complex).im == -1);
        REQUIRE(dummy_complex.real() == 1);
        REQUIRE(real(dummy_complex) == 1);
        REQUIRE(dummy_complex.imag() == 1);
        REQUIRE(imag(dummy_complex) == 1);
    }
    SECTION("Complex polar conversion") {
        fill_dummy_complex(1,1);
        Complex<MyType> temp_method;
        Complex<MyType> temp_global;
        temp_method.polar(sqrt(2),acos(-1)/4);
        temp_global = polar(sqrt(2), acos(-1)/4);
        
        REQUIRE(temp_method.im == Approx(dummy_complex.im));
        REQUIRE(temp_method.re == Approx(dummy_complex.re));

        REQUIRE(temp_global.im == Approx(dummy_complex.im));
        REQUIRE(temp_global.re == Approx(dummy_complex.re));
    }
    SECTION("Complex mathematical functions") {
        fill_dummy_complex(1,1);
        REQUIRE(multiply_by_i(dummy_complex) == (-1+1_i));
    }
}
