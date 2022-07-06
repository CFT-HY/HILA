#include "hila.h"
#include "catch.hpp"

using MyType = float;
struct TestComplex {
    Complex<MyType> dummy_complex;
};

TEST_CASE_METHOD(TestComplex, "Complex constructor", "[Complex]") {
    SECTION("Constructor from Complex") {
        dummy_complex.random();
        Complex<MyType> temporary_complex(dummy_complex);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Constructor from real") {
        dummy_complex = 1.0;
        Complex<MyType> temporary_complex(1.0);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Construction from zero") {
        dummy_complex = 0;
        Complex<MyType> temporary_complex(0);
        REQUIRE(temporary_complex == dummy_complex);  
    }
    SECTION("Constructor from real and imaginary") {
        dummy_complex.re = 1.0;
        dummy_complex.im = 1.0;
        Complex<MyType> temporary_complex(1.0,1.0);
        REQUIRE(temporary_complex == dummy_complex);
    }
}

TEST_CASE_METHOD(TestComplex, "Complex assignment", "[Complex]") {
    SECTION("Assignment from Complex") {
        dummy_complex.random();
        Complex<MyType> temporary_complex;
        temporary_complex = dummy_complex;
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Assignment from real") {
        dummy_complex = 1.0;
        Complex<MyType> temporary_complex(1.0);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Assignment from zero") {
        dummy_complex = 0;
        Complex<MyType> temporary_complex(0);
        REQUIRE(temporary_complex == dummy_complex);  
    }
    SECTION("Assignment from real and imaginary") {
        dummy_complex.re = 1.0;
        dummy_complex.im = 1.0;
        Complex<MyType> temporary_complex(1.0,1.0);
        REQUIRE(temporary_complex == dummy_complex);
    }
}
