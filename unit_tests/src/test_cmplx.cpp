#include "hila.h"
#include "catch.hpp"
#include <vector>

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
    std::vector<MyType> as_vector(Complex<MyType> z) {
        std::vector<MyType> vec({z.re, z.im});
        return vec;
    }
};

TEST_CASE_METHOD(ComplexTest, "Complex constructor", "[Complex]") {
    SECTION("Constructor from Complex") {
        fill_dummy_complex(dummy_re, dummy_im);
        Complex<MyType> temporary_complex(dummy_complex);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Constructor from real") {
        fill_dummy_complex(dummy_re, 0);
        Complex<MyType> temporary_complex(dummy_re);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Construction from zero") {
        fill_dummy_complex(0, 0);
        Complex<MyType> temporary_complex(0);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Constructor from real and imaginary") {
        fill_dummy_complex(dummy_re, dummy_im);
        Complex<MyType> temporary_complex(dummy_re, dummy_im);
        REQUIRE(temporary_complex == dummy_complex);
    }
    SECTION("Constructor from complex notation") {
        fill_dummy_complex(1, 1);
        Complex<MyType> temporary_complex((1 + 1_i));
        REQUIRE(temporary_complex == dummy_complex);
    }
}

TEST_CASE_METHOD(ComplexTest, "Complex assignment", "[Complex]") {
    SECTION("Assignment from Complex") {
        fill_dummy_complex(dummy_re, dummy_im);
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
        fill_dummy_complex(dummy_re, dummy_im);
        REQUIRE(dummy_complex.re == Approx(dummy_re));
        REQUIRE(dummy_complex.im == Approx(dummy_im));
    }
    SECTION("Assignment from Complex<A>") {
        fill_dummy_complex(dummy_re, dummy_im);
        Complex<double> temporary_complex;
        temporary_complex = dummy_complex;
        REQUIRE(temporary_complex.re == Approx(dummy_re));
        REQUIRE(temporary_complex.im == Approx(dummy_im));
    }
    SECTION("Assignment from complex notation") {
        fill_dummy_complex(1, 1);
        Complex<double> temporary_complex;
        temporary_complex = 1 + 1_i;
        REQUIRE(temporary_complex == dummy_complex);
    }
}

TEST_CASE_METHOD(ComplexTest, "Complex utilities", "[Complex]") {
    SECTION("Real and Imag function") {
        fill_dummy_complex(dummy_re, dummy_im);
        REQUIRE(dummy_complex.real() == Approx(dummy_re));
        REQUIRE(real(dummy_complex) == Approx(dummy_re));
        REQUIRE(dummy_complex.imag() == Approx(dummy_im));
        REQUIRE(imag(dummy_complex) == Approx(dummy_im));
    }
    // SECTION("Casting") {
    //     fill_dummy_complex(1.5, 1.5);
    //     Complex<int> temporary_complex;
    //     temporary_complex = dummy_complex.cast_to<int>();
    //     REQUIRE(temporary_complex == (1 + 1_i));
    // }
}

TEST_CASE_METHOD(ComplexTest, "Complex mathematical methods", "[Complex]") {
    SECTION("Complex arithmetic") {
        fill_dummy_complex(2, 2);
        REQUIRE(((1 + 1_i) += (1 + 1_i)) == dummy_complex);
        REQUIRE(((3 + 3_i) -= (1 + 1_i)) == dummy_complex);
        REQUIRE(((1 + 0_i) *= (2 + 2_i)) == dummy_complex);
        REQUIRE(((0 + 8_i) /= (2 + 2_i)) == dummy_complex);
        REQUIRE(((1 + 1_i) + (1 + 1_i)) == dummy_complex);
        REQUIRE(((3 + 3_i) - (1 + 1_i)) == dummy_complex);
        REQUIRE(((1 + 0_i) * (2 + 2_i)) == dummy_complex);
        REQUIRE(((0 + 8_i) / (2 + 2_i)) == dummy_complex);
    }
    SECTION("Constant arithmetic") {
        fill_dummy_complex(2, 0);
        REQUIRE(((1 + 0_i) += 1) == dummy_complex);
        REQUIRE(((3 + 0_i) -= 1) == dummy_complex);
        REQUIRE(((1 + 0_i) *= 2) == dummy_complex);
        REQUIRE(((4 + 0_i) /= 2) == dummy_complex);
        REQUIRE(((1 + 0_i) + 1) == dummy_complex);
        REQUIRE(((3 + 0_i) - 1) == dummy_complex);
        REQUIRE(((1 + 0_i) * 2) == dummy_complex);
        REQUIRE(((4 + 0_i) / 2) == dummy_complex);
    }
    SECTION("Complex overloaded functions and class methods") {
        INFO("Testing functions squarenorm, abs, arg, conj, dagger");
        fill_dummy_complex(1, 1);
        REQUIRE(dummy_complex.squarenorm() == 2);
        REQUIRE(squarenorm(dummy_complex) == 2);
        REQUIRE(dummy_complex.abs() == Approx(sqrt(2)));
        REQUIRE(abs(dummy_complex) == Approx(sqrt(2)));
        REQUIRE(dummy_complex.arg() == Approx(acos(-1) / 4));
        REQUIRE(arg(dummy_complex) == Approx(acos(-1) / 4));
        REQUIRE(dummy_complex.conj().im == -1);
        REQUIRE(conj(dummy_complex).im == -1);
        REQUIRE(dummy_complex.dagger().im == -1);
        REQUIRE(dagger(dummy_complex).im == -1);
        REQUIRE(dummy_complex.conj_mul(dummy_complex) == (2 + 0_i));
        REQUIRE(dummy_complex.mul_conj(dummy_complex) == (2 + 0_i));
    }

    SECTION("Complex polar conversion") {
        fill_dummy_complex(1, 1);
        Complex<MyType> temp_method;
        Complex<MyType> temp_global;
        temp_method.polar(sqrt(2), acos(-1) / 4);
        temp_global = polar(sqrt(2), acos(-1) / 4);

        REQUIRE(temp_method.im == Approx(dummy_complex.im));
        REQUIRE(temp_method.re == Approx(dummy_complex.re));

        REQUIRE(temp_global.im == Approx(dummy_complex.im));
        REQUIRE(temp_global.re == Approx(dummy_complex.re));
    }

    SECTION("Complex mathematical functions") {
        fill_dummy_complex(1, 1);
        REQUIRE(I*(dummy_complex) == (-1 + 1_i));
        REQUIRE_THAT(as_vector(exp(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.46869, 2.28736}));
        REQUIRE_THAT(as_vector(log(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.346574, 0.785398}));
        REQUIRE_THAT(as_vector(sqrt(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.09868, 0.45509}));
        REQUIRE_THAT(as_vector(cbrt(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.08422, 0.290515}));
        REQUIRE_THAT(as_vector(sin(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.29846, 0.634964}));
        REQUIRE_THAT(as_vector(cos(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.83373, -0.988898}));
        REQUIRE_THAT(as_vector(tan(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.271753, 1.08392}));
        REQUIRE_THAT(as_vector(sinh(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.634964, 1.29846}));
        REQUIRE_THAT(as_vector(cosh(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.83373, 0.988898}));
        REQUIRE_THAT(as_vector(tanh(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.08392, 0.271753}));
        REQUIRE_THAT(as_vector(atan(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.01722, 0.402359}));
        REQUIRE_THAT(as_vector(asin(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.666239, 1.06128}));
        REQUIRE_THAT(as_vector(acos(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.904557, -1.06128}));
        REQUIRE_THAT(as_vector(atanh(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.402359, 1.01722}));
        REQUIRE_THAT(as_vector(asinh(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.06128, 0.666239}));
        REQUIRE_THAT(as_vector(acosh(dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.06128, 0.904557}));
        REQUIRE_THAT(as_vector(pow(dummy_complex, dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.273957, 0.583701}));
        REQUIRE_THAT(as_vector(pow(dummy_complex, (MyType)1.5)),
                     Catch::Matchers::Approx(std::vector<MyType>{0.643594, 1.55377}));
        REQUIRE_THAT(as_vector(pow((MyType)1.5, dummy_complex)),
                     Catch::Matchers::Approx(std::vector<MyType>{1.37838, 0.591669}));
    }
}
