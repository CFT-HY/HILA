#include "hila.h"
#include "catch.hpp"
#include "catch_main.h"
#include <tuple>

using ArithmeticTypes =
    std::tuple<int, unsigned, long, unsigned long, float, double, long double>;

template <typename T>
class TestRealVarOps {
  public:
    T dummy_variable = hila::random();
};

TEMPLATE_LIST_TEST_CASE_METHOD(TestRealVarOps, "Test real variable operations",
                               "[TestRealVarOps]", ArithmeticTypes) {
    TestType temporary_value = TestRealVarOps<TestType>::dummy_variable;
    GIVEN("Arithmetic type") {
        CHECK(hila::is_arithmetic<TestType>::value);
        THEN("Test real value operations") {
            SECTION("real function") {
                REQUIRE(real(TestRealVarOps<TestType>::dummy_variable) ==
                        temporary_value);
            }
            SECTION("imag function") {
                REQUIRE(imag(TestRealVarOps<TestType>::dummy_variable) == 0);
            }
            SECTION("real function") {
                REQUIRE(conj(TestRealVarOps<TestType>::dummy_variable) ==
                        temporary_value);
            }
            SECTION("real function") {
                REQUIRE(squarenorm(TestRealVarOps<TestType>::dummy_variable) ==
                        temporary_value * temporary_value);
            }
            SECTION("real function") {
                REQUIRE(norm(TestRealVarOps<TestType>::dummy_variable) ==
                        Approx(sqrt(temporary_value * temporary_value)));
            }
        }
    }
}