#include "hila.h"
#include "catch.hpp"

template<typename T>
class TestRealVarOps {
public:
    T dummy_variable = 1.0;
};

TEMPLATE_TEST_CASE_METHOD(TestRealVarOps, "Test real variable operations", "[TestRealVarOps]", int, double, float) {
    REQUIRE(TestRealVarOps<TestType>::dummy_variable == 1.0);
}