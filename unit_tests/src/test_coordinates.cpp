#include "hila.h"
#include "catch.hpp"

using MyType = float;

class TestCoordinateVector {

public:
    CoordinateVector dummy_vector = {0,0,0};
};

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from CoordinateVector", "[CoordinateVectorCreation]") {
    CoordinateVector temporary_vector(dummy_vector);
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from vector", "[CoordinateVectorCreation]") {
    CoordinateVector temporary_vector({0,0,0});
    REQUIRE(temporary_vector == dummy_vector);
}


