#include "hila.h"
#include "catch.hpp"

using MyType = float;

struct TestCoordinateVector {

public:
    CoordinateVector dummy_vector = {0,0,0};
    Vector<3,int> vector = {0,0,0};
};

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from CoordinateVector", "[CoordinateVectorCreation]") {
    CoordinateVector temporary_vector(dummy_vector);
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from initializer list", "[CoordinateVectorCreation]") {
    CoordinateVector temporary_vector({0,0,0});
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from vector", "[CoordinateVectorCreation]") {
    CoordinateVector temporary_vector(vector);
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from zero", "[coordinateVectorCreation]") {
    CoordinateVector temporary_vector(0);
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor from direction", "[coordinateVectorCreation]") {
    CoordinateVector temporary_vector(e_x);
    dummy_vector = {1,0,0};
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector assignment from zero", "[coordinateVectorAssignment]") {
    CoordinateVector temporary_vector;
    temporary_vector = 0;
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector assignment from coordinateVector", "[coordinateVectorAssignment]") {
    CoordinateVector temporary_vector;
    temporary_vector = dummy_vector;
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector assignment from initializer list", "[coordinateVectorAssignment]") {
    CoordinateVector temporary_vector;
    temporary_vector = {0,0,0};
    REQUIRE(temporary_vector == dummy_vector);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector indexing with int", "[CoordinateVectorIndexing]") {
    CoordinateVector temporary_vector({0,1,2});
    for (int i = 0; i < NDIM; i++) {
        REQUIRE(temporary_vector[i] == i);
    }
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector indexing with direction", "[CoordinateVectorIndexing]") {
    CoordinateVector temporary_vector({0,1,2});
    CHECK(temporary_vector[e_x] == 0);
    CHECK(temporary_vector[e_y] == 1);
    REQUIRE(temporary_vector[e_z] == 2);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector parity check")

