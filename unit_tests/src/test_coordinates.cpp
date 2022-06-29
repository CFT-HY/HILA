#include "hila.h"
#include "catch.hpp"

using MyType = float;

struct TestCoordinateVector {
    CoordinateVector dummy_vector = {0,0,0};
    Vector<3,int> vector = {0,0,0};
};

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor", "[CoordinateVectorCreation]") {
    SECTION("Constructor from CoordinateVector") {
        CoordinateVector temporary_vector(dummy_vector);
        REQUIRE(temporary_vector == dummy_vector);
    }
    SECTION("Construction from initializer list") {
        CoordinateVector temporary_vector({0,0,0});
        REQUIRE(temporary_vector == dummy_vector);
    }
    SECTION("Construction from Vector<NDIM,T>") {
        CoordinateVector temporary_vector(vector);
        REQUIRE(temporary_vector == dummy_vector);        
    }
    SECTION("Construction from zero") {
        CoordinateVector temporary_vector(0);
        REQUIRE(temporary_vector == dummy_vector);  
    }
    SECTION("Construction from direction") {
        CoordinateVector temporary_vector(e_x);
        dummy_vector = {1,0,0};
        REQUIRE(temporary_vector == dummy_vector);      
    }
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector assignment", "[coordinateVectorAssignment]") {
    CoordinateVector temporary_vector;
    SECTION("Assignment from zero") {
        temporary_vector = 0;
        REQUIRE(temporary_vector == dummy_vector);
    }
    SECTION("Assignment from CoordinateVector") {
        temporary_vector = dummy_vector;
        REQUIRE(temporary_vector == dummy_vector);
    }
    SECTION("Assignment from initializer list") {
        temporary_vector = {0,0,0};
        REQUIRE(temporary_vector == dummy_vector);       
    }
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector indexing", "[CoordinateVectorIndexing]") {
    CoordinateVector temporary_vector({0,1,2});
    SECTION("Indexing with int") {
        for (int i = 0; i < NDIM; i++) {
            REQUIRE(temporary_vector[i] == i);
        }
    }
    SECTION("Indexing with direction") {
        CHECK(temporary_vector[e_x] == 0);
        CHECK(temporary_vector[e_y] == 1);
        REQUIRE(temporary_vector[e_z] == 2);
    }
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector parity check", "[CoordinateVectorIndexing]") {
    CoordinateVector temporary_vector(0);
    CHECK(temporary_vector.parity() == Parity::even);
    REQUIRE((temporary_vector+=e_y).parity() == Parity::odd);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector translation", "[CoordinateVectorTranslation]") {

}