#include "hila.h"
#include "catch.hpp"

struct TestCoordinateVector {
    CoordinateVector dummy_vector = {0,0,0};
    CoordinateVector x = {1,0,0};
    CoordinateVector y = {0,1,0};
    CoordinateVector z = {0,0,1};
};

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector constructor", "[CoordinateVector]") {
    SECTION("Constructor from CoordinateVector") {
        CoordinateVector temporary_vector(dummy_vector);
        REQUIRE(temporary_vector == dummy_vector);
    }
    SECTION("Construction from initializer list") {
        CoordinateVector temporary_vector({0,0,0});
        REQUIRE(temporary_vector == dummy_vector);
    }
    SECTION("Construction from Vector<NDIM,T>") {
        Vector<3,int> vector = {0,0,0};
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

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector assignment", "[coordinateVector]") {
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

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector indexing", "[CoordinateVector]") {
    CoordinateVector temporary_vector({0,1,2});
    SECTION("Indexing with int") {
        REQUIRE(temporary_vector[0] == 0);
        REQUIRE(temporary_vector[1] == 1);
        REQUIRE(temporary_vector[2] == 2);
    }
    SECTION("Indexing with direction") {
        CHECK(temporary_vector[e_x] == 0);
        CHECK(temporary_vector[e_y] == 1);
        REQUIRE(temporary_vector[e_z] == 2);
    }
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector parity check", "[CoordinateVector]") {
    CoordinateVector temporary_vector(0);
    REQUIRE(temporary_vector.parity() == Parity::even);
    REQUIRE((temporary_vector+=e_y).parity() == Parity::odd);
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector translation", "[CoordinateVector]") {
    SECTION("Translation with CoordinateVector") {
        REQUIRE((dummy_vector+=x) == x);
        REQUIRE((dummy_vector-=x) == dummy_vector);
    }
    SECTION("Translation with direction (unit vector)") {
        REQUIRE((dummy_vector+=e_x)==x);
        dummy_vector-=e_x;
        REQUIRE((dummy_vector+=e_y)==y);
        dummy_vector-=e_y;
        REQUIRE((dummy_vector+=e_z)==z);
        dummy_vector-=e_z;
    }
    SECTION("Translation with arithmetic") {
        REQUIRE((dummy_vector+x) ==   x );
        REQUIRE((dummy_vector-x) == (-x));
        REQUIRE((dummy_vector+y) ==   y );
        REQUIRE((dummy_vector-y) == (-y));
        REQUIRE((dummy_vector+z) ==   z );
        REQUIRE((dummy_vector-z) == (-z));
    }
    SECTION("Translation with direction arithmetic") {
        REQUIRE((dummy_vector+e_x) ==   x );
        REQUIRE((dummy_vector-e_x) == (-x));
        REQUIRE((dummy_vector+e_y) ==   y );
        REQUIRE((dummy_vector-e_y) == (-y));
        REQUIRE((dummy_vector+e_z) ==   z );
        REQUIRE((dummy_vector-e_z) == (-z));
        REQUIRE((e_x + e_x) == (x + x));
        REQUIRE((e_y + e_y) == (y + y));
        REQUIRE((e_z + e_z) == (z + z));
        REQUIRE((e_x - e_x) == dummy_vector );
        REQUIRE((e_y - e_y) == dummy_vector );
        REQUIRE((e_z - e_z) == dummy_vector );
    }
    SECTION("Scaling with product") {
        REQUIRE((2*e_x) == 2*x);
        REQUIRE((2*e_y) == 2*y);
        REQUIRE((2*e_z) == 2*z);
    }
}

TEST_CASE_METHOD(TestCoordinateVector, "CoordinateVector mathematical methods", "[CoordinateVector]") {
    SECTION("CoordinateVector dot product") {
        REQUIRE((x.dot(x)) == 1);
        REQUIRE((x.dot(y)) == 0);
    }
    // SECTION("CoordinateVector mod") {
    //     CoordinateVector scaled_vector = 2*x;
    //     std::cout << scaled_vector.mod(x) << '\n';
    // }
}