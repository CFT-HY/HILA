#include "hila.h"
#include "catch.hpp"

#define N 10
#define M 10

class ArrayTest {

public:
    Array<N,M,double> dummy_array;

    void fill_dummy_array() {
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < 10; j++){
                dummy_array.e(i,j) = hila::random();
            }
        }  
    };
};

TEST_CASE_METHOD(ArrayTest, "Array constructor", "[Array]") {
    SECTION("Constructor from CoordinateVector") {
        fill_dummy_array();
        auto temporary_array(dummy_array);
        REQUIRE(temporary_array==dummy_array);
    }
    SECTION("Construcor from constant") {
        dummy_array = 1;
        Array<N,M,double> temporary_array(1);
        REQUIRE(temporary_array==dummy_array);
    }
    SECTION("Constructor from zero") {
        dummy_array = 0;
        Array<N,M,double> temporary_array(0);
        REQUIRE(temporary_array==dummy_array); 
    }
    SECTION("Constructor from initializer list") {
        dummy_array = 1;
        Array<N,M,double> temporary_array = {1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1,
                                             1,1,1,1,1,1,1,1,1,1};
        REQUIRE(temporary_array==dummy_array);
    }
}