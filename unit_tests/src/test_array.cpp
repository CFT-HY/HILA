#include "hila.h"
#include "catch.hpp"

#define N 10
#define M 10
using MyType = double;

class ArrayTest {

public:
    Array<N,M,MyType> dummy_array;

    template<typename T>
    void fill_dummy_array(T assign_value) {
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < 10; j++){
                dummy_array.e(i,j) = assign_value;
            }
        }  
    };

    template<typename T>
    Array<N,M,MyType> generate_temporary_array(T assign_value) {
        Array<N,M,MyType> temp(assign_value);
        return  temp;
    }
};

TEST_CASE_METHOD(ArrayTest, "Array constructor", "[Array]") {
    SECTION("Constructor from CoordinateVector") {
        dummy_array.random();
        auto temporary_array(dummy_array);
        REQUIRE(temporary_array==dummy_array);
    }
    SECTION("Construcor from constant") {
        dummy_array = 1;
        Array<N,M,MyType> temporary_array(1);
        REQUIRE(temporary_array==dummy_array);
    }
    SECTION("Constructor from zero") {
        dummy_array = 0;
        Array<N,M,MyType> temporary_array(0);
        REQUIRE(temporary_array==dummy_array); 
    }
    SECTION("Constructor from initializer list") {
        dummy_array = 1;
        Array<N,M,MyType> temporary_array = {1,1,1,1,1,1,1,1,1,1,
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

TEST_CASE_METHOD(ArrayTest, "Array assignment", "[Array]") {
    Array<N,M,MyType> temporary_array;
    SECTION("Assignment from vector") {
        dummy_array.random();
        temporary_array = dummy_array;
        REQUIRE(temporary_array == dummy_array);
    }
    SECTION("Assignment from scalar") {
        fill_dummy_array(1);
        temporary_array = 1;
        REQUIRE(temporary_array == dummy_array);
    }
    SECTION("Assignment from initializer list") {
        fill_dummy_array(1);
        temporary_array = {1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1,
                           1,1,1,1,1,1,1,1,1,1};
        REQUIRE(temporary_array == dummy_array);
    }
}

TEST_CASE_METHOD(ArrayTest, "Array indexing", "[Array]") {
    Array<N,M,MyType> temporary_array(0);
    fill_dummy_array(1);
    int index_i, index_j;
    temporary_array.e(2,1) = 1;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (temporary_array.e(i,j) == dummy_array.e(i,j)) {
                index_i = i;
                index_j = j;
            }
        } 
    }
    REQUIRE(index_i == 2);
    REQUIRE(index_j == 1);
}

TEST_CASE_METHOD(ArrayTest, "Array mathematical operations", "[Array]") {
    fill_dummy_array(2);
    Array<N,M,MyType> temporary_array(1);
    SECTION("Array arithmetic") {
        REQUIRE((generate_temporary_array(1) += generate_temporary_array(1)) == dummy_array);
        REQUIRE((generate_temporary_array(3) -= generate_temporary_array(1)) == dummy_array);
        REQUIRE((generate_temporary_array(1) *= generate_temporary_array(2)) == dummy_array);
        REQUIRE((generate_temporary_array(4) /= generate_temporary_array(2)) == dummy_array);
        REQUIRE((generate_temporary_array(1) +  generate_temporary_array(1)) == dummy_array);
        REQUIRE((generate_temporary_array(3) -  generate_temporary_array(1)) == dummy_array);
        REQUIRE((generate_temporary_array(1) *  generate_temporary_array(2)) == dummy_array);
        REQUIRE((generate_temporary_array(4) /  generate_temporary_array(2)) == dummy_array);
    }
    SECTION("Constant arithmetic") {
        REQUIRE((generate_temporary_array(1) += 1) == dummy_array);
        REQUIRE((generate_temporary_array(3) -= 1) == dummy_array);
        REQUIRE((generate_temporary_array(1) *= 2) == dummy_array);
        REQUIRE((generate_temporary_array(4) /= 2) == dummy_array);
        REQUIRE((generate_temporary_array(1) +  1) == dummy_array);
        REQUIRE((generate_temporary_array(3) -  1) == dummy_array);
        REQUIRE((generate_temporary_array(1) *  2) == dummy_array);
        REQUIRE((generate_temporary_array(4) /  2) == dummy_array);
    }
    SECTION("Unary operators") {
        fill_dummy_array(-1);
        REQUIRE((+-temporary_array)==dummy_array);
    }
    SECTION("Operators"){
        Array<N,M,Complex<MyType>> temporary_complex_array;
        Matrix<N,M,Complex<MyType>> temporary_complex_matrix;
        temporary_complex_array.gaussian_random();
        temporary_complex_matrix.gaussian_random();
        std::cout << temporary_complex_matrix.real();
        Array<N,M,Complex<MyType>> temp = temporary_complex_array.conj()+temporary_complex_array;
        std::cout << temp.real();
    }
 }