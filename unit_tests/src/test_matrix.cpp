#include "hila.h"
#include "catch.hpp"

#define N 10
#define M 10
using MyType = double;

class MatrixTest {

  public:
    Matrix<N, M, MyType> dummy_matrix;
    Vector<N, MyType> dummy_vector;

    template <typename T>
    void fill_dummy_matrix(T assign_value) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                dummy_matrix.c[i * M + j] = assign_value;
            }
        }
    };

    template <typename T>
    void diag_dummy_matrix(T assign_value) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (i==j) dummy_matrix.c[i * M + j] = assign_value;
                else dummy_matrix.c[i * M + j] = 0;
            }
        }
    };

    template <typename T>
    void fill_dummy_vector(T assign_value) {
        for (int i = 0; i < N; i++) {
            dummy_vector.c[i] = assign_value;
        }
    };

    template <typename T>
    Matrix<N, M, MyType> generate_temporary_matrix(T assign_value) {
        Matrix<N, M, MyType> temp(assign_value);
        return temp;
    }
};

TEST_CASE_METHOD(MatrixTest, "Matrix constructor", "[Matrix]") {
    SECTION("Constructor from CoordinateVector") {
        dummy_matrix.random();
        auto temporary_matrix(dummy_matrix);
        REQUIRE(temporary_matrix == dummy_matrix);
    }
    SECTION("Construcor from constant") {
        diag_dummy_matrix(1);
        Matrix<N, M, MyType> temporary_matrix(1);
        REQUIRE(temporary_matrix == dummy_matrix);
    }
    SECTION("Constructor from zero") {
        diag_dummy_matrix(0);
        Matrix<N, M, MyType> temporary_matrix(0);
        REQUIRE(temporary_matrix == dummy_matrix);
    }
    SECTION("Constructor from initializer list") {
        fill_dummy_matrix(1);
        Matrix<N, M, MyType> temporary_matrix = {
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        REQUIRE(temporary_matrix == dummy_matrix);
    }
}

TEST_CASE_METHOD(MatrixTest, "Matrix assignment", "[Matrix]") {
    Matrix<N, M, MyType> temporary_matrix;
    SECTION("Assignment from matrix") {
        dummy_matrix.random();
        temporary_matrix = dummy_matrix;
        REQUIRE(temporary_matrix == dummy_matrix);
    }
    SECTION("Assignment from scalar") {
        diag_dummy_matrix(1);
        temporary_matrix = 1;
        REQUIRE(temporary_matrix == dummy_matrix);
    }
    SECTION("Assignment from initializer list") {
        fill_dummy_matrix(1);
        temporary_matrix = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        REQUIRE(temporary_matrix == dummy_matrix);
    }
}

TEST_CASE_METHOD(MatrixTest, "Indexing", "[Matrix]") {
    SECTION("Matrix indexing") {
        Matrix<N, M, MyType> temporary_matrix(0);
        fill_dummy_matrix(1);
        int index_i, index_j;
        temporary_matrix.c[2 * M + 1] = 1;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (temporary_matrix.e(i, j) == dummy_matrix.e(i, j)) {
                    index_i = i;
                    index_j = j;
                }
            }
        }
        REQUIRE(index_i == 2);
        REQUIRE(index_j == 1);
    }
    SECTION("Vector indexing") {
        Vector<N, MyType> temporary_vector(0);
        int index_e, index_normal;
        fill_dummy_vector(1);
        temporary_vector.c[3] = 1;
        for (int i = 0; i < N; i++) {
            if (temporary_vector.e(i) == dummy_vector.e(i)) index_e = i;
            if (temporary_vector[i] == dummy_vector[i]) index_normal = i;
        }
        REQUIRE(index_e == 3);
        REQUIRE(index_normal == 3);
    }
}

TEST_CASE_METHOD(MatrixTest, "Matrix mathematical operations", "[Matrix]") {
    SECTION("Matrix arithmetic") {
        diag_dummy_matrix(2);
        REQUIRE((generate_temporary_matrix(1) += generate_temporary_matrix(1)) ==
                dummy_matrix);
        REQUIRE((generate_temporary_matrix(3) -= generate_temporary_matrix(1)) ==
                dummy_matrix);
        REQUIRE((generate_temporary_matrix(1) *= generate_temporary_matrix(2)) ==
                dummy_matrix);
        REQUIRE((generate_temporary_matrix(1) + generate_temporary_matrix(1)) ==
                dummy_matrix);
        REQUIRE((generate_temporary_matrix(3) - generate_temporary_matrix(1)) ==
                dummy_matrix);
        REQUIRE((generate_temporary_matrix(1) * generate_temporary_matrix(2)) ==
                dummy_matrix);
    }
    SECTION("Constant arithmetic") {
        diag_dummy_matrix(2);
        REQUIRE((generate_temporary_matrix(1) += 1) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(3) -= 1) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(1) *= 2) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(4) /= 2) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(1) + 1) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(3) - 1) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(1) * 2) == dummy_matrix);
        REQUIRE((generate_temporary_matrix(4) / 2) == dummy_matrix);
    }
    SECTION("Unary operators") {
        Matrix<N, M, MyType> temporary_matrix(1);
        dummy_matrix = -1;
        REQUIRE((+-temporary_matrix) == dummy_matrix);
    }
    SECTION("Complex operators") {
        INFO("Testing Matrix:: imag, real, conj functions")
        fill_dummy_matrix(0);
        Matrix<N, M, Complex<MyType>> temporary_complex_matrix;
        temporary_complex_matrix.gaussian_random();
        GIVEN("A complex matrix") {
            WHEN("Conjugate is take") {
                THEN("Sum with original will have 0 imaginary part") {
                    Matrix<N, M, Complex<MyType>> temp =
                        temporary_complex_matrix.conj() + temporary_complex_matrix;
                    REQUIRE(temp.imag() == dummy_matrix);
                }
                THEN("Subracted with original will have 0 real part") {
                    Matrix<N, M, Complex<MyType>> temp =
                        temporary_complex_matrix.conj() - temporary_complex_matrix;
                    REQUIRE(temp.real() == dummy_matrix);
                }
            }
        }
    }
    SECTION("Squarenorm") {
        fill_dummy_matrix(1);
        REQUIRE((dummy_matrix.squarenorm() / (N * M)) == 1);
    }
}