#include "hila.h"
#include "catch.hpp"

#define N 10
#define M 10
using MyType = double;

class MatrixTest {

  public:
    Matrix<N, M, MyType> dummy_matrix = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                    10,11,12,13,14,15,16,17,18,19,
                    20,21,22,23,24,25,26,27,28,29,
                    30,31,32,33,34,35,36,37,38,39,
                    40,41,42,43,44,45,46,47,48,49,
                    50,51,52,53,54,55,56,57,58,59,
                    60,61,62,63,64,65,66,67,68,69,
                    70,71,72,73,74,75,76,77,78,79,
                    80,81,82,83,84,85,86,87,88,89,
                    90,91,92,93,94,95,96,97,98,99};
    Vector<N, MyType> dummy_vector = {0,10,20,30,40,50,60,70,80,90};
    HorizontalVector<M, MyType> dummy_vector_t = {0,1,2,3,4,5,6,7,8,9};

    template <typename T>
    void fill_dummy_matrix(T assign_value) {
        for (int i = 0; i < N * M; i++)
            dummy_matrix.c[i] = assign_value;
    };

    template <typename T>
    void diag_dummy_matrix(T assign_value) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                if (i == j)
                    dummy_matrix.c[i * M + j] = assign_value;
                else
                    dummy_matrix.c[i * M + j] = 0;
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
        temporary_matrix = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        REQUIRE(temporary_matrix == dummy_matrix);
    }
}

TEST_CASE_METHOD(MatrixTest, "Utilities", "[Matrix]") {
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
            if (temporary_vector.e(i) == dummy_vector.e(i))
                index_e = i;
            if (temporary_vector[i] == dummy_vector[i])
                index_normal = i;
        }
        REQUIRE(index_e == 3);
        REQUIRE(index_normal == 3);
    }
    SECTION("Get row and column") {
        REQUIRE(dummy_matrix.row(0) == dummy_vector_t);
        REQUIRE(dummy_vector.column(0) == dummy_vector);
    }
    SECTION("Set row and column") {
        dummy_matrix.set_row(1,dummy_vector_t);
        REQUIRE(dummy_matrix.row(1) == dummy_vector_t);
        dummy_matrix.set_column(1,dummy_vector);
        REQUIRE(dummy_matrix.column(1) == dummy_vector);
    }
    SECTION("Numpy style fill") {
        fill_dummy_matrix(0);
        Matrix<N,M,MyType> temporary_matrix;
        temporary_matrix.fill(0);
        REQUIRE(dummy_matrix == temporary_matrix);
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

    SECTION("Matrix operators") {
        INFO("Testing Matrix:: trace, mul_trace, dot, transpose") 
            
        Matrix<N,M,double> temporary_matrix = { 0,10,20,30,40,50,60,70,80,90,
                                                1,11,21,31,41,51,61,71,81,91,
                                                2,12,22,32,42,52,62,72,82,92,
                                                3,13,23,33,43,53,63,73,83,93,
                                                4,14,24,34,44,54,64,74,84,94,
                                                5,15,25,35,45,55,65,75,85,95,
                                                6,16,26,36,46,56,66,76,86,96,
                                                7,17,27,37,47,57,67,77,87,97,
                                                8,18,28,38,48,58,68,78,88,98,
                                                9,19,29,39,49,59,69,79,89,99};

        REQUIRE(temporary_matrix.transpose() == dummy_matrix);

        double trace_sum = 495;
        double mul_trace_sum = 261525;
        REQUIRE(temporary_matrix.trace() == trace_sum);
        REQUIRE(dummy_matrix.mul_trace(dummy_matrix) == mul_trace_sum);

        dummy_vector = {1,1,1,1,1,1,1,1,1,1};
        dummy_vector_t = {1,1,1,1,1,1,1,1,1,1};
        Vector<N,MyType> temporary_vector = {1,1,1,1,1,1,1,1,1,1};
        REQUIRE(dummy_vector.dot(temporary_vector)==10);
        REQUIRE(dummy_vector_t == temporary_vector.transpose());
    }

    SECTION("Complex operators") {
        INFO("Testing Matrix:: imag, real, conj, dagger and adjoint")
        fill_dummy_matrix(0);
        Matrix<N, M, Complex<MyType>> temporary_complex_matrix;
        temporary_complex_matrix.gaussian_random();
        GIVEN("A complex matrix") {
            WHEN("Conjugate is taken") {
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
            WHEN("Dagger is taken") {
                THEN("Sum with transpose will have 0 imaginary part"){
                    Matrix<N, M, Complex<MyType>> temp =
                        temporary_complex_matrix.dagger() + temporary_complex_matrix.transpose();
                    REQUIRE(temp.imag() == dummy_matrix);
                }
                THEN("Subracted with transpose will have 0 real part") {
                    Matrix<N, M, Complex<MyType>> temp =
                        temporary_complex_matrix.dagger() - temporary_complex_matrix.transpose();
                    REQUIRE(temp.real() == dummy_matrix);
                }
            }
            WHEN("Adjoing is taken") {
                THEN("Produces same result as dagger") {
                    REQUIRE(temporary_complex_matrix.dagger() == temporary_complex_matrix.adjoint());
                }
            }
        }
    }
    
    SECTION("Squarenorm and norm") {
        fill_dummy_matrix(1);
        REQUIRE((dummy_matrix.squarenorm() / (N * M)) == 1);
        REQUIRE(dummy_matrix.norm() == 10);
    }

    SECTION("Determinant") {
        INFO("Det function has 4 overloads with naive versions up to 3,3 sized matrices")
        diag_dummy_matrix(1);

        REQUIRE(det(dummy_matrix)  == 1);
        REQUIRE(det_lu(dummy_matrix)  == 1);

        fill_dummy_matrix(1);
        REQUIRE(det(dummy_matrix) == 0);
        
        Matrix<1,1,MyType> temp_1(1);
        REQUIRE(det(temp_1)  == 1);
        Matrix<2,2,MyType> temp_2(1);
        REQUIRE(det(temp_2)  == 1);
        Matrix<3,3,MyType> temp_3(1);
        REQUIRE(det(temp_3)  == 1);

        temp_1.fill(1);
        REQUIRE(det(temp_1) == 1);
        temp_2.fill(1);
        REQUIRE(det(temp_2) == 0);
        temp_3.fill(1);
        REQUIRE(det(temp_3) == 0);

    }

}