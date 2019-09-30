#include<iostream>
#include "cmplx.h"
#include "general_matrix.h"

//simple tests to find bugs and errors in general_matrix operations
//(could also implement unit testing on each function and compare
//with some known results)

int main(){

    matrix<6,6,cmplx<float>> A = 1, D = 2;
    matrix<6,6,cmplx<float>> B;
    matrix<6,6,cmplx<float>> C;
    matrix<6,1,cmplx<float>> vector, vector2;
    B.fill(2.0);
    vector.fill(1.0);
    vector2.fill(-1.0);

    C=A*B;
    A+=D;
    A-=D;
    A*=D;
    A*=1;

    A=2;
    A+=B;

    cmplx<float> val = det(A);
    std::cout << "determinant: " << val.re << ' ' << val.im<< '\n';
    A.conjugate();
    std::cout << "trace: " << A.trace().re << '\n';

   //linear transformation

    vector = A*vector2;

    //dot product
    cmplx<float> result = vector*vector2;
    return 0;
}
