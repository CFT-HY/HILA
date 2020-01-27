#include<iostream>
#include "cmplx.h"
#include "general_matrix.h"
#include "sun.h"
//simple tests to find bugs and errors in general_matrix operations
//(could also implement unit testing on each function and compare
//with some known results)

int main(){
    seed_mersenne(12);
    SU<20, double> a;
    for (int i = 3; i < 30; i++){
        a.random(i);
        std::cout << i << '\t' << a.det_lu().abs() << '\n';
    }
    return 0;
}
