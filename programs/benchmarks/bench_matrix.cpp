#include "bench.h"

int main(){
    bench_setup();
    int n_runs = 100;
    double msecs;
    clock_t init, end;

    field<matrix<2,2,double> > matrix1;
    field<matrix<2,2,double> > matrix2;
    field<matrix<2,2,double> > matrix3;

    matrix1[ALL] = 1; 
    matrix2[ALL] = 1;

    init = clock();

    for( int i=0; i<n_runs; i++){
        matrix3[ALL] = matrix1[X]*matrix2[X];
    }

    end = clock();

    float timing = ((float)(end - init)) *1000 / (CLOCKS_PER_SEC) / n_runs;
    output0 << "Matrix * Matrix: " << timing << " ms \n";

    return 0;
}
