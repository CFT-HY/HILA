#include "bench.h"

int main(){
    bench_setup();
    int n_runs = 100;
    double msecs;
    clock_t init, end;

    field<matrix<3,3,double> > matrix1;
    field<matrix<1,3,double> > vector1;
    field<matrix<1,3,double> > vector2;

    matrix1[ALL] = 1; 
    onsites(ALL){
        vector1[X].c[0][0]=1;
        vector1[X].c[0][1]=1;
        vector1[X].c[0][2]=1;
    }

    init = clock();

    for( int i=0; i<n_runs; i++){
        vector2[ALL] = vector1[X]*matrix1[X];
    }

    end = clock();

    float timing = ((float)(end - init)) *1000 / (CLOCKS_PER_SEC) / n_runs;
    output0 << "Matrix * Matrix: " << timing << " ms \n";

    return 0;
}
