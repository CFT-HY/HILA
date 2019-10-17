#include "bench.h"

int main(){
    bench_setup();
    int n_runs = 100;
    double msecs;
    clock_t init, end;

    field<matrix<1,3,double> > vector1;

    onsites(ALL){
        vector1[X].fill(1);
    }

    double sum=0;
    init = clock();

    for( int i=0; i<n_runs; i++){
        sum += vector1[X].sq_sum();
    }

    end = clock();

    float timing = ((float)(end - init)) *1000 / (CLOCKS_PER_SEC) / n_runs;
    output0 << "Matrix * Matrix: " << timing << " ms \n";

    return 0;
}
