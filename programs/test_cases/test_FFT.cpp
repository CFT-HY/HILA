#include "test.h"

int main(int argc, char **argv){

    test_setup(argc, argv);
    
    field<cmplx<double>> f;
    f[ALL] = 1;

    f.FFT();

    finishrun();
}
