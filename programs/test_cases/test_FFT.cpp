#include "test.h"

int main(int argc, char **argv){
  
  test_setup(argc, argv);
    
  field<cmplx<double>> f;
  f[ALL] = 1;

  f.FFT();

  for(int Index=0; Index<lattice->local_volume(); Index++){
    coordinate_vector c = lattice->site_coordinates(Index);
    int csum = 0;
    foralldir(dir){
      csum += c[dir];
    }
    cmplx<double> elem = f.get_value_at(Index);
    if(csum == 0){
      assert(elem.re == lattice->volume());
      assert(elem.im == 0);
    } else {
      assert(elem.re == 0);
      assert(elem.im == 0);
    }
  }

  

  finishrun();
}
