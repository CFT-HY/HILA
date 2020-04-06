#include "test.h"

int main(int argc, char **argv){

  test_setup(argc, argv);
    
  field<cmplx<double>> f, p;
  f[ALL] = 1;

  FFT_field(f, p);

  for(int Index=0; Index<lattice->local_volume(); Index++){
    coordinate_vector c = lattice->site_coordinates(Index);
    int csum = 0;
    foralldir(dir){
      csum += c[dir];
    }
    cmplx<double> elem = p.get_value_at(Index);
    if(csum == 0){
      assert(elem.re == lattice->volume() && "first fft");
      assert(elem.im == 0 && "first fft");
    } else {
      assert(elem.re == 0 && "first fft");
      assert(elem.im == 0 && "first fft");
    }
  }

  FFT_field(p, f);

  for(int Index=0; Index<lattice->local_volume(); Index++){
    coordinate_vector c = lattice->site_coordinates(Index);
    int csum = 0;
    foralldir(dir){
      csum += c[dir];
    }
    cmplx<double> elem = f.get_value_at(Index);
    assert(elem.re == lattice->volume() && "second fft");
    assert(elem.im == 0 && "first fft");
    
  }

  finishrun();
}
