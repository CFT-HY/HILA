#include "test.h"

int main(int argc, char **argv){

  using T = matrix<1,1,cmplx<double>>;

  test_setup(argc, argv);
    
  field<T> f, p;
  f[ALL] = 1;

  FFT_field(f, p);

  for(int Index=0; Index<lattice->local_volume(); Index++){
    coordinate_vector c = lattice->site_coordinates(Index);
    int csum = 0;
    foralldir(dir){
      csum += c[dir];
    }
    T elem = p.get_value_at(Index);
    if(csum == 0){
      assert(elem.c[0][0].re == lattice->volume() && "first fft");
      assert(elem.c[0][0].im == 0 && "first fft");
    } else {
      assert(elem.c[0][0].re == 0 && "first fft");
      assert(elem.c[0][0].im == 0 && "first fft");
    }
  }

  FFT_field(p, f);

  for(int Index=0; Index<lattice->local_volume(); Index++){
    coordinate_vector c = lattice->site_coordinates(Index);
    int csum = 0;
    foralldir(dir){
      csum += c[dir];
    }
    T elem = f.get_value_at(Index);
    assert(elem.c[0][0].re == lattice->volume() && "second fft");
    assert(elem.c[0][0].im == 0 && "first fft");
    
  }

  // Test reading and writing a field
  onsites(ALL){
    f[X].random();
  }

  f.write_to_file("test_config_filename");
  p.read_from_file("test_config_filename");

  double sum=0;
  onsites(ALL) {
    sum += (f[X]-p[X]).norm_sq();
  }

  assert(sum==0 && "Write and read field");


  finishrun();
}
