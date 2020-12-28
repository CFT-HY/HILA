#include "test.h"

int main(int argc, char **argv){

  using T = Matrix<2,2,cmplx<double>>;
  // using T = cmplx<double>;

  test_setup(argc, argv);
    
  field<T> f, f2, p, p2;
  double sum = 0;

  for (int iter=0; iter<3; iter++) {

  // Start with unit field
  f = 1;
  
  // After one FFT the field is 0 except at coord 0
  p2 = 0;
  T m = lattice->volume();
  coordinate_vector c;
  foralldir(d){
    c[d] = 0;
  }
  p2.set_element(m, c);


  FFT_field(f, p);
  
  onsites(ALL) {
    sum += (p[X]-p2[X]).norm_sq();
  }
  output0 << "Sum " << sum << '\n';
  assert(fabs(sum) < 1e-10  && "First FFT\n");
  

  // After two applications the field should be back to a constant * volume
  f2[ALL] = lattice->volume();
  
  FFT_field(p, f, fft_direction::back);
  
  sum = 0;
  double tnorm = 0;
  onsites(ALL) {
    sum += (f[X]-f2[X]).norm_sq();
    tnorm += f[X].norm_sq();
  }
  output0 << "Norm " << sum/tnorm << '\n';
  assert(fabs(sum/tnorm) < 1e-10 && "Second FFT\n");

  }

  // Test reading and writing a field
  onsites(ALL){
    f[X].random();
  }

  // write_fields("test_config_filename", p, f);
  // read_fields("test_config_filename", p, f2);

  //sum=0;
  //onsites(ALL) {
  //  sum += (f2[X]-f[X]).norm_sq();
  //}

  //assert(sum==0 && "Write and read field");


  hila::finishrun();
}
