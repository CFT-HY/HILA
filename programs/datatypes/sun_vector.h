#ifndef SU_VEC
#define SU_VEC

#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"


static double gaussian_ran() 
{
  static double VARIANCE = 0.5;
  static double second;
  static bool draw_new = true;
  if( draw_new ) {
    double first, phi, urnd, r;
    phi = 2.0 * 3.141592654 * (double) hila_random();
    urnd = (double)(1.0-hila_random());
    r  = sqrt( -log(urnd) * (2.0 * VARIANCE) );
    first = (double)(r*sin(phi));
    second = (r*cos(phi));

    draw_new = false;
    return first;
  } else {
    draw_new = true;
    return second;
  }
}




template<int n, typename radix>
class SU_vector : public matrix<1,n,cmplx<radix> >{
  public:
  using base_type = typename base_type_struct<radix>::type;
  
  SU_vector() = default;

  SU_vector(matrix<1,n,cmplx<radix>> m) {
    for (int i=0; i<n; i++){
      this->c[i][0] = m.c[i][0];
    }
  }

  void gaussian(){ 
    for (int i = 0; i < n; i++) {
      (*this).c[i][0].re = gaussian_ran();
      (*this).c[i][0].im = gaussian_ran();
    }
  }

  

};

#endif