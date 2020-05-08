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
      this->c[0][i] = m.c[0][i];
    }
  }

  void gaussian(){ 
    for (int i = 0; i < n; i++) {
      (*this).c[0][i].re = gaussian_ran();
      (*this).c[0][i].im = gaussian_ran();
    }
  }


  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
  #pragma transformer loop_function
  SU_vector & operator= (const scalart rhs) {
    for (int i=0; i<n; i++){
      this->c[i][i] = rhs;
    }
    return *this;
  }

};

#endif