
#include "../plumbing/field.h"


template <typename T>
class c {
public:
  template <typename U>
  class d {
  public:
    static U sum(U & a, T & b) { 
      T r;
      r[ALL] = a[X] + b[X]; 
      return r;
    }
  };
  
  template <typename U>
  U sumt(T& a, U& b) {
    T r;
    r[ALL] = a[X] + b[X]; 
    return r;
  }
  
  class f {
  public:
    inline static field<double> sum2(field<double> & a, field<double> & b) {
      field<double> r;
      r[ALL] = a[X] + b[X];
      return r;
    }
  };
};


transformer_ctl(dump_ast);
template <typename T>
inline T sum(const T a, const T b) {
  T r;
  r[ALL] = a[X] + b[X];

  return r;
}

// template <>  
// field<double> sum<field<double>>(const field<double> a, const field<double> b) {
//   field<double> r;
//   r[ALL] = a[X] + b[X];
// 
//   return r;
// }


// extern field<int> glob;

int main()
{
  field<double> a(2.0),x;
  field<int> y;
  c<field<int>>::f fv;
  c<field<double>> cv;
  double dd;

//  int j = c<int>::sum( 1, 2);
  
  x = cv.sumt(a,x);
  x = fv.sum2(a,x);
  
  extern int kissa;
  x = c<field<double>>::d<field<double>>::sum( a, x);
  
  x = sum(a,x);
  y = sum(y,y);
  

  x = x+x+2.0;
  

  // is.c = 2;

  
  x[ALL] = x[X]+x[X];
  
  onsites(EVEN) {
    for (int k=0; k<NDIM; k++) {
      x[X] = a[X+direction::xup] + dd*x[X];
    }
  }
  
  
  return 0;
}

