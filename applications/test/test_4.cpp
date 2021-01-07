
#include "plumbing/field.h"


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
    inline static Field<double> sum2(Field<double> & a, Field<double> & b) {
      Field<double> r;
      r[ALL] = a[X] + b[X];
      return r;
    }
  };
  
};


// test also dependent functions

template <typename V>
class h : c<V> {
public:
  V var;
  V func(V & a, V & b) { return this->sumt(a,b); }
};


template <typename T>
inline T sum(const T a, const T b) {
  T r;
  r[ALL] = a[X] + b[X];

  return r;
}

// template <>  
// Field<double> sum<Field<double>>(const Field<double> a, const Field<double> b) {
//   Field<double> r;
//   r[ALL] = a[X] + b[X];
// 
//   return r;
// }


// extern Field<int> glob;

int main()
{
  Field<double> a(2.0),x;
  Field<int> y;
  c<Field<int>>::f fv;
  c<Field<double>> cv;
  double dd;
  h<Field<double>> hv;

//  int j = c<int>::sum( 1, 2);
  
  // x = cv.sumt(a,x);
  x = fv.sum2(a,x);
  
  x = hv.func(a,x);
  
  extern int kissa;
  // x = c<Field<double>>::d<Field<double>>::sum( a, x);
  
  x = sum(a,x);
  y = sum(y,y);
  

  x = x+x+2.0;
  

  // is.c = 2;

  
  x[ALL] = x[X]+x[X];
  
  onsites(EVEN) {
    for (int k=0; k<NDIM; k++) {
      x[X] = a[X+e_x] + dd*x[X];
    }
  }
  
  
  return 0;
}

