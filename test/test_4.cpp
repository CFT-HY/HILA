
#include "field.h"

template <typename T>
class c {
public:
  static T sum(T a, T b) { return a + b; }
};


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
  field<double> a,x;
  field<int> y;

#pragma xyz kissa
  x = c<field<double>>::sum( a, x);
  
  double z = c<double>::sum( 2.0, 3.0);
  
  x = sum(a,x);
  
  x = x+x+2.0;
  

  // is.c = 2;

  
  x[ALL] = x[X]+x[X];
  
  
  return 0;
}

