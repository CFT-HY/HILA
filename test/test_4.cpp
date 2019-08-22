
#include "field.h"

template <typename T>
field<T> sum(const field<T> a, const field<T> b) {
  field<T> r;
  r[ALL] = a[X] + b[X];

  return r;
}

template <>  
field<double> sum<double>(const field<double> a, const field<double> b) {
  field<double> r;
  r[ALL] = a[X] + b[X];

  return r;
}



// extern field<int> glob;

int main()
{
  field<double> a,x;
  field<int> y;
  
  x = sum(a,x);
  y = sum(y,y);
  
  x = x+x+2.0;
  

  // is.c = 2;

  
  x[ALL] = x[X]+x[X];
  
  
  return 0;
}

