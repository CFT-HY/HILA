
#include "field.h"

template <typename T, typename A>
field<T> sum(const field<T> a, const field<A> b) {
  field<T> r;
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
  
  //x = x+x;
  

  // is.c = 2;

  
  x[ALL] = x[X]+x[X];
  
  
  return 0;
}

