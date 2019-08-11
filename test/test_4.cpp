
#include "field.h"

template <typename T>
field<T> sum(field<T> a, field<T> b) {
  field<T> r;
  r[ALL] = a[X] + b[X];
  return r;
}


// extern field<int> glob;

int main()
{
  field<double> x;
  field<int> y;
  
  x = sum(x,x);
  // y = sum(y,y);

  x = x+x;
  

  // is.c = 2;

  
  // i[ALL] = i[X]+j[X];
  
  
  return 0;
}

