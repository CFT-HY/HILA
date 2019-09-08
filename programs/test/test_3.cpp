
#include <iostream>
#include <string>
#include <math.h>

#include "field.h"

template <template <typename> class D, typename T=double>
struct F {
  D<T> d;
};

template <typename T=double> 
struct vector {
  T data[10];
};

template <class T, class A>
class c {

  void sub2(field<T> &a, const field<double> &b, parity p)
  {
    // transformer_control("dump-ast");
    a[p] -= b[X];
  }
};


//#include "sub.h"

int main() 
{
  field<vector<>> a,b;
  c<double,double> t;
  F<vector> vd;
  F<vector,float> vf;

  
    
  
  t.sub2(a,b,EVEN);
  
  return 0;
}

