
#include "../plumbing/defs.h"
#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
#pragma transformer ast dump
cmplx<double> f(const cmplx<double> & x) { return e(x);}

class tmp {
  int i,j;
};


using ft = cmplx<double>;

template <typename T>
field<T> ret(field<T> & a) {
  field<T> b;
  b[ALL] = a[X];
  return b;
}


int main()
{
  
  field<cmplx<double>> a,b,c;
#pragma transformer ast dump
  int i;
  field<double> t(1.0);

  coordinate_vector 
  v = XUP - 2*YUP;
  
  parity p = ODD;

  a = b.shift(v);

  #pragma transformer ast dump
  onsites(p) {

    #pragma transformer ast dump
    a[X] = f(b[X]);
    c[X] = a[X];
   
    //    c[X] = a[X];
  }
  
  return 0;
}

