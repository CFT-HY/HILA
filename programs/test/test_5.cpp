
#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../plumbing/field.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
cmplx<double> f(cmplx<double> x) {return e(x);}

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
  field<double> t(1.0);

  coordinate_vector v = XUP - 2*YUP;
  
  parity p = ODD;

  a = b.shift(v);
  
  onsites(p) {

    a[X] = a[X+opp_dir(XUP)];
    b[X] = a[X];
    c[X] = a[X];
   
    //    c[X] = a[X];
  }
  
  return 0;
}

