
#include "../plumbing/defs.h"
#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
cmplx<double> f(const cmplx<double> & x) { return e(x);}

transformer_ctl_dump_ast();
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
  transformer_ctl_dump_ast();
  int i;
  transformer_ctl_dump_ast();
  field<double> t(1.0);

  coordinate_vector v = XUP - 2*YUP;
  
  transformer_ctl_dump_ast();
  parity p = ODD;
  
  a = b.shift(v);
  
  onsites(p) {

    a[X] = f(b[X]);
    c[X] = a[X];
   
    //    c[X] = a[X];
  }
  
  return 0;
}

