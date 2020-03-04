
#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../plumbing/field.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
cmplx<double> f(cmplx<double> x) {return e(x);}

using ft = cmplx<double>;


int main()
{
  
  field<cmplx<double>> a,b;
  field<double> t(1.0);

  coordinate_vector v = XUP - 2*YUP;
  coordinate_vector y;
  y = v;

  parity p = ODD;
  
  onsites(ALL) {
    ft d(2,2);
    coordinate_vector dd = XUP + YUP;
    element<ft> t;
    a[X] = a[X+XUP];
    b[X] = t;
  }
  
  return 0;
}

