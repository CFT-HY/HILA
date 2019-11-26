
#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../plumbing/field.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
cmplx<double> f(cmplx<double> x) {return e(x);}

int main()
{
  
  field<cmplx<double>> a[2];
  field<double> t(1.0);

  onsites(ALL) {
    a[0][X] += 2.0;
  }
  
  return 0;
}

