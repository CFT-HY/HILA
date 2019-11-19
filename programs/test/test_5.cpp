
#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../plumbing/field.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
cmplx<double> f(cmplx<double> x) {return e(x);}

int main()
{
  
  field<cmplx<double>> a[10];
  field<double> t(1.0);
  
  // transformer_ctl(dump_ast);
  // x[EVEN] = f(a[X]) + x[X];

  onsites(ALL) {
    for (int i=0; i<2; i++) {
      a[i][X] = a[1-i][X];
    }
  }
  
  return 0;
}

