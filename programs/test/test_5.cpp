
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
  
  // transformer_ctl(dump_ast);
  // x[EVEN] = f(a[X]) + x[X];

  
  
  a[0][ALL] = a[1][X];
  
  
  return 0;
}

