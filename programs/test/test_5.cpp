
#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"


// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
cmplx<double> f(cmplx<double> x) {return e(x);}

int main()
{
  
  field<cmplx<double>> a,x;
  field<double> t(1.0);
  
  // transformer_ctl(dump_ast);
  // x[EVEN] = f(a[X]) + x[X];

  
  
  x = t;
  
  
  return 0;
}

