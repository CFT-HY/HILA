
#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"


// extern field<int> glob;

double d(double x) {return x;}

int main()
{
  
  field<double> a,x;
  field<double> t(1.0);
  
  transformer_ctl(dump_ast);
  x[EVEN] = a[X] + x[X];

  // x = a + x;
  
  
  return 0;
}

