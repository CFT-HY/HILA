
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

#include "plumbing/field.h"


template<typename T>
class testclass {
public:
  field<T> a;
  testclass(field<T> _a) : a(_a) {}

  void sub(field<T> b){
    a[ALL] -= b[X];
  }

  void add(field<T> b){
    onsites(ALL){
      a[X] += b[X];
    }
  }
};


int main()
{
  field<double> a,b,c;

  a[ALL]=1; b[ALL]=1;

  testclass<double> tc(a);
  tc.sub(b);
  tc.add(b);
  
  return 0;
}

