
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

#include "plumbing/field.h"


template<typename T>
class testclass {
public:
  Field<T> a;
  testclass(Field<T> _a) : a(_a) {}

  void sub(Field<T> b){
    a[ALL] -= b[X];
  }

  void add(Field<T> b){
    onsites(ALL){
      a[X] += b[X];
    }
  }
};


int main()
{
  Field<double> a,b,c;

  a[ALL]=1; b[ALL]=1;

  testclass<double> tc(a);
  tc.sub(b);
  tc.add(b);
  
  return 0;
}

