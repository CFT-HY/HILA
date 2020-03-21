
#include "../plumbing/defs.h"
#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
#pragma transformer ast dump
cmplx<double> f(const cmplx<double> & x) { return e(x);}


using ft = cmplx<double>;

template <typename T>
class tmp {
private:
  T in;
public:
  #pragma transformer loop_function
  T ret(const T & a) {
    return a+in;
  }
};


int main()
{
  
  field<cmplx<double>> a,b,c;
  int i;
  field<double> t(1.0);
  tmp<double> x;
  tmp<cmplx<double>> y;
  
  coordinate_vector 
  v = XUP - 2*YUP;
  
  parity p = ODD;

  a = b.shift(v);

  onsites(p) {

    a[X] = y.ret(b[X]);
    t[X] = x.ret(t[X]);
    
    #pragma transformer ast dump
    a[X] = f(b[X]);
    c[X] = a[X];
   
    //    c[X] = a[X];
  }
  
  return 0;
}

