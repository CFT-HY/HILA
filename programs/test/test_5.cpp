
#include "../plumbing/defs.h"
#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
// #pragma transformer ast dump
cmplx<double> f(const cmplx<double> & x) { return e(x);}



using ft = cmplx<double>;

template <typename T>
class v2 {
public:
  using base_type = typename base_type_struct<T>::type;
  cmplx<T> a[2];
};


#pragma kissa


int main()
{
  
  field<cmplx<double>> a,b,c;
  int i;
  field<double> t(1.0);

  field<v2<double>> A;  
  
  coordinate_vector 
  v = XUP - 2*YUP;
  
  parity p = ODD;

  // a[ALL] = b[X+2*XUP+YUP];
  
  a = b.shift(v);

  direction d = XUP, d2 = YUP;
  
  A[ALL] = { cmplx(1,0), cmplx(0,0) };
  
  onsites(p) {
    double dv = t[X];

    double t2;

    double t4 = 5.0;

    t4 = t2;
    
    t2 = t[X];
    
    double vv = (double)X.parity();

    if (t2 < 4) t[X] += 4;
    
    A[X].a[0] = 0;

    v2<double> vvv = A[X];
    // c[X] = b[X+(XUP+YUP)];
    // #pragma transformer ast dump
    c[X] = a[X];
    
    //    c[X] = a[X];
  }
  
  return 0;
}

