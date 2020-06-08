
#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"

#include "../plumbing/field.h"

// extern field<int> glob;

cmplx<double> d(cmplx<double> x) {return x;}
cmplx<double> e(cmplx<double> x) {return d(x);}
// #pragma hila ast dump
cmplx<double> f(const cmplx<double> & x) { return e(x);}



template <typename T>
T xyz( output_only T & v) {
  v = 1.3;
  return sin(v);
}


using ft = cmplx<double>;

template <typename T>
class v2 {
public:
  using base_type = typename base_type_struct<T>::type;
  cmplx<T> a[2];
  
  void setter() output_only  { a[0]=1; }
};


//template <typename g>
//field<g> f2(const field<g> &f, int i);


double dv( double * d ) {
  return *d + 1;
}



template <typename g>
void pf(field<g> & s) {
  s[ALL] = 1;
}

int main()
{
  
  field<cmplx<double>> a,b,c;
  int i;
  field<double> t(1.0),s;
  field<cmplx<float>> kissa;

  auto y = X;
  
  // pf(t);
  
  field<v2<double>> A;  

  
  coordinate_vector 
  v = XUP - 2*YUP;
    
  parity p = ODD;

  // a[ALL] = b[X+2*XUP+YUP];
  
  a = b.shift(v);

  direction d = XUP, d2 = YUP;
  
  A[ALL] = { cmplx(1,0), cmplx(0,0) };
  
  double dvar;
  onsites(p) {

    kissa[X] = 2;
    A[X].setter();
    auto tv = t[X];
    s[X] = xyz(tv);
    // t[X] = tmp;
    
  }
  
  return 0;
}

