
#include "plumbing/defs.h"
#include "datatypes/cmplx.h"

#include "plumbing/field.h"

// extern field<int> glob;

Cmplx<double> d(Cmplx<double> x) {return x;}
Cmplx<double> e(Cmplx<double> x) {return d(x);}
// #pragma hila ast dump
Cmplx<double> f(const Cmplx<double> & x) { return e(x);}



template <typename T>
T xyz( output_only T & v) {
  v = 1.3;
  return sin(v);
}


using ft = Cmplx<double>;

template <typename T>
class v2 {
public:
  using base_type = typename base_type_struct<T>::type;
  Cmplx<T> a[2];
  
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
  
  field<Cmplx<double>> a,b,c;
  int i;
  field<double> t(1.0),s;
  field<Cmplx<float>> kissa;

  auto y = X;
  
  // pf(t);
  
  field<v2<double>> A;  

  
  coordinate_vector 
  v = e_x - 2*e_y;
    
  parity p = ODD;

  // a[ALL] = b[X+2*e_x+e_y];
  
  a = b.shift(v);

  direction d = e_x, d2 = e_y;
  
  A[ALL] = { Cmplx(1,0), Cmplx(0,0) };
  
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

