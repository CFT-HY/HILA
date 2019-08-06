#include "field.h"


void add(field<int> &a, field<int> &b, parity par)
{
  a[par] = b[X + direction::xup];
  
}

//#include "sub.h"

template <typename T>
void sub(field<T> &a, field<T> &b, parity p)
{
  a[p] -= b[X];
}

void f(field<double> &d, double verylongvariablename) {
  d[ALL] = verylongvariablename;
}



class A {
private:
  field<double> d;
public:
  A(double var) {
    d[ALL] = var  ;
  }

  A& add(double var) {
    d[ALL] = var + d[X];         
    return *this;
  }
};



int main() 
{
  int i;
  field<double> lf, tf[2], *p;
  field<int> intfield, if2;
  field<double> dd = 2.0;
  double da[N], t,x, *dp;
  parity par(ODD);
  direction dir;
  A a(0);

  a.add(1);

  
  std::cout << "Starting..\n";
  
  for (i=0; i<N; i++) lf[i] = i;

  dp = &t;
  p = &tf[0];
  (*p) = lf;

  //a[EVEN];
  
  lf[EVEN] = tf[0][X+dir] + sin(da[0]) + (*p)[X];


  std::cout << "After 1st\n";

  if2 = if2 + intfield;
  lf = tf[0] + tf[1];
  lf[ALL] = tf[0][X] + tf[1][X];
  
  lf = tf[0] + (*dp) + 3*5.6;

  
  tf[1][par] = t + 1 + tf[0][parity::x] + 2*t;

  lf[parity::even] = t*tf[1][parity::x] + tf[0][X] * dd[X];

  lf[EVEN] = tf[1][X] + tf[0][X] * dd[X];

  lf[ALL] = tf[1][X] + tf[0][X] * dd[X];

  lf = tf[1] + tf[0] * dd;

  lf[ALL] = exp(lf[X]);
  
  
  
  sub( lf, tf[0], EVEN );

  // sub( intfield, if2, ALL );

  add( intfield, if2, ALL );

  add( intfield, if2, ALL );
  

  // for (i=0; i<N; i++) std::cout << tf[i] << " ";


  
  return 0;
}

