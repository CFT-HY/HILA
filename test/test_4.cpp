
#include "field.h"



template <int n, typename T>
struct s {
  T c[n];
};

using vec = s<3,double>;


template <>
class field<int>  {
  int j;
};


// extern field<int> glob;

int main()
{
  field<double> lf;
  field<double> dd;
  field<double> t;
  field<int> k;
  field<vec> sa;
  // struct s<int> is;
  
  lf[ALL] = dd[X]+t[X];

  // is.c = 2;

  
  // i[ALL] = i[X]+j[X];
  
  
  return 0;
}

