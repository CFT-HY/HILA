#include "field.h"

double f(field<double> t);

class cclass {
  double d;
public:
  double get() {return d;}
};

int main() 
{
  int i;
  field<double> lf;
  field<double> dd;
  double dp[10], somevar, b, c;
  cclass tst;
  
  // cout << "Starting..\n";
  
  //a[EVEN];

  parity p=EVEN;
  lf[EVEN] = dd[X+direction::xup] + exp(dp[1] + somevar +  dp[2] + somevar );

  onsites(p)  {
    for (int k=0; k<10; k++) {
      lf[X] = dd[X+direction::xup] + a + b + c;
    }
    auto t = lf[X];
  }
  
  std::cout << "After 1st\n";

 
  
  return 0;
}

