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
  double d[10];
  cclass tst;
  
  // cout << "Starting..\n";
  
  //a[EVEN];

  parity p=EVEN;
  onsites(p)  {
    lf[X] = dd[X+direction::xup] + exp(d[0]+d[1]);
    
  }
  
  std::cout << "After 1st\n";

 
  
  return 0;
}

