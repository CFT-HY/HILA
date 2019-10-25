
#include <iostream>
#include <string>
#include <math.h>


#include "../plumbing/field.h"
#include "../datatypes/matrix.h"

//#include "sub.h"

int main() 
{
  field<cmplx<double>> df;
  cmplx<double> cd;
  cmplx<int> ci;
  
  
  cd = +cd;
  ci = cmplx<double>(0,1);
  df = 2 + 3.0_i;
  
  df[ALL] = df[X] + df[X];
    
  df = df + df;
  
  
  return 0;
}

