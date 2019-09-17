
#include <iostream>
#include <string>
#include <math.h>

#include "../plumbing/field.h"
#include "../datatypes/matrix.h"

//#include "sub.h"

int main() 
{
  field<cmplx<double>> df;
  
  df = 1;
  
  df[ALL] = df[X] + df[X];
    
  df = df + df;
  
  
  return 0;
}

