#include <iostream>
#include <string>
#include <math.h>

#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"


int main() 
{
  field<cmplx<double>> spin;
  
  spin[ALL] = 1;
  
  return 0;
}
