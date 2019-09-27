#define NDIM 2

#include <iostream>
#include <string>
#include <math.h>

#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"


int main() 
{
  lattice_struct * lattice;
  lattice->setup( 8, 8 );

  field<cmplx<double>> spin;
  
  spin[ALL] = 1;
  
  return 0;
}
