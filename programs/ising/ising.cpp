#define NDIM 2

#include <iostream>
#include <string>
#include <math.h>

#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"

lattice_struct * lattice;
std::ostream &hila::output = std::cout;


int main() 
{
  lattice->setup( 8, 8 );

  field<cmplx<double>> spin;
  
  spin[ALL] = 1;
  
  return 0;
}
