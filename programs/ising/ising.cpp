#define NDIM 2

#include <iostream>
#include <string>
#include <math.h>

#include "../plumbing/field.h"
#include "../datatypes/cmplx.h"

std::ostream &hila::output = std::cout;

lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;

int main() 
{
  lattice->setup( 8, 8 );

  field<cmplx<double>> spin;
  
  spin[ALL] = 1;

  cmplx<double> sum=0;
  onsites(ALL){
    sum += spin[X];
  }

  printf(" %f\n", sum.re);
  printf(" %f\n", spin[12].re);
  
  return 0;
}
