#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>

class lattice_struct {
  std::array<unsigned,NDIM> size;
  
};

extern lattice_struct * current_lattice;

#endif
