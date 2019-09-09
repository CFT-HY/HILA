#ifndef LATTICE_H
#define LATTICE_H

#include <iostream>
#include <fstream>

// TODO: assertion moved somewhere where basic params
#undef NDEBUG
#include <assert.h>

void * allocate_aligned


// move these somewhere - use consts?
#define NDIM 4
#define NDIRS (2*NDIM)


class lattice_struct {
private:
  unsigned size_p[NDIM];
  
  unsigned sites_on_node_p;
  unsigned node_fullsize_p;
  
public:
  unsigned size(int d) { return size_p[d]; }
  unsigned volume() { 
    unsigned v = size_p[0];
    for (int i=1; i<NDIM; i++) v *= size_p[i];
    return v;
  }

  unsigned sites_on_node() { return sites_on_node_p; }
  unsigned node_fullsize() { return node_fullsize_p; }
};

extern lattice_struct * current_lattice;

#endif
