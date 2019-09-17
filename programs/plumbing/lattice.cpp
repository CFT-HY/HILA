
#include "lattice.h"

void lattice_struct::setup(int siz[NDIM]) {
  volume = 1;
  for (int i=0; i<NDIM; i++) {
    size[i] = siz[i];
    volume *= siz[i];
  }
    
}


#if NDIM==4
void lattice_struct::setup(int nx, int ny, int nz, int nt) {
  int s[NDIM] = {nx, ny, nz, nt};
  setup(s);
}
#elif NDIM==3
void lattice_struct::setup(int nx, int ny, int nz) {
  int s[NDIM] = {nx, ny, nz};
  setup(s);
}
#elif NDIM==2
void lattice_struct::setup(int nx, int ny) {
  int s[NDIM] = {nx, ny};
  setup(s);
}
#elif NDIM==1
void lattice_struct::setup(int nx) {
  int s[NDIM] = {nx};
  setup(s);
}



#endif

