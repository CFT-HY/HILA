//////////////////////////////////////////////////////////////////////////////
/// Test the standard gather here
//////////////////////////////////////////////////////////////////////////////

#include "defs.h"
#include "coordinates.h"
#include "lattice.h"
#include "field.h"


template <typename T>
// #pragma transformer dump ast
struct test_struct {
  T r[NDIM];

  using base_type = typename base_type_struct<T>::type;
};

using test_int = test_struct<int>;
using test_double = test_struct<double>;

void test_std_gathers()
{

  extern lattice_struct * lattice;
  field<test_int> t;
  field<double> f;
  
  onsites(ALL) {
    coordinate_vector v = X.coordinates();
    foralldir(d) {
      t[X].r[d] = v[d];
    }
    // foralldir(d)
    //   hila::output << t[X].r[d] << ' ';
    // hila::output << '\n';
  }

  for (parity p : {EVEN,ODD,ALL}) {

    foralldir(d) {
      direction d2;
      for (d2=d; is_up_dir(d2); d2=-d) {
      
        int diff = 0;
        int add;
        if (is_up_dir(d2)) add = 1; else add = -1;
        onsites(p) {
                  
          int j = t[X+d2].r[d];
//        if (j==0) j = 0;
          int s = (t[X].r[d] + add + lattice->size(d)) % lattice->size(d);

          int lv = s-j;
          int a = 0;
          foralldir(dir) if (dir != d) a+= t[X+d2].r[dir] - t[X].r[dir];
          
          if (lv != 0 || a != 0) {
            hila::output << "diff != 0! at " << X.coordinates() << " direction " << d2 
                         << " parity " << (int)p << '\n';
            hila::output << "t[X+d2].r[d] = " << j << " should be " << s << " a is " << a << '\n';
            for (int loop=0; loop<NDIM; loop++) hila::output << t[X+d2].r[loop] << ' ';
            hila::output << " - ";
            for (int loop=0; loop<NDIM; loop++) hila::output << t[X].r[loop] << ' ';
            
            hila::output << '\n';

 
            // exit(-1);
          }
          
          int i = s - j;

#if (0 && !defined(VANILLA) && !defined(TRANSFORMER))
          for (int k=0; k<8; k++) if (i[k] != 0) {
            hila::output << "Error!  node " << mynode() << " parity " 
                         << parity_name(p) << " direction " << (unsigned)d2 << '\n';
            hila::output << "\nCoordinate was ";
            for (int l=0; l<8; l++) hila::output << j[l] << ' ';
            hila::output << "\nShould be ";
            for (int l=0; l<8; l++) hila::output << s[l] << ' ';
            hila::output << '\n';
            exit(-1);
          }
#endif
        }

        t.mark_changed(ALL);  // foorce fetching

#ifdef VECTORIZED
        // above is not vectorized, so do it also in vec way

        diff = 0;
        onsites(p) {
          int j = t[X+d2].r[d];
          int s = (t[X].r[d] + add + lattice->size(d)) % lattice->size(d);

          diff += s-j;
        }     
     
        if (diff != 0) {
          hila::output << "Std gather test error! Node " << mynode() 
                       << " Parity " << parity_name(p) << " direction " << (unsigned)d2 << '\n';
          exit(-1);
        }

        t.mark_changed(ALL);
#endif

      }
    }
  }
}

