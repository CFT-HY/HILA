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


template <typename T>
void gather_test() {

  extern lattice_struct * lattice;
  field<test_struct<T>> t;
  
  onsites(ALL) {
    coordinate_vector v = X.coordinates();
    foralldir(d) {
      t[X].r[d] = v[d];
    }
  }

  for (parity p : {EVEN,ODD,ALL}) {

    foralldir(d) {
      direction d2;
      for (d2=d; is_up_dir(d2); d2=-d) {
      
        T diff = 0;
        int add;
        if (is_up_dir(d2)) add = 1; else add = -1;
        onsites(p) {
          T j = t[X+d2].r[d];
//        if (j==0) j = 0;
          T s;
          s = (t[X].r[d] + add + lattice->size(d)) % lattice->size(d);

          T lv = s-j;
          T a = 0;
          foralldir(dir) if (dir != d) a+= t[X+d2].r[dir] - t[X].r[dir];
          
          if (lv != 0 || a != 0) {
            hila::output << "Error in gather test at " << X.coordinates() << " direction " << d2 
                         << " parity " << (int)p << '\n';
            hila::output << "Fetched element t[X+d2].r[d] = " << j << " should be " << s << " perp diff is " << a << '\n';
            hila::output << "This element - neighbour element:  ";
            for (int loop=0; loop<NDIM; loop++) hila::output << t[X].r[loop] << ' ';
            hila::output << " - ";
            for (int loop=0; loop<NDIM; loop++) hila::output << t[X+d2].r[loop] << ' ';
            
            hila::output << '\n';

            exit(-1);
          }
        }

        t.mark_changed(ALL);  // foorce fetching, test it too

#ifdef VECTORIZED
        // above is not vectorized, so do it also in vec way

        
        diff = 0;
        onsites(p) {
          T j = t[X+d2].r[d];
          T s = (t[X].r[d] + add + lattice->size(d)) % lattice->size(d);

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



void test_std_gathers()
{
  gather_test<int>();
  // gather_test<int64_t>();
}

