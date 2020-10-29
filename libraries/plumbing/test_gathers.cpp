//////////////////////////////////////////////////////////////////////////////
/// Test the standard gather here
//////////////////////////////////////////////////////////////////////////////

#include "defs.h"
#include "coordinates.h"
#include "lattice.h"

// Declare the test class here, element type definitions must be before field.h
template <typename T>
// #pragma hila dump ast
struct test_struct {
  T r[NDIM];

  using base_type = typename base_type_struct<T>::type;

  // define unary - for antiperiodic b.c.
  #pragma hila loop_function
  test_struct<T> operator-() const { test_struct<T> t; foralldir(d) t.r[d] = -r[d]; return t; }

};

#include "field.h"
#include "datatypes/general_matrix.h"




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

#ifdef SPECIAL_BOUNDARY_CONDITIONS
  for (boundary_condition_t bc : {boundary_condition_t::PERIODIC, boundary_condition_t::ANTIPERIODIC}) {
#if NDIM > 3
      t.set_boundary_condition(TUP,bc);
#endif
#endif
  for (parity p : {EVEN,ODD,ALL}) {

    foralldir(d) {
      // Find size here, cannot call the function in a CUDA loop
      int size_d = lattice->size(d);
      int size_t = lattice->size(TUP);
      for (direction d2 : {d,-d}) {
      
        T diff = 0;
        int add;
        double sum1 = 0, sum2 = 0;  // use double to accumulate ints, should be accurate
        if (is_up_dir(d2)) add = 1; else add = -1;
        onsites(p) {
          auto n = t[X+d2];
#if defined(SPECIAL_BOUNDARY_CONDITIONS) && NDIM > 3
          if (bc == boundary_condition_t::ANTIPERIODIC &&
              (( X.coordinates()[TUP] == 0 && d2 == TDOWN) || 
               ( X.coordinates()[TUP] == size_t-1 && d2 ==TUP))) {
            n = -n;
          }
#endif

          T j = n.r[d];
          T t_r = t[X].r[d];
          T s = ((int)(t_r + add + size_d)) % size_d;

          sum2 += n.r[d] - size_d/2.0;
          sum1 += t[X].r[d] - size_d/2.0;

          T lv = s-j;
          T a = 0;
          foralldir(dir) if (dir != d) a+= n.r[dir] - t[X].r[dir];
          
          #ifndef CUDA
          if (lv != 0 || a != 0) {
            hila::output << "Error in gather test at " << X.coordinates() << " direction " << d2 
                         << " parity " << (int)p << '\n';
            hila::output << "Fetched element t[X+d2].r[d] = " << j << " should be " << s << " perp diff is " << a << '\n';
            hila::output << "This element - neighbour element:  ";
            for (int loop=0; loop<NDIM; loop++) hila::output << t[X].r[loop] << ' ';
            hila::output << " - ";
            for (int loop=0; loop<NDIM; loop++) hila::output << n.r[loop] << ' ';
            
            hila::output << '\n';
          }
          #endif
          assert(lv == 0 || a == 0 && "Test gathers");
        }

        double s_result;
        if (p == ALL) 
          s_result = lattice->volume()/2.0;
        else 
          s_result = lattice->volume()/4.0;

        if (sum1 + s_result != 0.0) {
          output0 << "Error in sum reduction!  answer " << sum1 + s_result << " should be 0, parity " << (int)p << ", direction " << (int)d2 << ", bc " << (int)bc << "\n";
          exit(-1);
        }

        if (sum2 + s_result != 0.0) {
          output0 << "Error in neighbour sum2 reduction!  answer " << sum2 + s_result << " should be 0, parity " << (int)p << ", direction " << (int)d2 << ", bc " << (int)bc << "\n";
          exit(-1);
        }


        t.mark_changed(ALL);  // force fetching, test it too

#ifdef VECTORIZED
        // above is not vectorized, so do it also in vec way

        
        diff = 0;
        sum1 = sum2 = 0;
        onsites(p) {
          int a = t[X+d2].r[d];
          T j = abs(a);
          T s = ((int)(t[X].r[d] + add + lattice->size(d))) % lattice->size(d);

          diff += s-j;
          sum1 += t[X].r[d] - lattice->size(d)/2;
          sum2 += j - lattice->size(d)/2;
        }     
     
        if (diff != 0) {
          hila::output << "Vectorized std gather test error! Node " << mynode() 
                       << " Parity " << parity_name(p) << " direction " << (unsigned)d2 << '\n';
          exit(-1);
        }

        if (sum1 + s_result != 0.0) {
          output0 << "Error in vector sum reduction!  answer " << sum1 + s_result << " should be 0\n";
          exit(-1);
        }

        if (sum2 + s_result != 0.0) {
          output0 << "Error in vector neighbour sum reduction!  answer " << sum2 + s_result << " should be 0\n";
          exit(-1);
        }
        
        t.mark_changed(ALL);
#endif
      }
    }
  }
  #ifdef SPECIAL_BOUNDARY_CONDITIONS
  }
  #endif
}



void test_std_gathers()
{
  gather_test<double>();
  gather_test<int>();
}

