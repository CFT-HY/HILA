//////////////////////////////////////////////////////////////////////////////
/// Test the standard gather here
//////////////////////////////////////////////////////////////////////////////

#include "defs.h"
#include "coordinates.h"
#include "lattice.h"

#if 0
// Declare the test class here, element type definitions must be before field.h
template <typename T>
// #pragma hila dump ast
struct test_struct {
  T r[NDIM];

  using base_type = typename base_type_struct<T>::type;

  // define unary - for antiperiodic b.c.
  #pragma hila loop_function  //TODO
  test_struct<T> operator-() const { test_struct<T> t; foralldir(d) t.r[d] = -r[d]; return t; }

};

#endif

#include "field.h"
#include "datatypes/matrix.h"


template <typename T>
using test_struct = Vector<NDIM,T>;

constexpr direction last_dir = (direction)(NDIM-1);

template <typename T>
void gather_test() {

  extern lattice_struct * lattice;
  Field<test_struct<T>> t;
  
  onsites(ALL) {
    t[X] = X.coordinates();
  }

  coordinate_vector lsize = lattice->size();
#if defined(SPECIAL_BOUNDARY_CONDITIONS)
  for (boundary_condition_t bc : {boundary_condition_t::PERIODIC, boundary_condition_t::ANTIPERIODIC}) {
      t.set_boundary_condition(last_dir,bc);
#endif
  for (parity p : {EVEN,ODD,ALL}) {

    foralldir(d) {
      int size_d = lattice->size(d);
      int size_t = lattice->size(last_dir);
      for (direction d2 : {d,-d}) {
      
        T diff = 0;
        int add;
        double sum1 = 0, sum2 = 0;  // use double to accumulate ints, should be accurate
        if (is_up_dir(d2)) add = 1; else add = -1;
        onsites(p) {
          auto n = t[X+d2];
#if defined(SPECIAL_BOUNDARY_CONDITIONS)
          if (bc == boundary_condition_t::ANTIPERIODIC &&
              (( X.coordinate(last_dir) == 0 && d2 == -last_dir) || 
               ( X.coordinate(last_dir) == lattice->size(last_dir)-1 && d2 == last_dir))) {
            n = -n;
          }
#endif

          T j = n[d];
          T t_r = t[X][d];
          T s = mod(t_r + add, size_d);
          
          sum2 += n[d] - size_d/2.0;
          sum1 += t[X][d] - size_d/2.0;

          T lv = s-j;
          T a = 0;
          foralldir(dir) if (dir != d) a+= n[dir] - t[X][dir];
          
          #ifndef CUDA
          if (lv != 0 || a != 0) {
            hila::output << "Error in gather test at " << X.coordinates() << " direction " << d2 
                         << " parity " << (int)p << '\n';
            hila::output << "Fetched element t[X+d2][d] = " << j << " should be " << s << " perp diff is " << a << '\n';
            hila::output << "This element - neighbour element:  ";
            for (int loop=0; loop<NDIM; loop++) hila::output << t[X][loop] << ' ';
            hila::output << " - ";
            for (int loop=0; loop<NDIM; loop++) hila::output << n[loop] << ' ';
            
            hila::output << '\n';
          }
          #endif
          assert((lv == 0 || a == 0) && "Test gathers");
        }

        double s_result;
        if (p == ALL) 
          s_result = lattice->volume()/2.0;
        else 
          s_result = lattice->volume()/4.0;

        if (sum1 + s_result != 0.0) {
          output0 << "Error in sum reduction!  answer " << sum1 + s_result 
                  << " should be 0, parity " << (int)p << ", direction " << (int)d2;
          #ifdef SPECIAL_BOUNDARY_CONDITIONS
          output0 << ", boundary condition " << (int)bc;
          #endif
          output0 << '\n';
          hila::terminate(-1);
        }

        if (sum2 + s_result != 0.0) {
          output0 << "Error in neighbour sum2 reduction!  answer " << sum2 + s_result 
                  << " should be 0, parity " << (int)p << ", direction " << (int)d2;
          #ifdef SPECIAL_BOUNDARY_CONDITIONS
          output0 << ", boundary condition " << (int)bc;
          #endif
          output0 << '\n';
          hila::terminate(-1);
        }


        t.mark_changed(ALL);  // force fetching, test it too

      }
#ifdef VECTORIZED
      // above is not vectorized, so do it also in vec way
      // make it simple to keep it vectorized
        
      T difT = 0;
      Field<T> f = 1;
      onsites(p) {
        difT += f[X+d] + f[X-d] - 2*f[X];
      }     
     
      if (difT != 0) {
        hila::output << "Vectorized std gather test error! Node " << hila::myrank() 
                     << " Parity " << parity_name(p) << " direction " << (unsigned)d << '\n';
        exit(1);
      }
    
      t.mark_changed(ALL);
#endif
    
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

  timestamp("Communication tests done");
  print_dashed_line();
}

