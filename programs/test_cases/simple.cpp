#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <immintrin.h>

#include "../vectorclass/vectorclass.h"

class c {
  private:
  public:
    __m256d _c;
    c(){};
    c(__m256d x){
      _c=x;
    }
    c(double x){
      _c = _mm256_broadcast_sd(&x);
    }
    c & operator=(__m256d const &x){
      this->_c = x;
      return *this;
    }
    c & operator=(double const &x){
      this->_c = _mm256_broadcast_sd(&x);
      return *this;
    }
};


template< class T >
struct is_arithmetic : std::integral_constant<
  bool,
  std::is_arithmetic<T>::value ||
  std::is_same<T,Vec8i>::value ||
  std::is_same<T,c>::value
> {};

#include "../datatypes/general_matrix.h"


class fs{
  public:

    // Array of structures implementation
    matrix<3, 3, cmplx<Vec8i> > * fieldbuf = NULL;

    void allocate_field( const int field_alloc_size ) {
      fieldbuf = (matrix<3, 3, cmplx<Vec8i> >*) aligned_alloc( 32, sizeof(matrix<3, 3, cmplx<Vec8i> >) * field_alloc_size);
      #pragma acc enter data create(fieldbuf)
    }

    void free_field() {
      #pragma acc exit data delete(fieldbuf)
      free((void *)fieldbuf);
      fieldbuf = nullptr;
    }

    //-- #pragma transformer loop_function
    matrix<3, 3, cmplx<Vec8i> > get(const int i, const int field_alloc_size) const
    {
      // There is some problem with directly assigning intrinsic vectors, at least.
      // This is a universal workaround, but could be fixed by assigning element
      // by element
      matrix<3, 3, cmplx<Vec8i> > value = fieldbuf[i];
      //matrix<3, 3, cmplx<Vec8i> > value;
      //std::memcpy( &value, fieldbuf+i, sizeof(matrix<3, 3, cmplx<Vec8i> >) );
      return value;
    }

    //-- #pragma transformer loop_function
    void set(const matrix<3, 3, cmplx<Vec8i> > value, const int i, const int field_alloc_size) 
    {
      //std::memcpy( fieldbuf+i, &value, sizeof(T) );
      fieldbuf[i] = value;
    }
};


int main(int argc, char **argv)
{
  fs a;
  matrix<3, 3, cmplx<Vec8i> > x;

  a.allocate_field(1024);
  
  for(int i=0;i<1024;i++){
    x = a.get(i,1024);
    x = x+x;
    a.set(x,i,1024);
  }
 
  return 0;
}

