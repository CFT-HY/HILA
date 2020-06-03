#ifndef REPRESENTATIONS_M
#define REPRESENTATIONS_M

#include "sun.h"


// The adjoint representation
template<int N, typename radix>
class adjoint : public squarematrix<N*N-1, radix> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    // Info on group generators
    constexpr static int size = N*N-1;
    static constexpr sun generator(int i){
      return sun::generator(i);
    }


    adjoint() = default;

    adjoint(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = 2*(generator(i)*m.conjugate()*generator(j)*m).trace().re;
      }
    }

    adjoint represent(sun &m){
      adjoint r(m);
      return r;
    }
};



template<int N, typename radix>
class antisymmetric : public squarematrix<N*(N-1)/2, cmplx<radix>> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    // Info on group generators
    constexpr static int size = N*(N-1)/2;
    static constexpr sun generator(int ng){
      sun generator = 0;
      int k=0;
      for( int m1=0; m1<N; m1++) for( int m2=m1+1; m2<N; m2++){
        if( ng == k ){
          generator.c[m1][m2].re = 1;
          generator.c[m2][m1].re =-1;
        }
        k++;
      }
      return generator;
    }


    antisymmetric() = default;

    antisymmetric(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = 2*(generator(i)*m.conjugate()*generator(j)*m).trace();
      }
    }

    antisymmetric represent(sun &m){
      antisymmetric r(m);
      return r;
    }
};


template<int N, typename radix>
class symmetric : public squarematrix<N*(N+1)/2, cmplx<radix>> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    // Info on group generators
    constexpr static int size = N*(N+1)/2;
    static constexpr sun generator(int ng){
      sun generator = 0;
      if(ng < N){
        generator.c[ng][ng].re = sqrt(2.0);
      }
      int k=N+1;
      for( int m1=0; m1<N; m1++) for( int m2=m1+1; m2<N; m2++){
        if( ng == k ){
          generator.c[m1][m2].re = 1;
          generator.c[m2][m1].re = 1;
        }
        k++;
      }
      return generator;
    }


    symmetric() = default;

    symmetric(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = 2*(generator(i)*m.conjugate()*generator(j)*m).trace();
      }
    }

    symmetric represent(sun &m){
      symmetric r(m);
      return r;
    }
};





#endif