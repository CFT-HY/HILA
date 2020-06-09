#ifndef REPRESENTATIONS_M
#define REPRESENTATIONS_M

#include "sun.h"


// The adjoint representation
template<int N, typename radix>
class adjoint : public squarematrix<N*N-1, radix> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    using squarematrix<N*N-1, radix>::squarematrix;


    // Info on group generators
    constexpr static int size = N*N-1;
    static constexpr sun generator(int i){
      return sun::generator(i);
    }


    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = -0.5*(generator(i)*m.conjugate()*generator(j)*m).trace().re;
      }
    }
};






template<int N, typename radix>
class antisymmetric : public squarematrix<N*(N-1)/2, cmplx<radix>> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    using squarematrix<N*(N-1)/2, cmplx<radix>>::squarematrix;
    using squarematrix<N*(N-1)/2, cmplx<radix>>::operator =;

    // Info on group generators
    constexpr static int size = N*(N-1)/2;
    static constexpr sun generator(int ng){
      sun generator = 0;
      int k=0;
      for( int m1=0; m1<N; m1++) for( int m2=m1+1; m2<N; m2++){
        if( ng == k ){
          generator.c[m1][m2].re = 1.0/sqrt(2.0);
          generator.c[m2][m1].re =-1.0/sqrt(2.0);
        }
        k++;
      }
      return generator;
    }


    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = -(generator(i)*m*generator(j)*m.transpose()).trace();
      }
    }
};



template<int N, typename radix>
class symmetric : public squarematrix<N*(N+1)/2, cmplx<radix>> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    using squarematrix<N*(N+1)/2, cmplx<radix>>::squarematrix;
    using squarematrix<N*(N+1)/2, cmplx<radix>>::operator =;
    // Info on group generators
    constexpr static int size = N*(N+1)/2;
    static constexpr sun generator(int ng){
      sun generator = 0;
      if(ng < N){
        generator.c[ng][ng].re = 1;
      }
      int k=N;
      for( int m1=0; m1<N; m1++) for( int m2=m1+1; m2<N; m2++){
        if( ng == k ){
          generator.c[m1][m2].re = 1.0/sqrt(2.0);
          generator.c[m2][m1].re = 1.0/sqrt(2.0);
        }
        k++;
      }
      return generator;
    }

    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = (generator(i)*m*generator(j)*m.transpose()).trace();
      }
    }
};




/*
template<typename sun, typename representation>
void represent_gauge_field(
  field<sun> (&u)[NDIM],
  field<representation> (&v)[NDIM]
){
  foralldir(dir){
    onsites(ALL){
      v[dir][X].represent(u[dir][X]);
    }
  }
}


template<typename sun, typename representation>
void project_representation_force(
  field<representation> (&rmom)[NDIM],
  field<sun> (&fmom)[NDIM]
){
  foralldir(dir){
    onsites(ALL){
      fmom[dir][X] = 0;
    }
  }
}
*/

#endif