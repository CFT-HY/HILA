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
    static sun generator(int i){
      return sun::generator(i);
    }

    /* Return a SU(N) generator in the adjoint representation,
     * multiplied by I */
    static adjoint represented_generator_I(int i){
      static bool initialize = true;
      static adjoint r_generators[N*N-1];
      if(initialize) for(int g=0;g<N*N-1;g++){
        r_generators[g] = 0;
        sun tg = sun::generator(g);
        for(int j=0; j<N*N-1; j++){
          sun tj = generator(j);
          for(int k=0; k<N*N-1; k++){
            sun tk = generator(k);
            cmplx<radix> tr1 = (tg*tj*tk).trace();
            cmplx<radix> tr2 = (tj*tg*tk).trace();
            r_generators[g].c[j][k] = 2*(tr1.im - tr2.im);
          }
        }
        initialize = false;
      }
      return r_generators[i];
    }

    /* Project a matrix into the adjoint representation */
    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = 2*(generator(i)*m.conjugate()*generator(j)*m).trace().re;
      }
    }

    /* Project a complex adjoint matrix into the algebra and
     * represent as a complex NxN (momentum) matrix */
    static squarematrix<N, cmplx<radix>> project_force(
      squarematrix<size, cmplx<radix>> rforce
    ){
      squarematrix<N, cmplx<radix>> fforce=0;
      for(int g=0; g<N*N-1; g++){
        adjoint rg = represented_generator_I(g);
        radix C = (rg.transpose()*rforce).trace().re;
        fforce += C*sun::generator(g);
      }
      cmplx<radix> ct(0,-2);
      fforce = fforce*ct;
      project_antihermitean(fforce);
      return fforce;
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

    /* Return a SU(N) generator in the antisymmetric representation */
    static constexpr antisymmetric represented_generator(int i){
      antisymmetric g;
      sun ti = sun::generator(i);
      for(int j=0; j<N*(N-1)/2; j++){
        sun tj = generator(j);
        for(int k=0; k<N*(N-1)/2; k++){
          sun tk = generator(k);
          
          cmplx<radix> tr1 = (ti*tj*tk).trace();
          cmplx<radix> tr2 = (tj*ti*tk).trace();
          g.c[j][k] = tr1 - tr2;
        }
      }
      return g;
    }


    /* Project a matrix into the antisymmetric representation */
    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = -(generator(i)*m*generator(j)*m.transpose()).trace();
      }
    }

    /* Project a complex antisymmetric matrix into the algebra and
     * represent as a complex NxN (momentum) matrix */
    static squarematrix<N, cmplx<radix>> project_force(
      squarematrix<size, cmplx<radix>> rforce
    ){
      squarematrix<N, cmplx<radix>> fforce=0;
      for(int g=0; g<N*N-1; g++){
        antisymmetric rg = represented_generator(g);
        radix C = -(rg.conjugate()*rforce).trace().im;
        fforce += C*sun::generator(g);
      }
      cmplx<radix> ct(0,2);
      fforce = fforce*ct;
      project_antihermitean(fforce);
      return fforce;
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

    /* Return a SU(N) generator in the symmetric representation */
    static constexpr symmetric represented_generator(int i){
      symmetric g;
      sun ti = sun::generator(i);
      for(int j=0; j<N*(N+1)/2; j++){
        sun tj = generator(j);
        for(int k=0; k<N*(N+1)/2; k++){
          sun tk = generator(k);
          
          cmplx<radix> tr1 = (ti*tj*tk).trace();
          cmplx<radix> tr2 = (tj*ti*tk).trace();
          g.c[j][k] = tr1 - tr2;
        }
      }
      return g;
    }

    /* Project a matrix into the symmetric representation */
    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).c[i][j] = (generator(i)*m*generator(j)*m.transpose()).trace();
      }
    }

    /* Project a complex symmetric matrix into the algebra and
     * represent as a complex NxN (momentum) matrix */
    static squarematrix<N, cmplx<radix>> project_force(
      squarematrix<size, cmplx<radix>> rforce
    ){
      squarematrix<N, cmplx<radix>> fforce=0;
      for(int g=0; g<N*N-1; g++){
        symmetric rg = represented_generator(g);
        radix C = -(rg.conjugate()*rforce).trace().im;
        fforce += C*sun::generator(g);
      }
      cmplx<radix> ct(0,2);
      fforce = fforce*ct;
      project_antihermitean(fforce);
      return fforce;
    }
};




#endif