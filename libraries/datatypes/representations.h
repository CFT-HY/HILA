#ifndef REPRESENTATIONS_H_
#define REPRESENTATIONS_H_

#include "sun.h"

// The adjoint representation
template<int N, typename radix>
class adjointRep : public SquareMatrix<N*N-1, radix> {
  public:
    static_assert(is_arithmetic<radix>::value, "adjointRep<type>: type has to be real");

    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    // Info on group generators
    constexpr static int size = N*N-1;
    static sun generator(int i){
      return sun::generator(i);
    }

    // std ops required for triviality
    adjointRep() = default;
    ~adjointRep() = default;
    adjointRep(const adjointRep &a) = default;

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    #pragma hila loop_function
    adjointRep(const scalart m) : SquareMatrix<size,radix>(m) {}

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    adjointRep(const adjointRep<N,scalart> m) : SquareMatrix<size,scalart>(m) {}

    /* Return a SU(N) generator in the adjoint representation,
     * multiplied by I */
    static adjointRep represented_generator_I(int i){
      static bool initialize = true;
      static adjointRep r_generators[size];
      if(initialize) for(int g=0;g<size;g++){
        r_generators[g] = 0;
        sun tg = sun::generator(g);
        for(int j=0; j<size; j++){
          sun tj = generator(j);
          for(int k=0; k<size; k++){
            sun tk = generator(k);
            cmplx<radix> tr1 = (tg*tj*tk).trace();
            cmplx<radix> tr2 = (tj*tg*tk).trace();
            r_generators[g].e(j,k) = 2*(tr1.im - tr2.im);
          }
        }
        initialize = false;
      }
      return r_generators[i];
    }

    /* Project a matrix into the adjoint representation */
    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).e(i,j) = 2*(m.adjoint()*generator(i)*m*generator(j)).trace().re;
      }
    }

    /* Project a complex adjoint matrix into the algebra and
     * represent as a complex NxN (momentum) matrix */
    static SquareMatrix<N, cmplx<radix>> project_force(
      SquareMatrix<size, cmplx<radix>> rforce
    ){
      SquareMatrix<N, cmplx<radix>> fforce=0;
      for(int g=0; g<size; g++){
        adjointRep rg = represented_generator_I(g);
        radix C = (rg.transpose()*rforce).trace().re;
        fforce += C*sun::generator(g);
      }
      cmplx<radix> ct(0,2.0);
      fforce = fforce*ct;
      project_antihermitean(fforce);
      return fforce;
    }
};






template<int N, typename radix>
class antisymmetric : public SquareMatrix<N*(N-1)/2, cmplx<radix>> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    constexpr static int size = N*(N-1)/2;

    using SquareMatrix<size, cmplx<radix>>::SquareMatrix;
    using SquareMatrix<size, cmplx<radix>>::operator =;

    // Info on group generators

    antisymmetric() = default;
    

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    antisymmetric(const scalart m) : SquareMatrix<size,cmplx<radix>>(m) {}

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    antisymmetric(const antisymmetric<N,scalart> m) {
      for (int j=0; j<size; j++) for (int i=0; i<size; i++){
        this->e(i,j) = m.e(i,j);
      }
    }

    static sun generator(int ng){
      static bool initialize = true;
      static sun generators[size];
      if(initialize) for(int g=0;g<size;g++){
        generators[g] = 0;
        int k=0;
        for( int m1=0; m1<N; m1++) for( int m2=m1+1; m2<N; m2++){
          if( g == k ){
            generators[g].c[m1][m2].re = 0.5;
            generators[g].c[m2][m1].re =-0.5;
          }
          k++;
        }
        initialize = false;
      }
      return generators[ng];
    }

    /* Return a SU(N) generator in the antisymmetric representation */
    static antisymmetric represented_generator_I(int i){
      static bool initialize = true;
      static antisymmetric r_generators[size];
      if(initialize) for(int g=0;g<size;g++){
        r_generators[g] = 0;
        sun tg = sun::generator(g);
        for(int j=0; j<N*(N-1)/2; j++){
          sun tj = generator(j);
          for(int k=0; k<N*(N-1)/2; k++){
            sun tk = generator(k);

            cmplx<radix> tr = (tj*tg*tk).trace();
            r_generators[g].c[j][k] = cmplx<radix>(0,4)*tr;
          }
        }
        initialize = false;
      }
      return r_generators[i];
    }


    /* Project a matrix into the antisymmetric representation */
    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).e(i,j) = 2*(generator(i)*m*generator(j).adjoint()*m.transpose()).trace();
      }
    }

    /* Project a complex antisymmetric matrix into the algebra and
     * represent as a complex NxN (momentum) matrix */
    static SquareMatrix<N, cmplx<radix>> project_force(
      SquareMatrix<size, cmplx<radix>> rforce
    ){
      SquareMatrix<N, cmplx<radix>> fforce=0;
      for(int g=0; g<size; g++){
        antisymmetric rg = represented_generator_I(g);
        radix C = (rg.adjoint()*rforce).trace().re;
        fforce += C*sun::generator(g);
      }
      cmplx<radix> ct(0,-2.0);
      fforce = fforce*ct;
      project_antihermitean(fforce);
      return fforce;
    }
};



template<int N, typename radix>
class symmetric : public SquareMatrix<N*(N+1)/2, cmplx<radix>> {
  public:
    using base_type = typename base_type_struct<radix>::type;
    using sun = SU<N,radix>;

    using SquareMatrix<N*(N+1)/2, cmplx<radix>>::SquareMatrix;
    using SquareMatrix<N*(N+1)/2, cmplx<radix>>::operator =;

    // Info on group generators
    constexpr static int size = N*(N+1)/2;

    symmetric() = default;

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    symmetric(const scalart m) {
      for (int j=0; j<size; j++) { 
        for (int i=0; i<size; i++){
          this->e(i,j) = 0;
        }
        this->c[j][j] = m;
      }
    }
    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    symmetric(const symmetric<N,scalart> m) {
      for (int j=0; j<size; j++) for (int i=0; i<size; i++){
        this->e(i,j) = m.e(i,j);
      }
    }


    static sun generator(int ng){
      static bool initialize = true;
      static sun generators[size];
      if(initialize) for(int g=0;g<size;g++){
        generators[g] = 0;
        if(g < N){
          generators[g].c[g][g].re = sqrt(0.5);
        }
        int k=N;
        for( int m1=0; m1<N; m1++) for( int m2=m1+1; m2<N; m2++){
          if( g == k ){
            generators[g].c[m1][m2].re = 0.5;
            generators[g].c[m2][m1].re = 0.5;
          }
          k++;
        }
        initialize = false;
      }
      return generators[ng];
    }

    /* Return a SU(N) generator in the symmetric representation */
    static symmetric represented_generator_I(int i){
      static bool initialize = true;
      static symmetric r_generators[size];
      if(initialize) for(int g=0;g<size;g++){
        r_generators[g] = 0;
        sun tg = sun::generator(g);
        for(int j=0; j<N*(N+1)/2; j++){
          sun tj = generator(j);
          for(int k=0; k<N*(N+1)/2; k++){
            sun tk = generator(k);

            cmplx<radix> tr = (tj*tg*tk).trace();
            r_generators[g].c[j][k] = cmplx<radix>(0,4)*tr;
          }
        }
        initialize = false;
      }
      return r_generators[i];
    }

    /* Project a matrix into the symmetric representation */
    void represent(sun &m){
      for(int i=0; i<size; i++) for(int j=0; j<size; j++){
        (*this).e(i,j) = 2*(generator(i)*m*generator(j)*m.transpose()).trace();
      }
    }

    /* Project a complex symmetric matrix into the algebra and
     * represent as a complex NxN (momentum) matrix */
    static SquareMatrix<N, cmplx<radix>> project_force(
      SquareMatrix<size, cmplx<radix>> rforce
    ){
      SquareMatrix<N, cmplx<radix>> fforce=0;
      for(int g=0; g<size; g++){
        symmetric rg = represented_generator_I(g);
        radix C = (rg*rforce).trace().re;
        fforce += C*sun::generator(g);
      }
      cmplx<radix> ct(0,-2.0);
      fforce = fforce*ct;
      project_antihermitean(fforce);
      return fforce;
    }
};




#endif