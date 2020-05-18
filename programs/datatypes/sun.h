#ifndef SUN_M
#define SUN_M

#include "../plumbing/defs.h"
#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/vector.h"
#include "../plumbing/mersenne.h" //has to be included
#include <cmath>

// Macros for sped-up operations 

#define CMULJJ(a,b,c) do { (c).re =  (a).re*(b).re - (a).im*(b).im; \
   		        (c).im = -(a).re*(b).im - (a).im*(b).re; } while(0)
#define CMULJJ_ADD(a,b,c) do { (c).re +=  (a).re*(b).re - (a).im*(b).im; \
   		            (c).im += -(a).re*(b).im - (a).im*(b).re; } while(0)
#define CMULJJ_SUB(a,b,c) do { (c).re -=  (a).re*(b).re - (a).im*(b).im; \
   	                (c).im += (a).re*(b).im + (a).im*(b).re; } while(0)

#define su3_vector_crossprod_conj(av, bv, res)  \
do {                                            \
                                                \
    CMULJJ((av).c[1], (bv).c[2], (res).c[0]);         \
    CMULJJ_SUB((av).c[2], (bv).c[1], (res).c[0]);     \
                                                \
    CMULJJ((av).c[2], (bv).c[0], (res).c[1]);         \
    CMULJJ_SUB((av).c[0], (bv).c[2], (res).c[1]);     \
                                                \
    CMULJJ((av).c[0], (bv).c[1], (res).c[2]);         \
    CMULJJ_SUB((av).c[1], (bv).c[0], (res).c[2]);     \
} while (0)

// SU2 matrix multiplication routines ------------------------

#define nn_a(x,y) (x.d*y.a + x.a*y.d - x.b*y.c + x.c*y.b)
#define nn_b(x,y) (x.d*y.b + x.b*y.d - x.c*y.a + x.a*y.c)
#define nn_c(x,y) (x.d*y.c + x.c*y.d - x.a*y.b + x.b*y.a)
#define nn_d(x,y) (x.d*y.d - x.a*y.a - x.b*y.b - x.c*y.c)

#define na_a(x,y) (-x.d*y.a + x.a*y.d + x.b*y.c - x.c*y.b)
#define na_b(x,y) (-x.d*y.b + x.b*y.d + x.c*y.a - x.a*y.c)
#define na_c(x,y) (-x.d*y.c + x.c*y.d + x.a*y.b - x.b*y.a)
#define na_d(x,y) ( x.d*y.d + x.a*y.a + x.b*y.b + x.c*y.c)

#define an_a(x,y) ( x.d*y.a - x.a*y.d + x.b*y.c - x.c*y.b)
#define an_b(x,y) ( x.d*y.b - x.b*y.d + x.c*y.a - x.a*y.c)
#define an_c(x,y) ( x.d*y.c - x.c*y.d + x.a*y.b - x.b*y.a)
#define an_d(x,y) ( x.d*y.d + x.a*y.a + x.b*y.b + x.c*y.c)

#define aa_a(x,y) (-x.d*y.a - x.a*y.d - x.b*y.c + x.c*y.b)
#define aa_b(x,y) (-x.d*y.b - x.b*y.d - x.c*y.a + x.a*y.c)
#define aa_c(x,y) (-x.d*y.c - x.c*y.d - x.a*y.b + x.b*y.a)
#define aa_d(x,y) ( x.d*y.d - x.a*y.a - x.b*y.b - x.c*y.c)


//--------------------------------------------------------

//////////////////////////////////////////////////////////
///
/// SU(N) class
/// Implementations are found in this header file, since 
/// non-specialized templates have to be visible to the 
/// tranlation units that use them.
///
//////////////////////////////////////////////////////////

template<int n, typename radix=double>
class SU : public matrix<n,n,cmplx<radix> >{
    public:
    using base_type = typename base_type_struct<radix>::type;

    SU() = default;

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    #pragma transformer loop_function
    SU(const scalart rhs) {
      for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
        if (i == j) this->c[i][j] = (rhs);
        else this->c[i][j] = (0);
      }
    }

    SU(matrix<n,n,cmplx<radix>> m) {
        for (int j=0; j<n; j++) for (int i=0; i<n; i++){
            this->c[i][j] = m.c[i][j];
        }
    }


    void reunitarize(){ //implement later
        make_unitary();
        fix_det();
    };

    void make_unitary(){};

    void fix_det()
    {
        cmplx<radix> d,factor;
        radix t;
        int i,j;

        d = det(*(this));
        t = d.arg() / static_cast<radix>(n);
        factor = cmplx<radix>( cos(-t), sin(-t) );
        for (j=0; j<n; j++) for (i=0; i<n; i++){
            this->c[j][i] = this->c[j][i]*factor;
        }
    }

    void exp(const int depth = 12){
        matrix<n,n,cmplx<radix>> A, An;
        radix factor = 1;
        A = *this;
        An = A;
        *this = 1;
        for (int k = 1; k<=depth; k++){
            factor = factor/static_cast<radix>(k);
            (*this) += An*factor;
            An *= A;
        }
    }

    //generate random SU(N) element by expanding exp(A), where A is a traceless hermitian matrix. 
    //more iterations are needed to generate larger elements: 12 works well for n < 10. 

    void random(const int depth = 12) {
        matrix<n,n,cmplx<radix>> A, An, res;
        An = 1; 
        res = 1;
        cmplx<radix> tr(1,0), factor(1, 0);
        for (int i = 0; i < n; i++) {
            A.c[i][i] = cmplx<radix>(hila_random(), 0.0);
            for (int j = 0; j < i; j++){
                cmplx<radix> a(static_cast<radix>(hila_random()/n), static_cast<radix>(hila_random()/n));
                A.c[i][j] = a;
                A.c[j][i] = a.conj();
            }
        }
        tr = A.trace()*(static_cast<radix>(1)/static_cast<radix>(n));
        for (int i = 0; i < n; i++){
            A.c[i][i] -= tr; 
        }
        An = A;
        for (int k = 1; k<=depth; k++){
            factor = factor*cmplx<radix>(0, 1)*(static_cast<radix>(1)/static_cast<radix>(k));
            res += An*factor;
            An *= A;
        }
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++){
            (*this).c[i][j] = res.c[i][j];
        }
    }


    // Generate a random algebra (antihermitean) matrix
    void gaussian_algebra(){
      for(int i=0; i<n; i++) {
        for(int j=0; j<i; j++) {
          double a = gaussian_ran();
          double b = gaussian_ran();
          (*this).c[i][j].re = a;
          (*this).c[j][i].re =-a;
          (*this).c[i][j].im = b;
          (*this).c[j][i].im = b;
        }
      }

      for(int i=0; i<n; i++) {
        (*this).c[i][i].re = 0;
        (*this).c[i][i].im = 0;
      }
      for(int i=1; i<n; i++) {
        double a = gaussian_ran()*sqrt(2.0/(i*(i+1)));
        for(int j=0; j<i; j++)
          (*this).c[j][j].im += a;
        (*this).c[i][i].im -= i*a;
      }
    }

    // The norm of an algebra matrix
    radix algebra_norm(){
      double thissum = 0;
      for(int i=0; i<n; i++) {
        for(int j=0; j<i; j++) {
          thissum += (*this).c[i][j].squarenorm();
        }
        double diag = (*this).c[i][i].im;
        thissum += 0.5*diag*diag;
      }
      return thissum;
    }


    // find determinant using LU decomposition. Algorithm: numerical Recipes, 2nd ed. p. 47 ff 
    cmplx<radix> det_lu(){

        int i, imax, j, k;
        radix big, d, temp, dum;
        cmplx<radix> cdum, csum, ctmp1;
        radix vv[n];
        cmplx<radix> a[n][n];
        cmplx<radix> one;

        one=cmplx<radix>(1,0);

        d=1;

        imax = -1;

        for (i=0; i<n; i++) for(j=0; j<n; j++) a[i][j] = this->c[i][j];

        for (i=0; i<n; i++) {
            big = 0;
            for(j=0; j<n; j++) {
            if ((temp = a[i][j].abs()) > big) big = temp;
            }
            if (big == 0.0) exit(1);
            vv[i] = 1.0/big;
        }

        for (j=0; j<n; j++) {
            for (i=0; i<j; i++) {
            csum = a[i][j];
            for (k=0; k<i; k++) {
                csum -= a[i][k]*a[k][j];
            }
            a[i][j] = csum;
            }

            big = 0;
            for (i=j; i<n; i++) {
            csum = a[i][j];
            for (k=0; k<j; k++) {
                csum -= a[i][k]*a[k][j];
            }
            a[i][j] = csum;
            if ((dum = vv[i]*csum.abs()) >= big) {
                big = dum;
                imax = i;
            }
            }

            if (j != imax) {
            for (k=0; k<n; k++) {
            cdum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = cdum;
            }
            d = -d;
            vv[imax] = vv[j];
            }

            if (a[j][j].abs() == static_cast<radix>(0.0)) a[j][j] = cmplx<radix>(1e-20,0);

            if (j != n-1) {
                cdum = one/a[j][j]; //check cmplx division
                for (i=j+1; i<n; i++) {
                    a[i][j] = a[i][j]*cdum;
                }
            }
        }

        csum = cmplx<radix>(d,0.0);
        for (j=0; j<n; j++) {
            csum = csum*a[j][j];
        }

        return (csum);
    }




    static constexpr SU generator(int ng){
      // SUN generators normalized as tr(T^2) = 2
      SU generator = 0;
      if(ng<n-1){
        // Diagonal generators
        double w = sqrt(2.0/((ng+1)*(ng+2)));
        for(int i=0; i<ng+1; i++ ){
          generator.c[i][i].im = w;
        }
        generator.c[ng+1][ng+1].im = -(ng+1)*w;
      } else {
        // Nondiagonal ones. Just run through the indexes and
        // count until they match...
        int k=n-1;
        for( int m1=0; m1<n; m1++) for( int m2=m1+1; m2<n; m2++){
          if( ng == k ){
            generator.c[m1][m2].im = 1;
            generator.c[m2][m1].im = 1;
          } else if( n == k+1 ){
            generator.c[m1][m2].re = 1;
            generator.c[m2][m1].re =-1;
          }
          k+=2;
        }
      }
      return generator;
    }

    static constexpr int generator_count(){
      return n*n-1;
    }

};




template<typename radix> 
class SU2; 

template<typename radix> 
class SU3; 

template<typename radix> 
class SU2adjoint; 

template<typename radix> 
class SU3 : public SU<3,radix> {
};

template<typename radix>
class SU2 { 
    public: 
        using base_type = typename base_type_struct<radix>::type;

        SU2() : a(0), b(0), c(0), d(1) {}

        SU2(radix * vals) : a(vals[0]), b(vals[1]), c(vals[2]), d(vals[3]) { normalize(); }

        SU2(const SU2<radix> & rhs){
            b = rhs.b;				
            a = rhs.a;				
            d = rhs.d;				
            c = rhs.c;           
        };

        friend SU2<radix> operator * (const SU2<radix> & x, const SU2<radix> & y){
            SU2<radix> r;
            r.a = nn_a(x,y); r.b = nn_b(x,y); 
            r.c = nn_c(x,y); r.d = nn_d(x,y);
            return r;
        }
 
        friend SU2<radix> operator * (const SU2<radix> & x, const SU2adjoint<radix> & y){
            SU2<radix> r;
            r.a = na_a(x,y.ref); r.b = na_b(x,y.ref); \
            r.c = na_c(x,y.ref); r.d = na_d(x,y.ref);
            return r;
        }
 
        friend SU2<radix> operator * (const SU2adjoint<radix> & x, const SU2<radix> & y){
            SU2<radix> r;
            r.a = an_a(x.ref,y); r.b = an_b(x.ref,y); \
            r.c = an_c(x.ref,y); r.d = an_d(x.ref,y);
            return r;
        }

        friend SU2<radix> operator * (const SU2adjoint<radix> & x, const SU2adjoint<radix> & y){
            SU2<radix> r;
            r.a = aa_a(x.ref, y.ref); r.b = aa_b(x.ref, y.ref); \
            r.c = aa_c(x.ref, y.ref); r.d = aa_d(x.ref, y.ref);
            return r;
        }

        SU2<radix> & operator = (const SU2<radix> &);
        SU2<radix> & operator = (const SU2adjoint<radix> &);
        SU2<radix> & normalize();
        SU2<radix> & reunitarize();  
        SU2<radix> & random(); 
        SU2<radix> & inv(); 
        SU2adjoint<radix> & adj();         
        radix sqr() const;
        radix tr() const;
        radix det() const; 
        SU2<radix> operator + (const SU2<radix> &); //basic operations. These versions return new matrix as result
        SU2<radix> operator - (const SU2<radix> &);
        SU2<radix> & operator += (const SU2<radix> &); //same ops as above, except store result in lhs matrix
        SU2<radix> & operator -= (const SU2<radix> &);
        SU2<radix> & operator *= (const SU2<radix> &);
        SU2<radix> & operator += (const SU2adjoint<radix> &); //same ops as above, except store result in lhs matrix
        SU2<radix> & operator -= (const SU2adjoint<radix> &);
        SU2<radix> & operator *= (const SU2adjoint<radix> &);
        SU2<radix> & operator *= (const radix &);
        SU2<radix> operator * (const radix &);
    private:
        radix a, b, c, d;  
};


template<typename radix>
class SU2adjoint {
    public:
        SU2adjoint (const SU2<radix> & rhs) : ref(rhs) {} ;
        const SU2<radix> & ref;
    private:
        SU2adjoint(){}
        SU2adjoint<radix> & operator = (const SU2adjoint & rhs){}
};

template<typename radix> 
inline SU2adjoint<radix> conj(SU2<radix> & ref){
  SU2adjoint<radix> result(ref);
  return result;
}

template<typename radix>
SU2adjoint<radix> & SU2<radix>::adj(){
    return SU2adjoint<radix>(*this);
}; 

template<typename radix>
radix SU2<radix>::sqr() const {
    return a*a + b*b + c*c + d*d;
}

template<typename radix>
radix SU2<radix>::det() const {
    return a*a + b*b + c*c + d*d;
}

template<typename radix>
radix SU2<radix>::tr() const {
    return 2*d;
}

template<typename radix>
SU2<radix> & SU2<radix>::normalize(){
    radix sq = sqrt(this->sqr());
    a /= sq;
    b /= sq;
    c /= sq;
    d /= sq;   
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::reunitarize(){
    return this->normalize();
}

template<typename radix>
SU2<radix> & SU2<radix>::random(){
    radix one, two;
    one = gaussian_ran2(&two);
    a = one;
    b = two;
    one = gaussian_ran2(&two);
    c = one;
    d = two;
    return this->normalize();
}

template<typename radix>
SU2<radix> & SU2<radix>::inv(){
    a *= static_cast<radix>(-1);
    b *= static_cast<radix>(-1);
    c *= static_cast<radix>(-1);
    d *= static_cast<radix>(-1);
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::operator = (const SU2<radix> & rhs){
    a = rhs.a;
    b = rhs.b;
    c = rhs.c;
    d = rhs.d;
    return *this;
};

template<typename radix>
SU2<radix> & SU2<radix>::operator = (const SU2adjoint<radix> & rhs){
    a = -rhs.a;
    b = -rhs.b;
    c = -rhs.c;
    d = rhs.d;
    return *this;
};

template<typename radix>
SU2<radix> & SU2<radix>::operator *= (const SU2<radix> & y){
    a = nn_a((*this),y); 
    b = nn_b((*this),y); 
    c = nn_c((*this),y); 
    d = nn_d((*this),y); 
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::operator *= (const SU2adjoint<radix> & y){
    a = na_a((*this),y); 
    b = na_b((*this),y); 
    c = na_c((*this),y); 
    d = na_d((*this),y); 
    return *this;
}

template<typename radix>
SU2<radix> & SU2<radix>::operator *= (const radix & rhs){ 
    a *= rhs;
    b *= rhs;
    c *= rhs;
    d *= rhs;
    return *this;
};

template<typename radix>
SU2<radix> SU2<radix>::operator * (const radix & rhs){ 
    SU2<radix> r;
    r.a = a*rhs;
    r.b = b*rhs;
    r.c = c*rhs;
    r.d = d*rhs;
    return r;
};

template<typename radix>
SU2<radix> SU2<radix>::operator + (const SU2<radix> & y){
    SU2<radix> r;
    r.a = a + y.a; r.b = b + y.b; \
    r.c = c + y.c; r.d = d + y.d;
    return r;
};

template<typename radix>
SU2<radix> SU2<radix>::operator - (const SU2<radix> & y){
    SU2<radix> r;
    r.a = a - y.a; r.b = b - y.b; \
    r.c = c - y.c; r.d = d - y.d;
    return r;
};






/// Project to the antihermitean part of a matrix
template<int N, typename radix>
void project_antihermitean(SU<N,radix> &matrix){
  double tr = 0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<i; j++) {
      double a = 0.5*(matrix.c[i][j].re - matrix.c[j][i].re);
      double b = 0.5*(matrix.c[i][j].im + matrix.c[j][i].im);
      matrix.c[i][j].re = a;
      matrix.c[j][i].re =-a;
      matrix.c[i][j].im = b;
      matrix.c[j][i].im = b;
    }
    tr += matrix.c[i][i].im;
    matrix.c[i][i].re = 0;
  }
  for(int i=0; i<N; i++) {
    matrix.c[i][i].im -= tr/N;
  }
};








template<int n, typename radix=double>
class SU_vector : public vector<n,cmplx<radix>>{
  public:
    using base_type = typename base_type_struct<radix>::type;

    SU_vector() = default;

    template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
    #pragma transformer loop_function
    SU_vector(const scalart rhs) {
      for(int i=0; i<n; i++){
        this->c[i] = (rhs);
      }
    }

    SU_vector(vector<n,cmplx<radix>> m) {
      for(int i=0; i<n; i++){
        this->c[i] = m.c[i];
      }
    }

    inline double rdot(const SU_vector &rhs) const {
      double r = 0;
      for (int i=0; i<n; i++) {
        r += this->c[i].re*rhs.c[i].re;
        r += this->c[i].im*rhs.c[i].im;
      }
      return r;
    }
    
};




#endif 