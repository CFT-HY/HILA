#include "cmplx.h"
#include "general_matrix.h"

//////////////////
///
/// Template class for SU matrices
/// 
/// Uses specialized format for cases n = 2 and n = 3 
/// (parametrized by pauli matrices) for efficiency
/// otherwise behaves like a regular complex matrix
///
//////////////////

//TODO: SU(2) vector class + conversion to and from matrix rep.

template<int n, typename radix>
class SU : public matrix<n,n,cmplx<radix>> {
    public: 
        void reunitarize();
};

template<typename radix>
class SU<2, radix> { 
    public: 

        SU();
        ~SU();
        void normalize(); //normalize elements
        void reunitarize();
        void random(); //generate random SU2 element 
        radix sqr(); //calculate square of all elements
        radix det(); //determinant
        radix tr(); //trace 
        SU<2,radix> & inv(); //invert matrix
        SU<2,radix> operator + (const SU<2,radix> &); //return copies of result 
        SU<2,radix> operator - (const SU<2,radix> &);
        SU<2,radix> operator * (const SU<2,radix> &);
        SU<2,radix> operator + (const radix &);
        SU<2,radix> operator - (const radix &);
        SU<2,radix> operator * (const radix &);
        SU<2,radix> & operator += (const SU<2,radix> &); //same ops as above, except store result in this matrix
        SU<2,radix> & operator -= (const SU<2,radix> &);
        SU<2,radix> & operator *= (const SU<2,radix> &);
        SU<2,radix> & operator += (const radix &);
        SU<2,radix> & operator -= (const radix &);
        SU<2,radix> & operator *= (const radix &);
        SU<2,radix> & operator = (const SU<2,radix> &);

    private:
    
        radix a, b, c, d;      
};
