#include "sun.h"
#include "cmplx.h"

template<typename T, int n>
class array {
    public:
        T c[n];
        array(){};
        void as_matrix();

        friend array<T, n> operator + (const array<T, n> & x, const array<T, n> & y){
            array<T, n> r;
            for (int i = 0; i < n; i++) r[i] =  typeConj(x[i])*y[i];
            return r;
        }

        friend array<T, n> operator - (const array<T, n> & x, const array<T, n> & y){
            array<T, n> r;
            return r;
        }

        friend array<T, n> operator * (const array<T, n> & x, const array<T, n> & y){
            array<T, n> r;
            for (int i = 0; i < n; i++) r[i] =  typeConj(x[i])*y[i];
            return r;
        }

        inline T & operator[](int index){
            static_assert(i < n, "array subscript out of bounds!");
            return c[i]; 
        }
}