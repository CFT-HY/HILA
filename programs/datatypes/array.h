

template<typename T, int n>
class array {
    public:
        T c[n];
        array() = default;
        void as_matrix();

        array & operator = (const array<T,n> & rhs){
            for (int i = 0; i < n; i++) r.c[i] = rhs.c[i];
        }

        friend array<T, n> operator + (const array<T, n> & x, const array<T, n> & y){
            array<T, n> r;
            for (int i = 0; i < n; i++) r.c[i] =  x.c[i] + y.c[i];
            return r;
        }

        friend array<T, n> operator - (const array<T, n> & x, const array<T, n> & y){
            array<T, n> r;
            for (int i = 0; i < n; i++) r.c[i] =  x.c[i] + y.c[i];
            return r;
        }

        friend array<T, n> operator * (const array<T, n> & x, const array<T, n> & y){
            array<T, n> r;
            for (int i = 0; i < n; i++) r.c[i] = conj(x.c[i])*y.c[i];
            return r;
        }

        inline T & operator[](int index){
            static_assert(i < n, "array subscript out of bounds!");
            return c[i]; 
        }
};