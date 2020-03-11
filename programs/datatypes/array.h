

template<int n, typename T>
class array {
    public:
        using base_type = typename base_type_struct<T>::type;
        T c[n];
        array() = default;
        void as_matrix();

        array & operator = (const array<n,T> & rhs){
            for (int i = 0; i < n; i++) c[i] = rhs.c[i];
        }

        friend array<n,T> operator + (const array<n,T> & x, const array<n,T> & y){
            array<n,T> r;
            for (int i = 0; i < n; i++) r.c[i] =  x.c[i] + y.c[i];
            return r;
        }

        friend array<n,T> operator - (const array<n,T> & x, const array<n,T> & y){
            array<n,T> r;
            for (int i = 0; i < n; i++) r.c[i] =  x.c[i] + y.c[i];
            return r;
        }

        friend array<n,T> operator * (const array<n,T> & x, const array<n,T> & y){
            array<n,T> r;
            for (int i = 0; i < n; i++) r.c[i] = conj(x.c[i])*y.c[i];
            return r;
        }

        inline T & operator[](int index){
            static_assert(index < n, "array subscript out of bounds!");
            return c[index]; 
        }
};