
#ifndef WVEC
#define WVEC 

#ifndef GAMMA_DIM
#define GAMMA_DIM 2
#endif

#ifndef COLORVEC_DIM
#define COLORVEC_DIM 2
#endif

template<typename T>
class vector {
    vector();
}; 

template<typename T> 
class wvector {
    vector<T> [GAMMA_DIM]d;
};

class halfwvector {
    vector<T> [GAMMA_DIM/2]d;
};

class colorwvector; //wilson_vector c[COLORVEC_DIM]
class spinwvector; //wilson_vector d[GAMMA_DIM]
class wmatrix; //color_wilson_vector d[GAMMA_DIM]
class wpropagator; //spin_wilson_vector c[COLORVEC_DIM]

#endif 