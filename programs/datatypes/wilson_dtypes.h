
#ifndef WVEC
#define WVEC 

///////////////////////////////////////////////////////
//                 Wilson vectors 
// e.g.                                              */
// wilson_propagator prop;                           */
// prop.c[ci].d[si].d[sf].c[cf]                      */
// ----------------------->    complex               */
// ----------------->          suN_vector            */
// ----------->                wilson_vector         */
// ----->                      spin_wilson_vector    */
// e.g.                                              */
// wilson_matrix matr;                               */
// matr.d[si].c[ci].d[sf].c[cf]                      */
// ----------------------->    complex               */
// ----------------->          suN_vector            */
// ----------->                wilson_vector         */
// ----->                      color_wilson_vector   */

#ifndef GAMMA_DIM
#define GAMMA_DIM 2
#endif

#ifndef COLORVEC_DIM
#define COLORVEC_DIM 2
#endif

template<typename vectortype> 
class wvector {
    using base_type = typename basetypeclass<vectortype>::type;
    vectortype [GAMMA_DIM]d;
};

template<typename vectortype> 
class halfwvector {
    using base_type = typename basetypeclass<vectortype>::type;
    vectortype [GAMMA_DIM/2]d;
};

template<typename T> 
class colorwvector {
    using base_type = typename basetypeclass<T>::type;
    wvector<T> c[COLORVEC_DIM]; 
};

template<typename T> 
class spinwvector {
    using base_type = typename basetypeclass<T>::type;
    wvector<T> d[GAMMA_DIM];
};

template<typename T>
class wmatrix {
    using base_type = typename basetypeclass<T>::type;
    colorwvector<T> d[GAMMA_DIM];
};

template<typename T>
class wpropagator{
    using base_type = typename basetypeclass<T>::type;
    spinwvector<T> c[COLORVEC_DIM];    
};

#endif 
