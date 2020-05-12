#ifndef WILSON_VECTOR_H
#define WILSON_VECTOR_H

#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/sun.h"
#include "../datatypes/sun_vector.h"
#include "../plumbing/random.h"



#if (NDIM==5 || NDIM==4) 
#define Gammadim 4
#define NGamma 5
enum class gamma_matrix : unsigned { g0 = 0, g1, g2, g3, g5};

constexpr gamma_matrix gamma0 = gamma_matrix::g0;
constexpr gamma_matrix gamma1 = gamma_matrix::g1;
constexpr gamma_matrix gamma2 = gamma_matrix::g2;
constexpr gamma_matrix gamma3 = gamma_matrix::g3;
constexpr gamma_matrix gamma5 = gamma_matrix::g5;

#elif (NDIM==3 || NDIM==2)
#define Gammadim 2
#define NGamma 3
enum class gamma_matrix : unsigned { g0 = 0, g1, g2};

constexpr gamma_matrix gamma0 = gamma_matrix::g0;
constexpr gamma_matrix gamma1 = gamma_matrix::g1;
constexpr gamma_matrix gamma2 = gamma_matrix::g2;

#else
static_assert(false, "Wilson fermions only implemented for 1 < NDIM < 6");
#endif




template<int n, typename radix=double>
class Wilson_vector {
  public:
  SU_vector<n, radix> c[Gammadim];
  using base_type = typename base_type_struct<radix>::type;
  
  Wilson_vector() = default;

  Wilson_vector(SU_vector<n, radix> m) {
    for (int i=0; i<Gammadim; i++){
      c[i] = m;
    }
  }

  Wilson_vector(radix m) {
    for (int i=0; i<Gammadim; i++){
      c[i] = m;
    }
  }

  void gaussian(){ 
    for (int i=0; i<Gammadim; i++){
      c[i].gaussian();
    }
  }


  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
  #pragma transformer loop_function
  Wilson_vector & operator= (const scalart rhs) {
    for (int i=0; i<Gammadim; i++){
      c[i] = rhs;
    }
    return *this;
  }

  radix norm_sq(){ 
    radix r=0;
    for (int i = 0; i < Gammadim; i++) {
      r += c[i].norm_sq();
    }
    return r;
  }

  std::string str() const {
    std::string text = "";
    for (int i=0; i<Gammadim; i++){
      text += c[i].str() + "\n"; 
    }
    text += "\n";
    return text;
  }
};



template<int n, typename radix>
cmplx<radix> operator*(const Wilson_vector<n,radix> lhs, const Wilson_vector<n,radix> rhs){
  cmplx<radix> r = (0.0);
  for (int i=0; i<Gammadim; i++) {
    r += lhs.c[i]*rhs.c[i];
  }
  return r;
}


template<int n, typename radix>
Wilson_vector<n,radix> operator+(const Wilson_vector<n,radix> lhs, const Wilson_vector<n,radix> rhs){
  Wilson_vector<n,radix>  r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs.c[i] + rhs.c[i];
  }
  return r;
}

template<int n, typename radix>
Wilson_vector<n,radix> operator-(const Wilson_vector<n,radix> lhs, const Wilson_vector<n,radix> rhs){
  Wilson_vector<n,radix>  r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs.c[i] - rhs.c[i];
  }
  return r;
}



template<int n, typename radix>
Wilson_vector<n,radix> operator*(const radix lhs, const Wilson_vector<n,radix> &rhs){
  Wilson_vector<n,radix>  r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs * rhs.c[i];
  }
  return r;
}

template<int n, typename radix>
Wilson_vector<n,radix>  operator*(const Wilson_vector<n,radix> lhs, const radix rhs){
  Wilson_vector<n,radix>  r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs.c[i] * rhs;
  }
  return r;
}


/// Returns an SUN matrix, which is the sum of the outer products
// of the SUN vectors c
template<int n, typename radix>
SU<n,radix> outer_product(const Wilson_vector<n,radix> lhs, const Wilson_vector<n,radix> rhs){
  SU<n,radix> r = 0;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] += lhs.c[i] * rhs.c[i];
  }
  return r;
}





// Multiplication with gamma matrices

#if (Gammadim==4) 

template<int n, typename radix>
Wilson_vector<n,radix> operator*(const gamma_matrix gamma, const Wilson_vector<n,radix> rhs){
  Wilson_vector<n,radix>  r;
  switch(gamma) {
    case gamma0:
      r.c[0] = rhs.c[2]; r.c[1] = rhs.c[3];
      r.c[2] = rhs.c[0]; r.c[3] = rhs.c[1];
      break;
    case gamma1:
      r.c[0] = cmplx(0,1)*rhs.c[3]; r.c[1] = cmplx(0,1)*rhs.c[2];
      r.c[2] = cmplx(0,-1)*rhs.c[1]; r.c[3] = cmplx(0,-1)*rhs.c[0];
      break;
    case gamma2:
      r.c[0] = -rhs.c[3]; r.c[1] = rhs.c[2];
      r.c[2] = rhs.c[1]; r.c[3] = -rhs.c[0];
      break;
    case gamma3:
      r.c[0] = cmplx(0,1)*rhs.c[2]; r.c[1] = cmplx(0,-1)*rhs.c[3];
      r.c[2] = cmplx(0,-1)*rhs.c[0]; r.c[3] = cmplx(0,1)*rhs.c[1];
      break;
    case gamma5:
      r.c[0] = rhs.c[0]; r.c[1] = rhs.c[1];
      r.c[2] = -rhs.c[2]; r.c[3] = -rhs.c[3];
      break;
  }
  return r;
}
  

#elif (Gammadim==2)



#endif






#endif