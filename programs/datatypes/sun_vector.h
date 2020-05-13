#ifndef SU_VEC
#define SU_VEC

#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../plumbing/random.h"

template<int n, typename radix=double>
class SU_vector {
  public:
  cmplx<radix> c[n];
  using base_type = typename base_type_struct<radix>::type;
  
  SU_vector() = default;

  SU_vector(matrix<1,n,cmplx<radix>> m) {
    for (int i=0; i<n; i++){
      c[i] = m.c[i];
    }
  }

  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
  #pragma transformer loop_function
  SU_vector & operator= (const scalart rhs) {
    for (int i=0; i<n; i++){
      c[i] = rhs;
    }
    return *this;
  }


  inline void gaussian(){ 
    for (int i = 0; i < n; i++) {
      (*this).c[i].re = gaussian_ran(1.0);
      (*this).c[i].im = gaussian_ran(1.0);
    }
  }

  inline radix norm_sq(){ 
    radix r=0;
    for (int i = 0; i < n; i++) {
      r += c[i].re * c[i].re;
      r += c[i].im * c[i].im;
    }
    return r;
  }


  inline SU_vector operator-() const {
    SU_vector<n,radix> r;
    for (int i = 0; i < n; i++) {
      r.c[i] = -c[i];
    }
    return r;
  }

  inline cmplx<radix> dot(const SU_vector &rhs) const {
    cmplx<radix> r = (0.0);
    for (int i=0; i<n; i++) {
      r += conj(c[i])*rhs.c[i];
    }
    return r;
  }

  inline radix rdot(const SU_vector &rhs) const {
    radix r = (0.0);
    for (int i=0; i<n; i++) {
      r += c[i].re*rhs.c[i].re + c[i].im*rhs.c[i].im;
    }
    return r;
  }

  std::string str() const {
    std::string text = "";
    for (int i=0; i<n; i++){
      text += c[i].str() + " "; 
    }
    return text;
  }

};




template<int n, typename radix>
SU_vector<n,radix>  operator*(SU_vector<n,radix> lhs, matrix<n,n,cmplx<radix>> rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = 0;
    for(int j=0; j<n; j++) {
      r.c[i] += lhs.c[j] * rhs.c[j][i];
    }
  }
  return r;
}

template<int n, typename radix>
SU_vector<n,radix>  operator*(matrix<n,n,cmplx<radix>>  lhs, SU_vector<n,radix> rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = 0;
    for(int j=0; j<n; j++) {
      r.c[i] += lhs.c[i][j] * rhs.c[j];
    }
  }
  return r;
}


template<int n, typename radix>
SU_vector<n,radix> operator*(const cmplx<radix> &lhs, const SU_vector<n,radix> &rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs*rhs.c[i];
  }
  return r;
}

template<int n, typename radix>
SU_vector<n,radix> operator*(const SU_vector<n,radix> &lhs, const cmplx<radix> &rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i]*rhs;
  }
  return r;
}



template<int n, typename radix>
SU_vector<n,radix>  operator+(SU_vector<n,radix> lhs, SU_vector<n,radix> rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i] + rhs.c[i];
  }
  return r;
}

template<int n, typename radix>
SU_vector<n,radix>  operator-(SU_vector<n,radix> lhs, SU_vector<n,radix> rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i] - rhs.c[i];
  }
  return r;
}



template<int n, typename radix>
SU_vector<n,radix>  operator*(radix lhs, SU_vector<n,radix> rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs * rhs.c[i];
  }
  return r;
}

template<int n, typename radix>
SU_vector<n,radix>  operator*(SU_vector<n,radix> lhs, radix rhs){
  SU_vector<n,radix>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i] * rhs;
  }
  return r;
}


template<int n, typename radix>
matrix<n,n,cmplx<radix>> outer_product(SU_vector<n,radix> lhs, SU_vector<n,radix> rhs){
  matrix<n,n,cmplx<radix>> r;
  for(int j=0; j<n; j++) for (int i=0; i<n; i++) {
    r.c[j][i] = lhs.c[i] * rhs.c[j].conj();
  }
  return r;
}



#endif