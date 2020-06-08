#ifndef SU_VEC
#define SU_VEC

#include "../datatypes/cmplx.h"
#include "../datatypes/matrix.h"
#include "plumbing/random.h"

template<int n, typename T>
class vector {
  public:
  T c[n];
  using base_type = typename base_type_struct<T>::type;
  
  vector() = default;

  vector(matrix<1,n,T> m) {
    for (int i=0; i<n; i++){
      c[i] = m.c[i];
    }
  }

  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
  #pragma hila loop_function
  vector & operator= (const scalart rhs) {
    for (int i=0; i<n; i++){
      c[i] = rhs;
    }
    return *this;
  }


  inline void gaussian(){ 
    for (int i = 0; i < n; i++) {
      (*this).c[i].re = gaussian_ran(0.5);
      (*this).c[i].im = gaussian_ran(0.5);
    }
  }

  inline void random(){ 
    (*this).gaussian();
  }

  inline auto norm_sq(){ 
    auto r=norm_squared(c[0]);
    for (int i = 1; i < n; i++) {
      r += norm_squared(c[i]);
    }
    return r;
  }


  inline vector operator-() const {
    vector<n,T> r;
    for (int i = 0; i < n; i++) {
      r.c[i] = -c[i];
    }
    return r;
  }

  #pragma hila loop_function
  vector & operator+=(const vector & rhs){
    for (int i = 0; i < n; i++){
      c[i] += rhs.c[i];
    }
    return *this;
  }
  
  #pragma hila loop_function
  vector & operator-=(const vector & rhs){
    for (int i = 0; i < n; i++){
      c[i] -= rhs.c[i];
    }
    return *this;
  }

  #pragma hila loop_function
  vector & operator-(){
    for (int i = 0; i < n; i++){
      c[i] = -c[i];
    }
    return *this;
  }

  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
  vector & operator*=(const scalart rhs){
    for (int i=0; i<n; i++) {
      c[i] *= rhs;
    }
    return *this;
  }

  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
  vector & operator/=(const scalart rhs){
    for (int i=0; i<n; i++) {
      c[i] /= rhs;
    }
    return *this;
  }


  inline T dot(const vector &rhs) const {
    T r = (0.0);
    for (int i=0; i<n; i++) {
      r += conj(c[i])*rhs.c[i];
    }
    return r;
  }

  inline squarematrix<n,T> outer_product(const vector<n,T> rhs) const {
    squarematrix<n,T> r;
    for(int j=0; j<n; j++) for (int i=0; i<n; i++) {
      r.c[i][j] = c[i] * rhs.c[j].conj();
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




template<int n, typename T>
vector<n,T>  operator*(vector<n,T> lhs, squarematrix<n,T> rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = 0;
    for(int j=0; j<n; j++) {
      r.c[i] += lhs.c[j] * rhs.c[j][i];
    }
  }
  return r;
}

template<int n, typename T>
vector<n,T>  operator*(squarematrix<n,T>  lhs, vector<n,T> rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = 0;
    for(int j=0; j<n; j++) {
      r.c[i] += lhs.c[i][j] * rhs.c[j];
    }
  }
  return r;
}


template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T>  operator*(vector<n,T> lhs, squarematrix<n,scalart> rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = 0;
    for(int j=0; j<n; j++) {
      r.c[i] += lhs.c[j] * rhs.c[j][i];
    }
  }
  return r;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T>  operator*(squarematrix<n,scalart>  lhs, vector<n,T> rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = 0;
    for(int j=0; j<n; j++) {
      r.c[i] += lhs.c[i][j] * rhs.c[j];
    }
  }
  return r;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T> operator*(const scalart &lhs, const vector<n,T> &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs*rhs.c[i];
  }
  return r;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T> operator*(const vector<n,T> &lhs, const scalart &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i]*rhs;
  }
  return r;
}



template<int n, typename T>
vector<n,T> operator*(const T &lhs, const vector<n,T> &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs*rhs.c[i];
  }
  return r;
}

template<int n, typename T>
vector<n,T> operator*(const vector<n,T> &lhs, const T &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i]*rhs;
  }
  return r;
}


template<int n, typename T>
vector<n,T>  operator+(vector<n,T> lhs, vector<n,T> rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i] + rhs.c[i];
  }
  return r;
}

template<int n, typename T>
vector<n,T>  operator-(vector<n,T> lhs, vector<n,T> rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i] - rhs.c[i];
  }
  return r;
}



template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T> operator+(const scalart &lhs, const vector<n,T> &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = rhs+rhs.c[i];
  }
  return r;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T> operator+(const vector<n,T> &lhs, const scalart &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i]+rhs;
  }
  return r;
}


template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T> operator-(const scalart &lhs, const vector<n,T> &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = rhs-rhs.c[i];
  }
  return r;
}

template <int n, typename T, typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 > 
vector<n,T> operator-(const vector<n,T> &lhs, const scalart &rhs){
  vector<n,T>  r;
  for (int i=0; i<n; i++) {
    r.c[i] = lhs.c[i]-rhs;
  }
  return r;
}



template<int n, typename T>
inline auto norm_squared(vector<n,T> & v){
  auto result = norm_squared(v.c[0]);
  for (int i=1; i<n; i++) {
    result += norm_squared(v.c[i]);
  }
  return result;
}




#endif