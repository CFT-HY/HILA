#ifndef WILSON_VECTOR_H
#define WILSON_VECTOR_H

#include "../datatypes/cmplx.h"
#include "../datatypes/general_matrix.h"
#include "../datatypes/sun.h"
#include "../datatypes/vector.h"
#include "../plumbing/coordinates.h"
#include "../plumbing/random.h"



#if (NDIM==5 || NDIM==4) 
#define Gammadim 4
#define NGamma 5
enum class gamma_matrix_type : unsigned { g0 = 0, g1, g2, g3, g5};

constexpr gamma_matrix_type gamma0 = gamma_matrix_type::g0;
constexpr gamma_matrix_type gamma1 = gamma_matrix_type::g1;
constexpr gamma_matrix_type gamma2 = gamma_matrix_type::g2;
constexpr gamma_matrix_type gamma3 = gamma_matrix_type::g3;
constexpr gamma_matrix_type gamma5 = gamma_matrix_type::g5;

// List order of directions for conveniance
gamma_matrix_type gamma_matrix[5] = {gamma1, gamma2, gamma3, gamma0, gamma5};

#elif (NDIM==3 || NDIM==2)
#define Gammadim 2
#define NGamma 3
enum class gamma_matrix_type : unsigned { g0 = 0, g1, g2};

constexpr gamma_matrix_type gamma0 = gamma_matrix_type::g0;
constexpr gamma_matrix_type gamma1 = gamma_matrix_type::g1;
constexpr gamma_matrix_type gamma2 = gamma_matrix_type::g2;

// List order of directions for conveniance
gamma_matrix_type gamma_matrix[5] = {gamma0, gamma1, gamma2};

#else
static_assert(false, "Wilson fermions only implemented for 1 < NDIM < 6");
#endif




template<typename vector_type>
class Wilson_vector {
  public:
  vector_type c[Gammadim];
  using base_type = typename base_type_struct<vector_type>::type;
  
  Wilson_vector() = default;

  Wilson_vector(vector_type m) {
    for (int i=0; i<Gammadim; i++){
      c[i] = m;
    }
  }

  template <typename scalart, std::enable_if_t<is_arithmetic<scalart>::value, int> = 0 >  
  #pragma transformer loop_function
  Wilson_vector(const scalart rhs) {
    for(int i=0; i<Gammadim; i++){
      c[i] = rhs;
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

  #pragma transformer loop_function
  Wilson_vector & operator+=(const Wilson_vector & rhs){
    for (int i = 0; i < Gammadim; i++){
      c[i] += rhs.c[i];
    }
    return *this;
  }

  Wilson_vector & operator-=(const Wilson_vector & rhs){
    for (int i = 0; i < Gammadim; i++){
      c[i] -= rhs.c[i];
    }
    return *this;
  }

  Wilson_vector operator-() const {
    Wilson_vector r;
    for (int i = 0; i < Gammadim; i++){
      r.c[i] = -c[i];
    }
    return r;
  }

  inline auto norm_sq(){ 
    auto r=c[0].norm_sq();
    for (int i = 1; i < Gammadim; i++) {
      r += c[i].norm_sq();
    }
    return r;
  }

  inline auto dot(const Wilson_vector &rhs) const {
    auto r = c[0].dot(rhs.c[0]);
    for (int i=1; i<Gammadim; i++) {
      r += c[i].dot(rhs.c[i]);
    }
    return r;
  }

  inline double rdot(const Wilson_vector &rhs) const {
    double r = (0.0);
    for (int i=0; i<Gammadim; i++) {
      r += c[i].rdot(rhs.c[i]);
    }
    return r;
  }

  /// Returns an SUN matrix, which is the sum of the outer products
  // of the SUN vectors c
  inline auto outer_product(const Wilson_vector rhs) const{
    auto r = c[0].outer_product(rhs.c[0]);
    for (int i=1; i<Gammadim; i++) {
      r += c[i].outer_product(rhs.c[i]);
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




template<typename vector, typename T>
Wilson_vector<vector> operator*(const T lhs, const Wilson_vector<vector> rhs){
  Wilson_vector<vector> r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs*rhs.c[i];
  }
  return r;
}

template<typename vector, typename T>
Wilson_vector<vector> operator*(const Wilson_vector<vector> lhs, const T rhs){
  Wilson_vector<vector> r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs.c[i]*rhs;
  }
  return r;
}


template<typename vector>
Wilson_vector<vector> operator+(const Wilson_vector<vector>  lhs, const Wilson_vector<vector>  rhs){
  Wilson_vector<vector> r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs.c[i] + rhs.c[i];
  }
  return r;
}

template<typename vector>
Wilson_vector<vector> operator-(const Wilson_vector<vector> lhs, const Wilson_vector<vector> rhs){
  Wilson_vector<vector> r;
  for (int i=0; i<Gammadim; i++) {
    r.c[i] = lhs.c[i] - rhs.c[i];
  }
  return r;
}









// Multiplication with gamma matrices

#if (Gammadim==4) 

template<typename vector>
Wilson_vector<vector> operator*(const gamma_matrix_type gamma, const Wilson_vector<vector> rhs){
  Wilson_vector<vector>  r;
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

template<typename vector>
Wilson_vector<vector> operator*(const gamma_matrix_type gamma, const Wilson_vector<vector> rhs){
  Wilson_vector<vector>  r;
  switch(gamma) {
    case gamma0:
      r.c[0] = rhs.c[1]; r.c[1] = rhs.c[0];
      break;
    case gamma1:
      r.c[0] = cmplx(0,-1)*rhs.c[1]; r.c[1] = cmplx(0,1)*rhs.c[0];
      break;
    case gamma2:
      r.c[0] = rhs.c[0]; r.c[1] = -rhs.c[1];
      break;
  }
  return r;
}

#endif




/// half_Wilson_vector is a Wilson vector projected by
/// 1 +- gamma_j and contains half the degrees of freedom
/* 
  (1 +- gamma_j) is a projection operator. We will apply the projection
  to a Wilson_vector and store the result in a half_Wilson_vector. This
  will store all necessary information in half the amount of data
  and reduce the effort of multiplying with gauge matrices by half.

  The constructor half_Wilson_vector(Wilson_vector<n, radix> w, direction d)
  will take the projection automatically. A positive direction d corresponds
  to 1+gamma_j and a negative direction d corresponds to 1-gamma_j.

  A half_Wilson_vector can be expanded to a full Wilson_vector using the
  method expand(direction d). Notice that this requires knowing the direction
  that was used to project it. The direction is not stored. 

  Note that the eigenvectors below are normalized to sqrt(2), or |v*v| = 2.
  This is why we don't explicitly multiply by 2 when expanding to full
  Wilson_vector.

 gamma(XUP) 			eigenvectors	eigenvalue
  0  0  0  i		( 1, 0, 0,-i)	  +1
  0  0  i  0		( 0, 1,-i, 0)	  +1
  0 -i  0  0		( 1, 0, 0,+i)	  -1
 -i  0  0  0		( 0, 1,+i, 0)	  -1

 gamma(YUP)			eigenvectors	eigenvalue
  0  0  0 -1		( 1, 0, 0,-1)	  +1
  0  0  1  0		( 0, 1, 1, 0)	  +1
  0  1  0  0		( 1, 0, 0, 1)	  -1
 -1  0  0  0		( 0, 1,-1, 0)	  -1

 gamma(ZUP)			eigenvectors	eigenvalue
  0  0  i  0		( 1, 0,-i, 0)	  +1
  0  0  0 -i		( 0, 1, 0,+i)	  +1
 -i  0  0  0		( 1, 0,+i, 0)	  -1
  0  i  0  0		( 0, 1, 0,-i)	  -1

 gamma(TUP)			eigenvectors	eigenvalue
  0  0  1  0		( 1, 0, 1, 0)	  +1
  0  0  0  1		( 0, 1, 0, 1)	  +1
  1  0  0  0		( 1, 0,-1, 0)	  -1
  0  1  0  0		( 0, 1, 0,-1)	  -1

 gamma(FIVE) 			eigenvectors	eigenvalue
  1  0  0  0    sq2( 1, 0, 0, 0)   +1
  0  1  0  0    sq2( 0, 1, 0, 0)   +1
  0  0 -1  0    sq2( 0, 0, 1, 0)   -1
  0  0  0 -1    sq2( 0, 0, 0, 1)   -1
*/

template<typename vector>
class half_Wilson_vector {
  public:
  vector c[Gammadim/2];
  using base_type = typename base_type_struct<vector>::type;
  
  half_Wilson_vector() = default;


  // This will take the projection 1 +- gamma_j
#if (Gammadim==4) 
  half_Wilson_vector(Wilson_vector<vector> w, direction dir, int sign) {
    switch(dir){
      case XUP:
        if(sign==1){
	        c[0] = w.c[0] + cmplx(0,1)*w.c[3];
	        c[1] = w.c[1] + cmplx(0,1)*w.c[2];
	      } else {
 	        c[0] = w.c[0] - cmplx(0,1)*w.c[3];
	        c[1] = w.c[1] - cmplx(0,1)*w.c[2];
        }
	      break;
      case YUP:
        if(sign==1){
	        c[0] = w.c[0] - w.c[3];
	        c[1] = w.c[1] + w.c[2];
	      } else {
	        c[0] = w.c[0] + w.c[3];
	        c[1] = w.c[1] - w.c[2];
        }
	      break;
      case ZUP:
        if(sign==1){
	        c[0] = w.c[0] + cmplx(0,1)*w.c[2];
	        c[1] = w.c[1] - cmplx(0,1)*w.c[3];
	      } else {
	        c[0] = w.c[0] - cmplx(0,1)*w.c[2];
	        c[1] = w.c[1] + cmplx(0,1)*w.c[3];
        }
	      break;
      case TUP:
        if(sign==1){
	        c[0] = w.c[0] + w.c[2];
	        c[1] = w.c[1] + w.c[3];
	      } else {
	        c[0] = w.c[0] - w.c[2];
	        c[1] = w.c[1] - w.c[3];
        }
	      break;
#if NDIM == 5
      case 4:
        if(sign==1){
  	      c[0] = sqrt(2.0)*w.c[0];
  	      c[1] = sqrt(2.0)*w.c[1];
	      } else {
          c[0] = sqrt(2.0)*w.c[2];
  	      c[1] = sqrt(2.0)*w.c[3];
        }
        break;
#endif
    }
  }

  Wilson_vector<vector> expand(direction dir, int sign) const{
    Wilson_vector<vector> r;
    switch(dir){
      case XUP:
        if(sign==1){
          r.c[0] = c[0]; r.c[1] = c[1];
          r.c[2] = cmplx(0,-1)*c[1];
          r.c[3] = cmplx(0,-1)*c[0];
	      } else {
          r.c[0] = c[0]; r.c[1] = c[1];
          r.c[2] = cmplx(0,1)*c[1];
          r.c[3] = cmplx(0,1)*c[0];
        }
	      break;
      case YUP:
        if(sign==1){
          r.c[0] = c[0]; r.c[1] = c[1];
          r.c[2] = c[1]; r.c[3] = -c[0];
	      } else {
          r.c[0] = c[0];  r.c[1] = c[1];
          r.c[2] = -c[1]; r.c[3] = c[0];
        }
        break;
      case ZUP:
        if(sign==1){
          r.c[0] = c[0]; r.c[1] = c[1];
          r.c[2] = cmplx(0,-1)*c[0];
          r.c[3] = cmplx(0, 1)*c[1];
	      } else {
          r.c[0] = c[0]; r.c[1] = c[1];
          r.c[2] = cmplx(0, 1)*c[0];
          r.c[3] = cmplx(0,-1)*c[1];
        }
        break;
      case TUP:
        if(sign==1){
          r.c[0] = c[0]; r.c[1] = c[1];
          r.c[2] = c[0]; r.c[3] = c[1];
	      } else {
          r.c[0] = c[0];  r.c[1] = c[1];
          r.c[2] = -c[0]; r.c[3] = -c[1];
        }
	      break;
#if NDIM == 5
      case 4:
        if(sign==1){
          r.c[0] = sqrt(2.0)*c[0]; r.c[1] = sqrt(2.0)*c[1];
          r.c[2] = 0; r.c[3] = 0;
	      } else {
          r.c[0] = 0; r.c[1] = 0;
          r.c[2] = sqrt(2.0)*c[0]; r.c[3] = sqrt(2.0)*c[1];
        }
        break;
#endif
    }
    return r;
  }

#elif (Gammadim==2)
/*
 gamma(XUP) 	 eigenvectors	 eigenvalue
   0  1		      ( 1, 1)	       +1
   1  0		      ( 1,-1)	       -1

 gamma(YUP)		 eigenvectors	 eigenvalue
   0  i	        ( 1, i)	       +1
  -i  0	  	    ( 1,-i)	       -1

 gamma(ZUP)		 eigenvectors  eigenvalue
   1  0	        ( 1, 0)	       +1
   0 -1	  	    ( 0, 1)	       -1
*/
  half_Wilson_vector(Wilson_vector<vector> w, direction dir, int sign) {
    switch(dir){
      case XUP:
        if(sign==1){
	        c[0] = w.c[0] + w.c[1];
	      } else {
	        c[0] = w.c[0] - w.c[1];
        }
	      break;
      case YUP:
        if(sign==1){
	        c[0] = w.c[0] - cmplx(0,1)*w.c[1];
	      } else {
 	        c[0] = w.c[0] + cmplx(0,1)*w.c[1];
        }
	      break;
#if NDIM == 3
      case ZUP:
        if(sign==1){
  	      c[0] = sqrt(2.0)*w.c[0];
	      } else {
          c[0] = sqrt(2.0)*w.c[1];
        }
        break;
#endif
    }
  }

  Wilson_vector<vector> expand(direction dir, int sign) const{
    Wilson_vector<vector> r;
    switch(dir){
      case XUP:
        if(sign==1){
          r.c[0] = c[0]; r.c[1] = c[0];
	      } else {
          r.c[0] = c[0]; r.c[1] = -c[0];
        }
	      break;
      case YUP:
          r.c[0] = c[0]; r.c[1] = cmplx(0,1)*c[0];
	      } else {
          r.c[0] = c[0]; r.c[1] = cmplx(0,-1)*c[0];
        }
        break;
#if NDIM == 3
      case ZUP:
          r.c[0] = sqrt(2.0)*c[0]; r.c[1] = 0;
	      } else {
          r.c[0] = 0; r.c[1] = sqrt(2.0)*c[0];
        }
        break;
#endif
    }
    return r;
  }

#endif


  /// Returns the norm squared of (1+-gamma_j) * wilson_vector.
  /// Thus the factor 2.
  inline auto norm_sq(){ 
    auto r=c[0].norm_sq();
    for (int i = 1; i < Gammadim; i++) {
      r += c[i].norm_sq();
    }
    return r;
  }

  #pragma transformer loop_function
  half_Wilson_vector & operator+=(const half_Wilson_vector & rhs){
    for (int i = 0; i < Gammadim/2; i++){
      c[i] += rhs.c[i];
    }
    return *this;
  }

  #pragma transformer loop_function
  half_Wilson_vector & operator-=(const half_Wilson_vector & rhs){
    for (int i = 0; i < Gammadim/2; i++){
      c[i] -= rhs.c[i];
    }
    return *this;
  }

  half_Wilson_vector operator-() const {
    half_Wilson_vector r;
    for (int i = 0; i < Gammadim/2; i++){
      r.c[i] = -c[i];
    }
    return r;
  }

  std::string str() const {
    std::string text = "";
    for (int i=0; i<Gammadim/2; i++){
      text += c[i].str() + "\n"; 
    }
    text += "\n";
    return text;
  }
};





template<typename vector, typename T>
half_Wilson_vector<vector> operator*(const T lhs, const half_Wilson_vector<vector> rhs){
  half_Wilson_vector<vector> r;
  for (int i=0; i<Gammadim/2; i++) {
    r.c[i] = lhs*rhs.c[i];
  }
  return r;
}

template<typename vector, typename T>
half_Wilson_vector<vector> operator*(const half_Wilson_vector<vector> lhs, const T rhs){
  half_Wilson_vector<vector> r;
  for (int i=0; i<Gammadim/2; i++) {
    r.c[i] = lhs.c[i]*rhs;
  }
  return r;
}


template<typename vector>
half_Wilson_vector<vector> operator+(const half_Wilson_vector<vector> lhs, const half_Wilson_vector<vector> rhs){
  half_Wilson_vector<vector>  r;
  for (int i=0; i<Gammadim/2; i++) {
    r.c[i] = lhs.c[i] + rhs.c[i];
  }
  return r;
}

template<typename vector>
half_Wilson_vector<vector> operator-(const half_Wilson_vector<vector> lhs, const half_Wilson_vector<vector> rhs){
  half_Wilson_vector<vector>  r;
  for (int i=0; i<Gammadim/2; i++) {
    r.c[i] = lhs.c[i] - rhs.c[i];
  }
  return r;
}









#endif