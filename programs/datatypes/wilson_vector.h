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
Wilson_vector<n,radix> operator*(const gamma_matrix_type gamma, const Wilson_vector<n,radix> rhs){
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

template<int n, typename radix=double>
class half_Wilson_vector {
  public:
  SU_vector<n, radix> c[Gammadim/2];
  using base_type = typename base_type_struct<radix>::type;
  
  half_Wilson_vector() = default;


  // This will take the projection 1 +- gamma_j
#if (Gammadim==4) 
  half_Wilson_vector(Wilson_vector<n, radix> w, direction dir) {
    switch(dir){
      case XUP:
	      c[0] = w.c[0] + cmplx(0,1)*w.c[3];
	      c[1] = w.c[1] + cmplx(0,1)*w.c[2];
	      break;
      case XDOWN:
 	      c[0] = w.c[0] - cmplx(0,1)*w.c[3];
	      c[1] = w.c[1] - cmplx(0,1)*w.c[2];
	      break;
      case YUP:
	      c[0] = w.c[0] - w.c[3];
	      c[1] = w.c[1] + w.c[2];
	      break;
      case YDOWN:
	      c[0] = w.c[0] + w.c[3];
	      c[1] = w.c[1] - w.c[2];
	      break;
      case ZUP:
	      c[0] = w.c[0] + cmplx(0,1)*w.c[2];
	      c[1] = w.c[1] - cmplx(0,1)*w.c[3];
	      break;
      case ZDOWN:
	      c[0] = w.c[0] - cmplx(0,1)*w.c[2];
	      c[1] = w.c[1] + cmplx(0,1)*w.c[3];
	      break;
      case TUP:
	      c[0] = w.c[0] + w.c[2];
	      c[1] = w.c[1] + w.c[3];
	      break;
      case TDOWN:
	      c[0] = w.c[0] - w.c[2];
	      c[1] = w.c[1] - w.c[3];
	      break;
#if NDIM == 5
      case 4:
  	    c[0] = sqrt(2.0)*w.c[0];
  	    c[1] = sqrt(2.0)*w.c[1];
        break;
      case 5:
        c[0] = sqrt(2.0);*w.c[2];
  	    c[1] = sqrt(2.0);*w.c[3];
        break;
#endif
    }
  }

  Wilson_vector<n, radix> expand(direction dir){
    Wilson_vector<n, radix> r;
    switch(dir){
      case XUP:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = cmplx(0,-1)*c[1];
        r.c[3] = cmplx(0,-1)*c[0];
	      break;
      case XDOWN:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = cmplx(0,1)*c[1];
        r.c[3] = cmplx(0,1)*c[0];
	      break;
      case YUP:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = c[1];
        r.c[3] = -c[0];
	      break;
      case YDOWN:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = -c[1];
        r.c[3] = c[0];
        break;
      case ZUP:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = cmplx(0,-1)*c[0];
        r.c[3] = cmplx(0, 1)*c[1];
	      break;
      case ZDOWN:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = cmplx(0, 1)*c[0];
        r.c[3] = cmplx(0,-1)*c[1];
        break;
      case TUP:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = c[0];
        r.c[3] = c[1];
	      break;
      case TDOWN:
        r.c[0] = c[0]; r.c[1] = c[1];
        r.c[2] = -c[0];
        r.c[3] = -c[1];
	      break;
#if NDIM == 5
      case 4:
        r.c[0] = sqrt(2.0)*c[0]; r.c[1] = sqrt(2.0)*c[1];
        r.c[2] = 0; r.c[3] = 0;
        break;
      case 5:
        r.c[0] = 0; r.c[1] = 0;
        r.c[2] = sqrt(2.0)*c[0]; r.c[3] = sqrt(2.0)*c[1];
        break;
#endif
    }
    return r;
  }

#elif (Gammadim==2)


#endif


  /// Returns the norm squared of (1+-gamma_j) * wilson_vector.
  /// Thus the factor 2.
  radix norm_sq(){ 
    radix r=0;
    for (int i = 0; i < Gammadim/2; i++) {
      r += c[i].norm_sq();
    }
    return 2*r;
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

















#endif