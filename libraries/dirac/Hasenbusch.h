#ifndef __HASENBUSCH_H__
#define __HASENBUSCH_H__

#include "staggered.h"
#include "wilson.h"

/* The Hasenbusch method for updating fermion fields:
 * Split the Dirac determinant into two parts, 
 * D_h1 = D + mh and
 * D_h2 = D * 1 / (D + mh)
 */

/// This operator applies D + mh. It is 
/// necessary for defining the inverse 1/(D+mh)
template<typename Dirac_type>
class Hasenbusch_operator {
  public:
    using vector_type = typename Dirac_type::vector_type;
    using type_flt = Hasenbusch_operator<typename Dirac_type::type_flt>;

    Dirac_type D;
    parity par;
    double h_parameter;
    Hasenbusch_operator(Dirac_type &d, double hasenbusch_parameter) : D(d) {
      h_parameter = hasenbusch_parameter;
      par = D.par;
    }
    Hasenbusch_operator(Hasenbusch_operator &h) : D(h.D) {
      h_parameter = h.h_parameter;
      par = D.par;
    }

    template<typename H_type, typename matrix>
    Hasenbusch_operator(H_type &h, matrix &mat) : D(h.D, mat) {
      h_parameter = h.h_parameter;
      par = D.par;
    }


    inline void apply( const Field<vector_type> & in, Field<vector_type> & out){
      D.apply(in, out);
      out[D.par] = out[X] + h_parameter*in[X];
    }

    inline void dagger( const Field<vector_type> & in, Field<vector_type> & out){
      D.dagger(in, out);
      out[D.par] = out[X] + h_parameter*in[X];
    }

    template<typename momtype>
    inline void force(const Field<vector_type> & chi, const Field<vector_type> & psi, Field<momtype> (&force)[NDIM], int sign){
      D.force(chi, psi, force, sign);
    }
};


#endif