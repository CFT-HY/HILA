#ifndef __HASENBUSCH_H__
#define __HASENBUSCH_H__

#include "staggered.h"
#include "wilson.h"

/* The Hasenbusch method for updating fermion fields:
 * Split the Dirac determinant into two parts, 
 * D_h1 = D + mh and
 * D_h2 = D * 1 / (D + mh)^dagger
 */

template<typename Dirac_type>
class Hasenbusch_operator {
  /* An operator that applies D + hm. This is 
   * necessary for defining the inverse 1/(D+hm) */
  public:
    using vector_type = typename Dirac_type::vector_type;
    Dirac_type D;
    parity par = D.par;
    double h_parameter;
    Hasenbusch_operator(Dirac_type &d, double hasenbusch_parameter) : D(d) {
      h_parameter = hasenbusch_parameter;
    }
    Hasenbusch_operator(Hasenbusch_operator &h) : D(h.D) {
      h_parameter = h.h_parameter;
    }

    inline void apply( const field<vector_type> & in, field<vector_type> & out){
      D.apply(in, out);
      out[D.par] = out[X] + h_parameter*in[X];
    }

    inline void dagger( const field<vector_type> & in, field<vector_type> & out){
      D.dagger(in, out);
      out[D.par] = out[X] + h_parameter*in[X];
    }

    template<typename momtype>
    inline void force(const field<vector_type> & chi, const field<vector_type> & psi, field<momtype> (&force)[NDIM], int sign){
      D.force(chi, psi, force, sign);
    }
};


#endif