#ifndef STOUT_SMEAR_H_
#define STOUT_SMEAR_H_

#include "hila.h"

#include "gauge/staples.h"

/////////////////////////////////////////////////////////////////////////////
/// Do stout smearing for the gauge fields
///  res = U * exp( Algebra (U * Staples) )

template <typename T>
void stout_smear1(const GaugeField<T> &U, Field<T> &s, float coeff, Direction d) {
    staplesum(U, s, d, ALL);
    onsites(ALL) {
        s[X] = exp(coeff * (s[X] * U[d][X].dagger()).project_to_algebra()) * U[d][X];
    }
}

template <typename T>
void stout_smear1(const GaugeField<T> &U, GaugeField<T> &stout, float coeff) {
    foralldir(d) stout_smear1(U,stout[d],coeff,d);
    // stout_smear1(U,stout[e_t],coeff,e_t);
}



template <typename T>
void stout_smear(const GaugeField<T> &U, GaugeField<T> &stout, float coeff, int iter) {
    stout = U;
    if (iter > 0) {
        stout_smear1(U,stout,coeff);
    
        GaugeField<T> tmp;
        for (int i=1; i<iter; i++) {
            stout_smear1(stout,tmp,coeff);
            stout = tmp;
//            std::swap(stout,tmp);
        }
    }
}



#endif