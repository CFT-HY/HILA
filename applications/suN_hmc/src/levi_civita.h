/** @file levi_civita.h */

#ifndef LEVI_CIVITA_H_
#define LEVI_CIVITA_H_

#include "hila.h"

/**
 * @brief class providing non-zero entries of N-dimensional 
 * Levi-Civita symbol in sparse form.
 */

template <int N>
struct Factorial {
    enum { value=N*Factorial<N-1>::value };
};

template <>
struct Factorial<0> {
    enum { value=1 };
};


template<int N,int Nfac=Factorial<N>::value>
class levi_civita {
    // template class to provide the non-zero elements of the N-dimensional 
    // Levi-Civita symbol.
    // (would be nice to make a[Nfac][N+1] static constexpr, but haven't
    //  figured out yet how to do it with C++17) 
public:
    static constexpr int n=N;
    static constexpr int nfac=Nfac;

    int a[Nfac][N+1];

    constexpr levi_civita() {
        int tarr0[N+1]; // inital sequence and permutation sign :
        // sequence {0,1,...} :
        for(int i=0; i<N; ++i) {
            tarr0[i]=i;
        }
        // sign 1 :
        tarr0[N]=1;

        int ind[N]{0}; //initial permutation vector is zero
        for(int j=0; j<Nfac; ++j) {
            // initialize a[j][]=tarr0[] :
            for(int i=0; i<n+1; ++i) {
                a[j][i]=tarr0[i];
            }
            // apply permutation given by ind[] to a[j][] :
            for(int i=0; i<n-1; ++i) {
                if(ind[i]!=0) {
                    // swap a[j][i] with a[j][i+ind[i]]:
                    int te=a[j][i];
                    int ti=i+ind[i];
                    a[j][i]=a[j][ti];
                    a[j][ti]=te;
                    // flip sign on a[j][]: 
                    a[j][N]=-a[j][N];
                }
            }

            // next permuation vector: 
            ++ind[N-1];
            for(int i=n-1; i>=1; --i) {
                if(ind[i]>n-1-i) {
                    ind[i]=0;
                    ++ind[i-1];
                }
            }
        }

    }

};
