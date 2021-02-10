#ifndef MRE_GUESS_H
#define MRE_GUESS_H

#include "gauge_field.h"
#include "dirac/Hasenbusch.h"
#include <cmath>

/// Builds an initial guess for a matrix inverter given a set of basis vectors
template <typename vector_type, typename DIRAC_OP>
void MRE_guess(Field<vector_type> &psi, Field<vector_type> &chi, DIRAC_OP D,
               std::vector<Field<vector_type>> old_chi_inv) {
    int MRE_size = old_chi_inv.size();
    double M[MRE_size][MRE_size];
    double v[MRE_size];
    Field<vector_type> basis[MRE_size];
    Field<vector_type> tmp;

    // Build an orthogonal basis from the previous solutions
    for (int i = 0; i < MRE_size; i++) {
        // Start with the original solution vector
        basis[i][ALL] = old_chi_inv[i][X];
        // Remove the projected components of all previous vectors
        for (int j = 0; j < i; j++) {
            double vdot = 0, norm = 0;
            onsites(D.par) {
                norm += basis[i][X].rdot(basis[i][X]);
                vdot += basis[j][X].rdot(basis[i][X]);
            }
            if (norm * norm > 1e-32) {
                onsites(D.par) { basis[i][X] = basis[i][X] - vdot / norm * basis[j][X]; }
            }
        }
    }

    // Build the projected matrix, M[i][j] = v[i].v[j]
    for (int i = 0; i < MRE_size; i++) {
        Field<vector_type> Dchi, DDchi;
        D.apply(basis[i], Dchi);
        D.dagger(Dchi, DDchi);
        for (int j = 0; j < MRE_size; j++) {
            double sum = 0;
            onsites(D.par) { sum += basis[j][X].rdot(DDchi[X]); }
            M[j][i] = sum;
        }
    }
    // And the projected vector
    for (int i = 0; i < MRE_size; i++) {
        double sum = 0;
        onsites(D.par) { sum += basis[i][X].rdot(chi[X]); }
        v[i] = sum;
    }

    // Now invert the small matrix M (Gaussian elimination)
    for (int i = 0; i < MRE_size; i++) {
        // Check that the diagonal element is nonzero
        if (M[i][i] * M[i][i] > 1e-32) {
            // Normalize the i:th row
            double diag_inv = 1.0 / M[i][i];
            for (int j = 0; j < MRE_size; j++) {
                M[i][j] *= diag_inv;
            }
            v[i] *= diag_inv;
            // Subtract from all other rows
            for (int k = 0; k < MRE_size; k++)
                if (k != i) {
                    double weight = M[k][i];
                    for (int j = 0; j < MRE_size; j++) {
                        M[k][j] -= weight * M[i][j];
                    }
                    v[k] -= weight * v[i];
                }
        } else {
            // In the matrix element is too small, just remove it from the basis
            v[i] = 0;
            for (int j = 0; j < MRE_size; j++) {
                M[j][i] = 0;
            }
        }
    }

    // Construct the solution in the original basis
    psi[ALL] = 0;
    for (int i = 0; i < MRE_size; i++)
        if (!isnan(v[i])) {
            onsites(D.par) { psi[X] = basis[i][X] + v[i] * basis[i][X]; }
        }
}

#endif