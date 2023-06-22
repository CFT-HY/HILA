/**
 * @file hila_simple_example.cpp
 * @author Aaron Haarti
 * @brief Simple HILA example application which computes laplacian of Gaussian random complex
 * field
 * @details The application generates a \f$ 32^2 \f$ size field \f$ f \f$ of random Complex numbers.
 * We proceed to compute the summed absolute value of the Laplacian of said field:
 *
 * \f{align}{g(X) &= |\nabla^2 f(X)| \\
 *       &= \sum_{d \in \hat{e}} |f(X + d) - 2f(X) + f(X-d)|, \f}
 *
 * where \f$\hat{e} = \{e_x,e_y,e_z\}\f$ is the set of unit vectors. Afterwards this quantity is
 * printed out.
 *
 */

#include "hila.h"
static_assert(NDIM == 3, "NDIM must be 3");

/**
 * @brief
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char *argv[]) {

    hila::initialize(argc, argv);

    // set up 32^3 lattice
    lattice.setup({32, 32, 32});

    // Random numbers are used here
    hila::seed_random(32345);

    Field<Complex<double>> f;
    Field<double> g = 0;

    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian_random();

    // calculate sum of 2nd derivatives of f in to g
    foralldir(d) {
        g[ALL] += abs(f[X + d] - 2 * f[X] + f[X - d]);
    }

    // get average of g
    double average = 0;
    onsites(ALL) {
        average += g[X];
    }

    average = average / lattice.volume();
    hila::out0 << "Average of g is " << average << '\n';

    // make a clean exit
    hila::finishrun();
}