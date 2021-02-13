#include "SUN.h"

// Define some parameters for the simulation
double beta = 8;
int n_measurements = 10;
int n_updates_per_measurement = 10;
long seed = 123456;
int NX = 16, NY = 16, NZ = 16, NT = 16;

void calc_staples(Field<Matrix<N, N, Cmplx<double>>> (&U)[NDIM],
                  Field<Matrix<N, N, Cmplx<double>>> &staple_sum, direction dir) {
    /* Calculate the sum of staples connected to links in direction
     * dir
     */
    static Field<Matrix<N, N, Cmplx<double>>> down_staple;
    staple_sum[ALL] = 0;
    foralldir(d2) {
        direction dir2 = (direction)d2;
        // Calculate the down side staple.
        // This will be communicated up.
        down_staple[ALL] = U[dir2][X].conjugate() * U[dir][X] * U[dir2][X + dir];
        // Forward staple
        staple_sum[ALL] +=
            U[dir2][X + dir] * U[dir][X + dir2].conjugate() * U[dir2][X].conjugate();
        // Add the two staples together
        staple_sum[ALL] += down_staple[X - dir2];
    }
}

template <typename T> void update(T &U, const T &staple, double beta) {
    monte(U, staple, beta);
}

void update(Matrix<2, 2, Cmplx<double>> &U, const Matrix<2, 2, Cmplx<double>> &staple,
            double beta) {
    Matrix<2, 2, Cmplx<double>> temp = -beta * staple;
    KennedyPendleton(U, temp);
}

int main(int argc, char **argv) {
    hila::initialize(argc, argv);
    const CoordinateVector nd{NX, NY, NZ, NT};
    lattice->setup(nd);

    // Define a field
    Field<Matrix<N, N, Cmplx<double>>> U[NDIM];
    Field<Matrix<N, N, Cmplx<double>>> staple;

    seed_random(seed);

    // Set to 1
    foralldir(d) { U[d][ALL] = 1; }

    // Run update-measure loop
    for (int i = 0; i < n_measurements; i++) {

        // Run a number of updates
        for (int j = 0; j < n_updates_per_measurement; j++) {
            foralldir(d) {
                // update direction dir
                direction dir = (direction)d;
                // First we need the staple sum
                calc_staples(U, staple, dir);

                // Now update, first even then odd
                parity p = EVEN;
                for (int par = 0; par < 2; par++) {
                    onsites(p) {
                        if (disable_avx[X] == 0) {
                        };
                        update(U[dir][X], staple[X], beta);
                    }
                    p = opp_parity(p);
                }
            }
        }

        // Measure plaquette
        double Plaq = 0;
        foralldir(d1) foralldir(d2) if (d1 != d2) {
            direction dir1 = (direction)d1, dir2 = (direction)d2;
            onsites(ALL) {
                element<Matrix<N, N, Cmplx<double>>> temp;
                temp = U[dir1][X] * U[dir2][X + dir1];
                temp *= U[dir1][X + dir2].conjugate();
                temp *= U[dir2][X].conjugate();
                Plaq += 1 - temp.trace().re / N;
            }
        }
        output0 << "Plaquette " << Plaq / (lattice->volume() * NDIM * (NDIM - 1)) << "\n";
    }

    hila::finishrun();
    return 0;
}
