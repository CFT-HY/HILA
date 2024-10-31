/**
 * @file bench.cpp
 * @author Kari Rummukainen
 * @brief simple benchmark application, measures performance of various operations
 */
#include "hila.h"

// unistd.h needed for isatty()
#include <unistd.h>


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// benchmark rnd generation

void bench_random() {

    hila::out0 << "\n-------------------------------------\n";
    hila::out0 << "Random number generation: ";

    Field<double> f;

    // warm it up
    f[ALL] = hila::random();

    // estimate loop number
    auto time = hila::gettime();
    for (int i = 0; i < 10; ++i) {
        f[ALL] = hila::random();
    }
    hila::synchronize();
    time = hila::gettime() - time;

    int n_loops = 50.0 / time; // gives 5s
    hila::out0 << n_loops << " iterations\n";

    time = hila::gettime();
    for (int i = 0; i < n_loops; ++i) {
        f[ALL] = hila::random();
    }
    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  In separate onsites loops: " << time << " seconds, " << time / n_loops
               << " per loop, " << time / n_loops / lattice.volume() << " per site\n";

    time = hila::gettime();

    onsites(ALL) {
        for (int i = 0; i < n_loops; ++i) {
            f[X] = hila::random();
        }
    }

    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  In a single onsites loop: " << time << " seconds, " << time / n_loops
               << " per loop, " << time / n_loops / lattice.volume() << " per site\n";
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Benchmark 3x3 matrix multiply

void bench_matrix() {

    hila::out0 << "\n-------------------------------------\n";
    hila::out0 << "3x3 Complex double matrix multiply: ";

    Field<Matrix<3, 3, Complex<double>>> f, g;

    // warm it up
    f = 1;
    g = 1;

    // estimate loop number
    auto time = hila::gettime();
    for (int i = 0; i < 10; ++i) {
        g[ALL] = f[X] * f[X];
    }
    hila::synchronize();
    time = hila::gettime() - time;

    int n_loops = 50.0 / time; // gives 5s
    hila::out0 << n_loops << " iterations\n";

    time = hila::gettime();
    for (int i = 0; i < n_loops; ++i) {
        f[ALL] = f[X] * f[X];
    }
    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  In separate onsites loops: " << time << " seconds, " << time / n_loops
               << " per loop, " << time / n_loops / lattice.volume() << " per site\n";

    time = hila::gettime();
    onsites(ALL) {
        for (int i = 0; i < n_loops; ++i) {
            f[X] = f[X] * f[X];
        }
    }
    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  In a single onsites loop: " << time << " seconds, " << time / n_loops
               << " per loop, " << time / n_loops / lattice.volume() << " per site\n";
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Benchmark mpi comms

void bench_comm() {

    constexpr int n_gathers = 300;

    hila::out0 << "\n-------------------------------------\n";
    hila::out0 << "Nearest neighbour communication: complex field\n";

    Field<Complex<double>> df = 0;

    // warm up a bit
    foralldir(d) {
        for (int i = 0; i < 3; i++) {
            df.gather(d);
            df.gather(-d);
            df.mark_changed(ALL);
        }
    }

    for (Direction d = e_x; d < NDIRECTIONS; ++d) {
        df.gather(d);
        df.mark_changed(ALL);

        auto time = hila::gettime();

        for (int i = 0; i < n_gathers; i++) {
            df.gather(d);
            df.mark_changed(ALL);
        }
        hila::synchronize();
        time = hila::gettime() - time;
        hila::out0 << "  Gather from direction " << hila::prettyprint(d) << ": " << time / n_gathers
                   << " s/gather\n";
    }
}

//--------------------------------------------------------------------------------

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Benchmark FFT

void bench_fft() {

    constexpr int n_fft = 20;

    hila::out0 << "\n-------------------------------------\n";
    hila::out0 << "FFT of a double complex field\n";

    Field<Complex<double>> df = 0, rf;

    df.gaussian_random();

    rf = df.FFT();

    auto time = hila::gettime();

    for (int i = 0; i < n_fft; ++i) {
        rf = df.FFT();
    }
    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  " << n_fft << " FFTs " << time << " sec, " << time / n_fft
               << " for single FFT\n";
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Benchmark simple smear-update

void bench_update() {

    constexpr int n_update = 500;

    hila::out0 << "\n-------------------------------------\n";
    hila::out0 << "NN-smear complex field " << n_update << " times\n";

    Field<Complex<double>> df, rf;

    df.gaussian_random();

    onsites(ALL) {
        rf[X] = df[X];
        foralldir(d) rf[X] += df[X + d] + df[X - d];
    }

    auto time = hila::gettime();
    for (int i = 0; i < n_update; i++) {
        onsites(ALL) {
            rf[X] = df[X];
            foralldir(d) rf[X] += df[X + d] + df[X - d];
        }
        df = rf;
    }
    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  Total time " << time << " s, one update " << time / n_update << ", per site "
               << time / n_update / lattice.volume() << '\n';
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Benchmark simple smear-update

void bench_matrix_update() {

    constexpr int n_update = 50;

    hila::out0 << "\n-------------------------------------\n";
    hila::out0 << "NN-mult SU(5) matrix field " << n_update << " times\n";

    Field<SU<5, double>> df, rf;

    df = 1;
    rf = 1;

    foralldir(d) onsites(ALL) {
        rf[X] = rf[X] * df[X + d] * df[X - d];
    }

    auto time = hila::gettime();
    for (int i = 0; i < n_update; i++) {
        onsites(ALL) {
            foralldir(d) rf[X] = df[X + d] * df[X - d];
        }
    }
    hila::synchronize();
    time = hila::gettime() - time;

    hila::out0 << "  Total time " << time << " s, one update " << time / n_update << ", per site "
               << time / n_update / lattice.volume() << '\n';
}


//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

    hila::initialize(argc, argv);

    hila::out0 << "HILA benchmark program\n";

    hila::input par("parameters");

    CoordinateVector lsize = par.get("lattice size"); // reads NDIM numbers
    long seed = par.get("random seed");

    par.close();

    // setting up the lattice is convenient to do after reading
    // the parameters
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);


    ///////////////////////////////////////////////////////////////
    // start tests

    hila::out0 << "###################################\n\n";

    bench_random();

    bench_matrix();

    bench_comm();

    bench_fft();

    bench_update();

    bench_matrix_update();

    hila::out0 << "\n##################################\n";

    hila::finishrun();
}
