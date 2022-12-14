#include "test.h"

#include "dirac/staggered.h"
#include "dirac/wilson.h"
#include "dirac/Hasenbusch.h"
#include "dirac/conjugate_gradient.h"

#define N 3

void test_gamma_matrices() {
    Wilson_vector<N, double> w1, w2, w3;
    SU<N> U;
    U.random();
    w1.gaussian_random();

#if NDIM == 4
    w2 = w1 - gamma5 * (gamma5 * w1);
    assert(w2.squarenorm() < 0.0001 && "g5*g5 = 1");
#endif

    foralldir(d) {
        half_Wilson_vector<N, double> h1;
        w2 = w1 - gamma_matrix[d] * (gamma_matrix[d] * w1);
        assert(w2.squarenorm() < 0.0001 && "gamma_d*gamma_d = 1");

        w2 = w1 + gamma_matrix[d] * w1;
        h1 = half_Wilson_vector(w1, d, 1);
        double diff = w2.squarenorm() - 2 * h1.squarenorm();
        assert(diff * diff < 0.0001 && "half_Wilson_vector projection +1 norm");

        w3 = (U * h1).expand(d, 1) - U * w2;
        assert(w3.squarenorm() < 0.0001 && "half_wilson_vector expand");

        w2 = w1 - gamma_matrix[d] * w1;
        h1 = half_Wilson_vector(w1, d, -1);
        diff = w2.squarenorm() - 2 * h1.squarenorm();
        assert(diff * diff < 0.0001 && "half_Wilson_vector projection -1 norm");

        w3 = (U * h1).expand(d, -1) - U * w2;
        assert(w3.squarenorm() < 0.0001 && "half_wilson_vector expand");
    }
}

int main(int argc, char **argv) {

#if NDIM == 1
    const CoordinateVector nd = {64};
#elif NDIM == 2
    const CoordinateVector nd = {32, 8};
#elif NDIM == 3
    const CoordinateVector nd = {16, 8, 8};
#elif NDIM == 4
    const CoordinateVector nd = {16, 8, 8, 8};
#endif
    hila::initialize(argc, argv);
    lattice.setup(nd);

    hila::seed_random(2);

    test_gamma_matrices();

    Field<SU<N>> U[NDIM];
    foralldir(d){onsites(ALL){U[d][X] = 1;
}
}

// Check conjugate of the staggered Dirac operator
{
    hila::out0 << "Checking with dirac_staggered\n";
    using dirac = dirac_staggered<SU<N>>;
    dirac D(0.1, U);
    Field<SU_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    onsites(ALL) {
        a[X].gaussian_random();
        b[X].gaussian_random();
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(ALL) {
        diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
        diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }

    assert(diffre * diffre < 1e-16 && "test dirac_staggered");
    assert(diffim * diffim < 1e-16 && "test dirac_staggered");

    // Now run CG on Ddaggera. b=1/D a -> Db = a
    CG<dirac> inverse(D);
    onsites(ALL) {
        a[X].gaussian_random();
    }
    b[ALL] = 0;
    D.dagger(a, Ddaggera);
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    diffre = 0;
    onsites(ALL) { diffre += squarenorm(a[X] - Db[X]); }
    assert(diffre * diffre < 1e-16 && "test D (DdgD)^-1 Ddg");

    // Run CG on a, multiply by DdaggerD and check the same
    inverse.apply(a, b);
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);

    diffre = 0;
    onsites(ALL) { diffre += squarenorm(a[X] - DdaggerDb[X]); }
    assert(diffre * diffre < 1e-16 && "test DdgD (DdgD)^-1");

    // The other way around, multiply first
    onsites(ALL) {
        b[X].gaussian_random();
    }
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);
    a[ALL] = 0;
    inverse.apply(DdaggerDb, a);

    diffre = 0;
    onsites(ALL) { diffre += squarenorm(a[X] - b[X]); }
    assert(diffre * diffre < 1e-16 && "test (DdgD)^-1 DdgD");
}

// Check conjugate of the wilson Dirac operator
{
    hila::out0 << "Checking with Dirac_Wilson\n";
    using dirac = Dirac_Wilson<SU<N, double>>;
    dirac D(0.05, U);
    Field<Wilson_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
#if NDIM > 3
    a.set_boundary_condition(e_t, hila::bc::ANTIPERIODIC);
    b.copy_boundary_condition(a);
    Db.copy_boundary_condition(a);
    Ddaggera.copy_boundary_condition(a);
    DdaggerDb.copy_boundary_condition(a);
#endif

    onsites(ALL) {
        a[X].gaussian_random();
        b[X].gaussian_random();
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(ALL) {
        diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
        diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }

    assert(diffre * diffre < 1e-16 && "test Dirac_Wilson");
    assert(diffim * diffim < 1e-16 && "test Dirac_Wilson");

    // Now run CG on Ddaggera. b=1/D a -> Db = a
    CG<dirac> inverse(D);
    b[ALL] = 0;
    onsites(ALL) {
        a[X].gaussian_random();
    }
    inverse.apply(a, b);
    D.dagger(b, Db);
    D.apply(Db, DdaggerDb);

    diffre = 0;
    onsites(ALL) { diffre += squarenorm(a[X] - DdaggerDb[X]); }
    assert(diffre * diffre < 1e-8 && "test DdgD (DdgD)^-1");

    // Now run CG on Ddaggera. b=1/D a -> Db = a
    b[ALL] = 0;
    onsites(ALL) {
        a[X].gaussian_random();
    }
    D.dagger(a, Ddaggera);
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    diffre = 0;
    onsites(ALL) { diffre += squarenorm(a[X] - Db[X]); }
    assert(diffre * diffre < 1e-8 && "test D (DdgD)^-1 Ddg");
}

// Check conjugate of the even-odd preconditioned staggered Dirac operator
{
    hila::out0 << "Checking with dirac_staggered_evenodd\n";
    dirac_staggered_evenodd D(5.0, U);
    Field<SU_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    onsites(ALL) {
        a[X].gaussian_random();
        b[X].gaussian_random();
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(EVEN) {
        diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
        diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }

    assert(diffre * diffre < 1e-16 && "test dirac_staggered_evenodd");
    assert(diffim * diffim < 1e-16 && "test dirac_staggered_evenodd");

    // Now run CG on Ddaggera. b=1/D a -> Db = a
    CG inverse(D);
    b[ALL] = 0;
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    onsites(EVEN) { diffre += squarenorm(a[X] - Db[X]); }
    assert(diffre * diffre < 1e-8 && "test CG");
}

// Check conjugate of the even-odd preconditioned wilson Dirac operator
{
    hila::out0 << "Checking with Dirac_Wilson_evenodd\n";
    using dirac = Dirac_Wilson_evenodd<SU<N>>;
    dirac D(0.12, U);
    Field<Wilson_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
#if NDIM > 3
    a.set_boundary_condition(e_t, hila::bc::ANTIPERIODIC);
    b.set_boundary_condition(e_t, hila::bc::ANTIPERIODIC);
    Db.set_boundary_condition(e_t, hila::bc::ANTIPERIODIC);
    Ddaggera.set_boundary_condition(e_t, hila::bc::ANTIPERIODIC);
    DdaggerDb.set_boundary_condition(e_t, hila::bc::ANTIPERIODIC);
#endif

    a[ODD] = 0;
    b[ODD] = 0;
    onsites(EVEN) {
        a[X].gaussian_random();
        b[X].gaussian_random();
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(EVEN) {
        diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
        diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }

    assert(diffre * diffre < 1e-16 && "test Dirac_Wilson_evenodd");
    assert(diffim * diffim < 1e-16 && "test Dirac_Wilson_evenodd");

    // Now check that D (DdgD)^-1 Ddg = 1
    CG<dirac> inverse(D);
    a[ALL] = 0;
    b[ALL] = 0;
    onsites(EVEN) {
        a[X].gaussian_random();
    }
    D.dagger(a, Ddaggera);
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    diffre = 0;
    onsites(EVEN) { diffre += squarenorm(a[X] - Db[X]); }
    hila::out0 << "D inv Dg " << diffre << "\n";
    assert(diffre * diffre < 1e-16 && "test D (DdgD)^-1 Ddg");

    // Run CG on a, multiply by DdaggerD and check the same
    onsites(EVEN) {
        a[X].gaussian_random();
    }
    b[ALL] = 0;
    inverse.apply(a, b);
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);

    diffre = 0;
    onsites(EVEN) { diffre += squarenorm(a[X] - DdaggerDb[X]); }
    assert(diffre * diffre < 1e-16 && "test DdgD (DdgD)^-1");

    // The other way around, multiply first
    onsites(EVEN) {
        b[X].gaussian_random();
    }
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);
    a[ALL] = 0;
    inverse.apply(DdaggerDb, a);

    diffre = 0;
    onsites(EVEN) { diffre += squarenorm(a[X] - b[X]); }
    assert(diffre * diffre < 1e-16 && "test (DdgD)^-1 DdgD");
}

// The the Hasenbusch operator
{
    hila::out0 << "Checking with Hasenbusch_operator\n";
    using dirac_base = Dirac_Wilson_evenodd<SU<N>>;
    dirac_base Dbase(0.1, U);
    using dirac = Hasenbusch_operator<dirac_base>;
    dirac D(Dbase, 0.1);
    Field<Wilson_vector<N, double>> a, b, Db, Ddaggera, DdaggerDb;
    Field<Wilson_vector<N, double>> sol;

    a[ODD] = 0;
    b[ODD] = 0;
    onsites(EVEN) {
        a[X].gaussian_random();
        b[X].gaussian_random();
        sol[X] = 0;
    }

    double diffre = 0, diffim = 0;
    D.apply(b, Db);
    D.dagger(a, Ddaggera);
    onsites(EVEN) {
        diffre += a[X].dot(Db[X]).re - Ddaggera[X].dot(b[X]).re;
        diffim += a[X].dot(Db[X]).im - Ddaggera[X].dot(b[X]).im;
    }

    assert(diffre * diffre < 1e-16 && "test Hasenbusch_operator");
    assert(diffim * diffim < 1e-16 && "test Hasenbusch_operator");

    // Now run CG on Ddaggera. b=1/D a -> Db = a
    CG<dirac> inverse(D);
    b[ALL] = 0;
    inverse.apply(Ddaggera, b);
    D.apply(b, Db);

    onsites(EVEN) { diffre += squarenorm(a[X] - Db[X]); }
    assert(diffre * diffre < 1e-16 && "test CG");

    // Run CG on a, multiply by DdaggerD and check the same
    onsites(EVEN) {
        a[X].gaussian_random();
    }
    b[ALL] = 0;
    inverse.apply(a, b);
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);

    diffre = 0;
    onsites(EVEN) { diffre += squarenorm(a[X] - DdaggerDb[X]); }
    assert(diffre * diffre < 1e-16 && "test DdgD (DdgD)^-1");

    // The other way around, multiply first
    onsites(EVEN) {
        b[X].gaussian_random();
    }
    D.apply(b, Db);
    D.dagger(Db, DdaggerDb);
    a[ALL] = 0;
    inverse.apply(DdaggerDb, a);

    diffre = 0;
    onsites(EVEN) { diffre += squarenorm(a[X] - b[X]); }
    assert(diffre * diffre < 1e-16 && "test (DdgD)^-1 DdgD");
}

hila::finishrun();
}
