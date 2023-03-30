#include "hila.h"
#include "gauge/staples.h"
#include "gauge/polyakov.h"
#include "gauge/stout_smear.h"

#include <fftw3.h>

using mygroup = SU<3, double>;

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    double beta;
    double deltab;
    double dt;
    int trajlen;
    int n_traj;
    int n_save;
    int n_profile;
    std::string config_file;
    double time_offset;
    bool poly_range_on;
    double poly_probability_coeff;
    Vector<2, double> poly_range;
    int n_smear;
    double smear_coeff;
    int z_smear;
    int n_surface;
    int n_dump_polyakov;
};

template <typename T>
T get_ch_inv(const T &U) {
    T tB[2];
    Complex<hila::scalar_type<T>> tc;
    int ip, iip;
    ip = 0;
    iip = 1;
    tB[ip] = 1.;
    tc = trace(U);
    for (int k = 2; k <= T::size(); ++k) {
        tB[iip] = U * tB[ip];
        tB[iip] -= tc;
        tc = trace(U * tB[iip]) / k;
        ip = iip;
        iip = (iip + 1) % 2;
    }
    return tB[ip] / tc;
}

template <typename T>
T get_bp_Amat(const T &U) {
    T tA1;
    T tA2;
    tA1 = 0.5 * U;
    tA1 += 0.5;
    tA2 = get_ch_inv(tA1);
    tA1 = tA2 * tA2.dagger();
    return tA1 * tA1 * tA2;
}

template <typename T>
T get_bp_iOsqmat(const T &U) {
    T tA1;
    T tA2;
    tA1 = 0.5 * U;
    tA1 += 0.5;
    tA2 = tA1.dagger() * tA1;
    tA1 = get_ch_inv(tA2);
    tA2 = tA1 * tA1;
    tA2 -= 1.;
    return tA2;
}

template <typename T>
void plaqpm(const GaugeField<T> &U, Field<T> &plaqp, Direction d1, Direction d2) {

    Field<T> lower;

    if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(ALL) {
            lower[X] = U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        onsites(ALL) {
            auto p1 = U[d1][X] * U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            auto p2 = U[d1][X] * lower[X - d2].dagger();
            plaqp[X] = p1 * get_bp_Amat(p1) + p2 * get_bp_Amat(p2);
        }
    }
}

template <typename T>
void plaqpm_db(const GaugeField<T> &U, Field<T> &plaqp, Direction d1, Direction d2, double deltab) {

    Field<T> lower;

    if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(ALL) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            auto p2 = U[d1][X + d2] * (U[d2][X].dagger() * U[d1][X] * U[d2][X + d1]).dagger();
            lower[X] = m * p2 * get_bp_Amat(p2);
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        onsites(ALL) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            auto p1 = U[d1][X] * U[d2][X + d1] * (U[d2][X] * U[d1][X + d2]).dagger();
            plaqp[X] = m * p1 * get_bp_Amat(p1) + lower[X - d2];
        }
    }
}


template <typename group>
void update_E_bp(const GaugeField<group> &U, VectorField<Algebra<group>> &E, const parameters &p,
                 double delta) {

    Field<group> plaqp;
    hila::scalar_type<group> eps = delta * 2.0 * p.beta / group::size();

    foralldir(d1) {
        foralldir(d2) if (d2 != d1) {
            if (p.deltab == 0) {
                plaqpm(U, plaqp, d1, d2);
            } else {
                plaqpm_db(U, plaqp, d1, d2, p.deltab);
            }

            onsites(ALL) {
                E[d1][X] -= eps * plaqp[X].project_to_algebra();
            }
        }
    }
}

template <typename T>
void staplesum_db(const GaugeField<T> &U, Field<T> &staples, Direction d1, double deltab) {

    Field<T> lower;

    // zero staples just inc case
    staples = 0;
    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, ALL);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(ALL) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            lower[X] = m * U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        onsites(ALL) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            staples[X] += m * U[d2][X] * U[d1][X + d2] * U[d2][X + d1].dagger() + lower[X - d2];
        }
    }
}

template <typename group>
void update_E(const GaugeField<group> &U, VectorField<Algebra<group>> &E, const parameters &p,
              double delta) {

    Field<group> staple;
    hila::scalar_type<group> eps = delta * p.beta / group::size();

    foralldir(d) {
        if (p.deltab != 0)
            staplesum_db(U, staple, d, p.deltab);
        else
            staplesum(U, staple, d);

        onsites(ALL) {
            E[d][X] -= eps * (U[d][X] * staple[X].dagger()).project_to_algebra();
        }
    }
}


template <typename group>
void update_U(GaugeField<group> &U, const VectorField<Algebra<group>> &E, double delta) {

    foralldir(d) {
        onsites(ALL) U[d][X] = exp(E[d][X] * delta) * U[d][X];
    }
}

template <typename group>
void regroup_gauge(GaugeField<group> &U) {
    foralldir(d) {
        onsites(ALL) U[d][X].reunitarize();
    }
}

template <typename group>
double measure_plaq_bp(const GaugeField<group> &U, double db = 0.0) {

    Reduction<double> plaq;
    plaq.allreduce(false).delayed(true);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        if (db == 0.0) {

            onsites(ALL) {
                plaq += real(trace(get_bp_iOsqmat(U[dir1][X] * U[dir2][X + dir1] *
                                                  (U[dir2][X] * U[dir1][X + dir2]).dagger()))) /
                        group::size();
            }
        } else {

            onsites(ALL) {

                double c;
                if (2 * X.z() < lattice.size(e_z))
                    c = 1 - db;
                else
                    c = 1 + db;

                plaq += c *
                        real(trace(get_bp_iOsqmat(U[dir1][X] * U[dir2][X + dir1] *
                                                  (U[dir2][X] * U[dir1][X + dir2]).dagger()))) /
                        group::size();
            }
        }
    }
    return plaq.value();
}


template <typename group>
double measure_plaq(const GaugeField<group> &U, double db = 0.0) {

    Reduction<double> plaq;
    plaq.allreduce(false);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        if (db == 0.0) {
            onsites(ALL) {
                plaq += 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                         U[dir1][X + dir2].dagger() * U[dir2][X].dagger())) /
                                  group::size();
            }
        } else {
            onsites(ALL) {
                double c;
                if (2 * X.z() < lattice.size(e_z))
                    c = 1 - db;
                else
                    c = 1 + db;

                plaq += c * (1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] *
                                              U[dir1][X + dir2].dagger() * U[dir2][X].dagger())) /
                                       group::size());
            }
        }
    }
    return plaq.value();
}

template <typename group>
double measure_e2(const VectorField<Algebra<group>> &E) {

    Reduction<double> e2 = 0;
    e2.allreduce(false).delayed(true);

    foralldir(d) {
        onsites(ALL) e2 += E[d][X].squarenorm();
    }
    return e2.value() / 2;
}

template <typename group>
double measure_action(const GaugeField<group> &U, const VectorField<Algebra<group>> &E,
                      const parameters &p) {

    auto plaq = measure_plaq_bp(U, p.deltab);
    auto e2 = measure_e2(E);

    return p.beta * plaq + e2 / 2;
}

/////////////////////////////////////////////////////////////////////////////
/// Measure polyakov "field"
///

template <typename T>
void measure_polyakov_field(const GaugeField<T> &U, Field<float> &pl) {
    Field<T> polyakov = U[e_t];

    // mult links so that polyakov[X.dir == 0] contains the polyakov loop
    for (int plane = lattice.size(e_t) - 2; plane >= 0; plane--) {

        // safe_access(polyakov) pragma allows the expression below, otherwise
        // hilapp would reject it because X and X+dir can refer to the same
        // site on different "iterations" of the loop.  However, here this
        // is restricted on single dir-plane so it works but we must tell it to hilapp.

#pragma hila safe_access(polyakov)
        onsites(ALL) {
            if (X.coordinate(e_t) == plane) {
                polyakov[X] = U[e_t][X] * polyakov[X + e_t];
            }
        }
    }

    onsites(ALL) if (X.coordinate(e_t) == 0) {
        pl[X] = real(trace(polyakov[X]));
    }
}

////////////////////////////////////////////////////////////////////

void smear_polyakov_field(Field<float> &pl, int nsmear, float smear_coeff) {

    if (nsmear > 0) {
        Field<float> pl2 = 0;

        for (int i = 0; i < nsmear; i++) {
            onsites(ALL) if (X.coordinate(e_t) == 0) {
                pl2[X] = pl[X] + smear_coeff * (pl[X + e_x] + pl[X - e_x] + pl[X + e_y] +
                                                pl[X - e_y] + pl[X + e_z] + pl[X - e_z]);
            }

            onsites(ALL) if (X.coordinate(e_t) == 0) {
                pl[X] = pl2[X] / (1 + 6 * smear_coeff);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

std::vector<float> measure_polyakov_profile(Field<float> &pl, std::vector<float> &pro1) {
    ReductionVector<float> p(lattice.size(e_z)), p1(lattice.size(e_z));
    p.allreduce(false);
    p1.allreduce(false);
    onsites(ALL) if (X.coordinate(e_t) == 0) {
        p[X.z()] += pl[X];
        if (X.x() == 0 && X.y() == 0)
            p1[X.z()] += pl[X];
    }
    pro1 = p1.vector();
    return p.vector();
}


///////////////////////////////////////////////////////////////////////////////////

template <typename group>
void measure_stuff(const GaugeField<group> &U, const VectorField<Algebra<group>> &E) {

    static bool first = true;

    if (first) {
        hila::out0 << "Legend MEAS: plaq  E^2  P.real  P.imag\n";
        first = false;
    }

    auto plaq = measure_plaq(U) / (lattice.volume() * NDIM * (NDIM - 1) / 2);

    auto plaqbp = measure_plaq_bp(U) / (lattice.volume() * NDIM * (NDIM - 1) / 2);

    auto e2 = measure_e2(E) / (lattice.volume() * NDIM);


    std::vector<Complex<double>> profile;

    auto poly = measure_polyakov(U, e_t);

    hila::out0 << "MEAS " << std::setprecision(8) << plaqbp << ' ' << plaq << ' ' << e2 << ' '
               << poly << '\n';
}

///////////////////////////////////////////////////////////////////////////////////

// helper function to get valid z-coordinate index

int z_ind(int z) {
    return (z + lattice.size(e_z)) % lattice.size(e_z);
}

///////////////////////////////////////////////////////////////////////////////////

template <typename group>
void measure_polyakov_surface(const GaugeField<group> &U, const parameters &p, int traj) {

    Field<float> pl;
    measure_polyakov_field(U, pl);

    std::vector<float> profile, profile_unsmear, prof1;
    profile_unsmear = measure_polyakov_profile(pl, prof1);

    smear_polyakov_field(pl, p.n_smear, p.smear_coeff);

    if (p.z_smear > 0) {
        Field<float> pl2 = 0;
        for (int i = 0; i < p.z_smear; i++) {
            onsites(ALL) if (X.coordinate(e_t) == 0) {
                pl2[X] = pl[X] + p.smear_coeff * (pl[X + e_z] + pl[X - e_z]);
            }
            onsites(ALL) if (X.coordinate(e_t) == 0) {
                pl[X] = pl2[X] / (1 + 2 * p.smear_coeff);
            }
        }
    }

    profile = measure_polyakov_profile(pl, prof1);

    hila::out0 << std::setprecision(5);
    double m = 1.0 / (lattice.size(e_x) * lattice.size(e_y));
    for (int i = 0; i < profile.size(); i++) {
        profile[i] *= m;
        hila::out0 << "PRO " << i << ' ' << profile[i] << ' ' << prof1[i] << ' '
                   << profile_unsmear[i] * m << '\n';
    }

    float min = 1e8, max = 0;
    int minloc, maxloc;
    for (int i = 0; i < profile.size(); i++) {
        if (min > profile[i]) {
            min = profile[i];
            minloc = i;
        }
        if (max < profile[i]) {
            max = profile[i];
            maxloc = i;
        }
    }

    // find the surface between minloc and maxloc
    float surface_level = max * 0.5; // assume min is really 0
    int area = lattice.size(e_x) * lattice.size(e_y);

    hila::out0 << "Surface level " << surface_level << '\n';

    int startloc;
    if (maxloc > minloc)
        startloc = (maxloc + minloc) / 2;
    else
        startloc = ((maxloc + minloc + lattice.size(e_z)) / 2) % lattice.size(e_z);

    std::vector<float> surf;
    if (hila::myrank() == 0)
        surf.resize(area);

    hila::out0 << std::setprecision(6);

    for (int y = 0; y < lattice.size(e_y); y++)
        for (int x = 0; x < lattice.size(e_x); x++) {
            auto line = pl.get_slice({x, y, -1, 0});
            if (hila::myrank() == 0) {
                // start search of the surface from the center between min and max
                int z = startloc;
                if (line[z] > surface_level) {
                    while (line[z_ind(z)] > surface_level && startloc - z < lattice.size(e_z) * 0.8)
                        z--;
                } else {
                    while (line[z_ind(z + 1)] <= surface_level &&
                           z - startloc < lattice.size(e_z) * 0.8)
                        z++;
                }


                // do linear interpolation
                surf[x + y * lattice.size(e_x)] =
                    z + (surface_level - line[z_ind(z)]) / (line[z_ind(z + 1)] - line[z_ind(z)]);

                if (p.n_surface > 0 && (traj + 1) % p.n_surface == 0) {
                    hila::out0 << "SURF " << x << ' ' << y << ' ' << surf[x + y * lattice.size(e_x)]
                               << '\n';
                }
            }
        }

    if (hila::myrank() == 0) {
        // do fft for the surface
        static bool first = true;
        static Complex<double> *buf;
        static fftw_plan fftwplan;

        if (first) {
            first = false;

            buf = (Complex<double> *)fftw_malloc(sizeof(Complex<double>) * area);

            // note: we had x as the "fast" dimension, but fftw wants the 2nd dim to be
            // the "fast" one. thus, first y, then x.
            fftwplan = fftw_plan_dft_2d(lattice.size(e_y), lattice.size(e_x), (fftw_complex *)buf,
                                        (fftw_complex *)buf, FFTW_FORWARD, FFTW_ESTIMATE);
        }

        for (int i = 0; i < area; i++) {
            buf[i] = surf[i];
        }

        fftw_execute(fftwplan);

        constexpr int pow_size = 200;

        std::vector<double> npow(pow_size);
        std::vector<int> hits(pow_size);

        for (int i = 0; i < area; i++) {
            int x = i % lattice.size(e_x);
            int y = i / lattice.size(e_x);
            x = (x <= lattice.size(e_x) / 2) ? x : (lattice.size(e_x) - x);
            y = (y <= lattice.size(e_y) / 2) ? y : (lattice.size(e_y) - y);

            int k = x * x + y * y;
            if (k < pow_size) {
                npow[k] += buf[i].squarenorm() / (area * area);
                hits[k]++;
            }
        }

        for (int i = 0; i < pow_size; i++) {
            if (hits[i] > 0)
                hila::out0 << "POW " << i << ' ' << npow[i] / hits[i] << ' ' << hits[i] << '\n';
        }
    }
}


////////////////////////////////////////////////////////////////

template <typename group>
void checkpoint(const GaugeField<group> &U, int trajectory, const parameters &p) {

    double t = hila::gettime();
    // save config
    U.config_write(p.config_file);

    // write run_status file
    if (hila::myrank() == 0) {
        std::ofstream outf;
        outf.open("run_status", std::ios::out | std::ios::trunc);
        outf << "trajectory  " << trajectory + 1 << '\n';
        outf << "seed        " << static_cast<uint64_t>(hila::random() * (1UL << 61)) << '\n';
        outf << "time        " << hila::gettime() << '\n';
        outf.close();
    }

    std::stringstream msg;
    msg << "Checkpointing, time " << hila::gettime() - t;
    hila::timestamp(msg.str().c_str());
}

////////////////////////////////////////////////////////////////

template <typename group>
bool restore_checkpoint(GaugeField<group> &U, int &trajectory, parameters &p) {

    uint64_t seed;
    bool ok = true;
    p.time_offset = 0;

    hila::input status;
    if (status.open("run_status", false, false)) {

        hila::out0 << "RESTORING FROM CHECKPOINT:\n";

        trajectory = status.get("trajectory");
        seed = status.get("seed");
        p.time_offset = status.get("time");
        status.close();

        hila::seed_random(seed);

        U.config_read(p.config_file);

        ok = true;
    } else {

        std::ifstream in;
        in.open(p.config_file, std::ios::in | std::ios::binary);

        if (in.is_open()) {
            in.close();

            hila::out0 << "READING initial config\n";

            U.config_read(p.config_file);

            ok = true;
        } else {

            ok = false;
        }
    }

    return ok;
}

////////////////////////////////////////////////////////////////

template <typename group>
void do_trajectory(GaugeField<group> &U, VectorField<Algebra<group>> &E, const parameters &p) {

    update_U(U, E, p.dt / 2);
    for (int n = 0; n < p.trajlen - 1; n++) {
        update_E_bp(U, E, p, p.dt);
        update_U(U, E, p.dt);
    }
    // and bring U and E to the same time value
    update_E_bp(U, E, p, p.dt);
    update_U(U, E, p.dt / 2);
    regroup_gauge(U);
}

/////////////////////////////////////////////////////////////////


int main(int argc, char **argv) {

    // hila::initialize should be called as early as possible
    hila::initialize(argc, argv);

    // hila provides an input class hila::input, which is
    // a convenient way to read in parameters from input files.
    // parameters are presented as key - value pairs, as an example
    //  " lattice size  64, 64, 64, 64"
    // is read below.
    //
    // Values are broadcast to all MPI nodes.
    //
    // .get() -method can read many different input types,
    // see file "input.h" for documentation

    parameters p;

    hila::input par("parameters");

    CoordinateVector lsize;
    lsize = par.get("lattice size"); // reads NDIM numbers

    p.beta = par.get("beta");
    // deltab sets system to different beta on different sides, by beta*(1 +- deltab)
    // use for initial config generation only
    p.deltab = par.get("delta beta fraction");
    p.dt = par.get("dt");
    // trajectory length in steps
    p.trajlen = par.get("trajectory length");
    p.n_traj = par.get("number of trajectories");
    // random seed = 0 -> get seed from time
    long seed = par.get("random seed");
    // save config and checkpoint
    p.n_save = par.get("trajs/saved");
    // measure surface properties and print "profile"
    p.config_file = par.get("config name");

    // if polyakov range is off, do nothing with
    if (par.get_item("polyakov range", {"off", "%f"}) == 1) {
        p.poly_range = par.get();
        p.poly_range_on = true;
        p.poly_probability_coeff = par.get("polyakov probability coeff");

    } else {
        p.poly_range = 0;
        p.poly_range_on = false;
    }

    p.n_profile = par.get("trajs/profile meas");

    if (p.n_profile) {
        p.n_smear = par.get("smearing steps");
        p.smear_coeff = par.get("smear coefficient");
        p.z_smear = par.get("z smearing steps");
        p.n_surface = par.get("trajs/surface");
        p.n_dump_polyakov = par.get("trajs/polyakov dump");
    }


    par.close(); // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);

    // Alloc gauge field and momenta (E)
    GaugeField<mygroup> U;
    VectorField<Algebra<mygroup>> E;

    U = 1;

    int start_traj = 0;
    if (!restore_checkpoint(U, start_traj, p)) {
        U = 1;
        if (p.n_profile > 0) {
            foralldir(d) onsites(ALL) {
                double mag;
                // mag =
                //     0.4 +
                //     0.3 * (1 + cos(2 * M_PI * (X.z() - lattice.size(e_z) / 4) / lattice.size(e_z)));
                if (X.z() <= lattice.size(e_z) / 2) {
                    mag = 1;
                } else {
                    mag = 0.2;
                }
                mygroup u;
                u.gaussian_random(mag);
                U[d][X] += u * 0.01;
                U[d][X].reunitarize();
            }
        }
    }

    double p_now = measure_polyakov(U, e_t).real();
    bool searching = true;

    GaugeField<mygroup> U_old;
    for (int trajectory = start_traj; trajectory < p.n_traj; trajectory++) {

        U_old = U;
        double ttime = hila::gettime();

        double p_old = p_now;

        foralldir(d) onsites(ALL) E[d][X].gaussian_random(sqrt(2.0));

        double act_old = measure_action(U, E, p);

        do_trajectory(U, E, p);

        double act_new = measure_action(U, E, p);

        bool poly_ok = true;
        if (p.poly_range_on) {
            p_now = measure_polyakov(U, e_t).real();

            double db = 0;
            if (p_now > p.poly_range[1]) {
                if (p_old > p.poly_range[1])
                    db = p_now - p_old;
                else
                    db = p_now - p.poly_range[1];
            } else if (p_now < p.poly_range[0]) {
                if (p_old < p.poly_range[0])
                    db = p_old - p_now;
                else
                    db = p.poly_range[0] - p_now;
            }

            poly_ok = (hila::random() < exp(-p.poly_probability_coeff * db));
            //   poly_ok = hila::broadcast(poly_ok);

            // if (p_now > p.poly_range[0] && p_now < p.poly_range[1]) {
            //     poly_ok = true;    // normal, nice branch
            //     searching = false; // turn off search
            // } else if ((p_old < p.poly_range[0] && p_now > p_old && p_now <
            // p.poly_range[1])
            // ||
            //            (p_old > p.poly_range[1] && p_now < p_old && p_now >
            //            p.poly_range[0]))
            //            {
            //     poly_ok = true; // this is when we "search" for the range
            // } else {
            //     poly_ok = false;
            // }
        }
        // if (!poly_ok && searching) {
        //     poly_ok = (hila::random() < 0.2);
        // }
        bool reject = hila::broadcast(!poly_ok || (exp(act_old - act_new) < hila::random()));

        hila::out0 << std::setprecision(12) << "HMC " << trajectory
                   << (reject ? " REJECT" : " ACCEPT") << " start " << act_old << " ds "
                   << std::setprecision(6) << act_new - act_old;
        if (p.poly_range_on)
            hila::out0 << " ploop " << p_now;

        hila::out0 << " time " << std::setprecision(3) << hila::gettime() - ttime << '\n';

        if (reject) {
            U = U_old;
            p_now = p_old;
        }

        hila::out0 << "Measure_start " << trajectory << '\n';

        measure_stuff(U, E);

        if (p.n_profile && (trajectory + 1) % p.n_profile == 0) {
            measure_polyakov_surface(U, p, trajectory);
        }

        if (p.n_dump_polyakov && (trajectory + 1) % p.n_dump_polyakov == 0) {
            Field<float> pl;
            std::ofstream poly;
            if (hila::myrank() == 0) {
                poly.open("polyakov", std::ios::out | std::ios::app);
            }
            measure_polyakov_field(U, pl);
            pl.write_slice(poly, {-1, -1, -1, 0});
        }


        hila::out0 << "Measure_end " << trajectory << std::endl;

        if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
            checkpoint(U, trajectory, p);
        }
    }

    hila::finishrun();
}
