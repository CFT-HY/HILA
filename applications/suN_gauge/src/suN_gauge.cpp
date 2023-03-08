#include "hila.h"
#include "gauge/staples.h"
#include "gauge/polyakov.h"
#include "gauge/stout_smear.h"

#include "gauge/sun_heatbath.h"
#include "gauge/sun_overrelax.h"

#include <fftw3.h>

using mygroup = SU<NCOLOR, double>;

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    double beta;
    double deltab;
    int n_overrelax;
    int n_update;
    int n_trajectories;
    int n_thermal;
    int n_save;
    int n_profile;
    std::string config_file;
    double time_offset;
    bool poly_range_on;
    double poly_min, poly_m2;
    int n_smear;
    double smear_coeff;
    int z_smear;
    int n_surface;
    int n_dump_polyakov;
};

template <typename T>
void staplesum_db(const GaugeField<T> &U, Field<T> &staples, Direction d1, Parity par,
                  double deltab) {

    Field<T> lower;

    // zero staples just in case
    staples[par] = 0;

    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        U[d2].start_gather(d1, ALL);
        U[d1].start_gather(d2, par);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(opp_parity(par)) {
            double m;
            if (2 * X.z() < lattice.size(e_z))
                m = 1.0 - deltab;
            else
                m = 1.0 + deltab;

            lower[X] = m * U[d2][X].dagger() * U[d1][X] * U[d2][X + d1];
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        onsites(par) {
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
void reunitarize_gauge(GaugeField<group> &U) {
    foralldir(d) {
        onsites(ALL) U[d][X].reunitarize();
    }
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
double measure_action(const GaugeField<group> &U, const VectorField<Algebra<group>> &E,
                      const parameters &p) {

    auto plaq = measure_plaq_bp(U, p.deltab);

    return p.beta * plaq;
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
void measure_stuff(const GaugeField<group> &U, const parameters &p) {

    static bool first = true;

    if (first) {
        hila::out0 << "Legend MEAS: plaq  P.real  P.imag";
        if (p.poly_range_on)
            hila::out0 << " V(polyakov)";
        hila::out0 << '\n';

        first = false;
    }

    hila::out0 << "MEAS " << std::setprecision(8);

    hila::out0 << measure_plaq(U) / (lattice.volume() * NDIM * (NDIM - 1) / 2);

    auto poly = measure_polyakov(U);
    auto poly_x = measure_polyakov(U,e_x);

    hila::out0 << ' ' << poly << ' ' << poly_x;

    if (p.poly_range_on) {
        hila::out0 << ' ' << polyakov_potential(p, poly.real());
    }

    hila::out0 << std::endl;

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
void checkpoint(const GaugeField<group> &U, int iteration, const parameters &p) {

    double t = hila::gettime();
    // save config
    U.config_write(p.config_file);

    // write run_status file
    if (hila::myrank() == 0) {
        std::ofstream outf;
        outf.open("run_status", std::ios::out | std::ios::trunc);
        outf << "iteration   " << iteration + 1 << '\n';
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

        trajectory = status.get("iteration");
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
void update(GaugeField<group> &U, const parameters &p, bool relax) {

    Field<group> staples;

    static hila::timer hb_timer("Heatbath");
    static hila::timer or_timer("Overrelax");
    static hila::timer staples_timer("Staplesum");

    foralldir(d) {
        for (Parity par : {EVEN, ODD}) {

            staples_timer.start();
            if (p.deltab != 0)
                staplesum_db(U, staples, d, par, p.deltab);
            else
                staplesum(U, staples, d, par);
            staples_timer.stop();

            if (relax) {

                or_timer.start();
                onsites(par) {
                    suN_overrelax(U[d][X], staples[X]);
                }
                or_timer.stop();

            } else {

                hb_timer.start();
                onsites(par) {
                    suN_heatbath(U[d][X], staples[X], p.beta);
                }
                hb_timer.stop();
            }
        }
    }
}

////////////////////////////////////////////////////////////////

template <typename group>
void do_trajectory(GaugeField<group> &U, const parameters &p) {

    for (int n = 0; n < p.n_update; n++) {
        for (int i = 0; i < p.n_overrelax; i++) {
            update(U, p, true);
        }
        update(U, p, false);
    }
    reunitarize_gauge(U);
}

////////////////////////////////////////////////////////////////

double polyakov_potential(const parameters &p, const double poly) {
    return p.poly_m2 * (sqr(poly - p.poly_min));
}


////////////////////////////////////////////////////////////////

bool accept_polyakov(const parameters &p, const double p_old, const double p_new) {

    double dpot = polyakov_potential(p, p_old) - polyakov_potential(p, p_new);

    bool accept = hila::broadcast(hila::random() < exp(-dpot));

    return accept;
}

////////////////////////////////////////////////////////////////

template <typename group>
double trajectory_with_range(GaugeField<group> &U, const parameters &p, double &poly) {

    double acc_or = 0;
    double acc_hb = 0;

    Field<group> Ut_old;
    Ut_old = U[e_t]; // need to store the t-direction field

    double p_now;

    for (int n = 0; n < p.n_update; n++) {
        for (int i = 0; i < p.n_overrelax; i++) {

            // spatial links always updated, t-links are conditional
            // This works only if the t-links are the last to be updated inside update routines!

            update(U, p, true);

            p_now = measure_polyakov(U).real();

            bool acc_update = accept_polyakov(p, poly, p_now);

            if (acc_update) {
                poly = p_now;
                acc_or++;
                Ut_old = U[e_t]; // prepare for next step
            } else {
                U[e_t] = Ut_old; // restore rejected
            }
        }

        // and then HB -- Ut_old is valid in any case

        update(U, p, false);

        double p_hb = measure_polyakov(U).real();
        bool acc_update = accept_polyakov(p, poly, p_hb);

        if (acc_update) {
            poly = p_hb;
            acc_hb++;
            Ut_old = U[e_t]; // prepare for next OR
        } else {
            U[e_t] = Ut_old;
        }
    }

    reunitarize_gauge(U);
    return (acc_or + acc_hb) / (p.n_update * (p.n_overrelax + 1));
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

    hila::out0 << "SU(" << mygroup::size() << ") heat bath + overrelax update\n";

    hila::input par("parameters");

    CoordinateVector lsize;
    lsize = par.get("lattice size"); // reads NDIM numbers

    p.beta = par.get("beta");
    // deltab sets system to different beta on different sides, by beta*(1 +- deltab)
    // use for initial config generation only
    p.deltab = par.get("delta beta fraction");
    // trajectory length in steps
    p.n_overrelax = par.get("overrelax steps");
    p.n_update = par.get("updates in trajectory");
    p.n_trajectories = par.get("trajectories");
    p.n_thermal = par.get("thermalization");

    // random seed = 0 -> get seed from time
    long seed = par.get("random seed");
    // save config and checkpoint
    p.n_save = par.get("traj/saved");
    // measure surface properties and print "profile"
    p.config_file = par.get("config name");

    // if polyakov range is off, do nothing with
    if (par.get_item("polyakov potential", {"off", "min"}) == 1) {
        p.poly_range_on = true;
        p.poly_min = par.get();
        p.poly_m2 = par.get("mass2");
    } else {
        p.poly_range_on = false;
    }

    if (par.get_item("updates/profile meas", {"off", "%f"}) == 1) {
        p.n_profile = par.get();
    } else {
        p.n_profile = 0;
    }

    if (p.n_profile) {
        p.n_smear = par.get("smearing steps");
        p.smear_coeff = par.get("smear coefficient");
        p.z_smear = par.get("z smearing steps");
        p.n_surface = par.get("traj/surface");
        p.n_dump_polyakov = par.get("traj/polyakov dump");
    }

    par.close(); // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);

    // Alloc gauge field
    GaugeField<mygroup> U;

    U = 1;

    // use negative trajectory for thermal
    int start_traj = -p.n_thermal;

    hila::timer update_timer("Updates");
    hila::timer measure_timer("Measurements");

    restore_checkpoint(U, start_traj, p);

    double p_now = measure_polyakov(U).real();

    for (int trajectory = start_traj; trajectory < p.n_trajectories; trajectory++) {

        double ttime = hila::gettime();

        update_timer.start();

        double acc = 0;
        if (p.poly_range_on) {
            acc += trajectory_with_range(U, p, p_now);
        } else {
            do_trajectory(U, p);
        }

        update_timer.stop();

        // trajectory is negative during thermalization
        if (trajectory >= 0) {
            measure_timer.start();

            hila::out0 << "Measure_start " << trajectory << '\n';

            measure_stuff(U, p);

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

            measure_timer.stop();
        }

        if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
            checkpoint(U, trajectory, p);
        }
    }

    hila::finishrun();
}
