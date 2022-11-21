#include "hila.h"
#include "gauge/staples.h"
#include "gauge/polyakov.h"

using mygroup = SU<3, double>;

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    double beta;
    double dt;
    int trajlen;
    int n_traj;
    int n_save;
    std::string measure_file;
    std::string config_file;
    double time_offset;
    bool poly_range_on;
    Vector<2, double> poly_range;
};


template <typename group>
void update_E(const GaugeField<group> &U, VectorField<Algebra<group>> &E, const parameters &p,
              double delta) {

    Field<group> staple;
    hila::number_type<group> eps = delta * p.beta / group::size();

    foralldir(d) {
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

template <typename T>
void get_ch_inv(const T& U, T& iU) {
    T tB[2];
    Complex< hila::number_type<T> > tc;
    int k,ip,iip;
    ip=0;
    iip=1;
    tB[ip]=1.;
    tc=trace(U);
    for(k=2; k<=T::size(); ++k) {
        tB[iip]=U*tB[ip];
        tB[iip]-=tc;
        tc=trace(U*tB[iip])/k;
        ip=iip;
        iip=(iip+1)%2;
    }
    Ui=tB[ip]/tc;
}

template <typename group>
double measure_plaq(const GaugeField<group> &U) {

    Reduction<double> plaq;
    plaq.allreduce(false).delayed(true);

    foralldir(dir1) foralldir(dir2) if (dir1 < dir2) {
        onsites(ALL) {
            plaq += 1.0 - real(trace(U[dir1][X] * U[dir2][X + dir1] * U[dir1][X + dir2].dagger() *
                                     U[dir2][X].dagger())) /
                              group::size();
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
    auto plaq = measure_plaq(U);
    auto e2 = measure_e2(E);

    return p.beta * plaq + e2 / 2;
}

template <typename group>
void measure_stuff(std::ofstream &out, const GaugeField<group> &U,
                   const VectorField<Algebra<group>> &E, int trajectory) {


    auto plaq = measure_plaq(U) / (lattice.volume() * NDIM * (NDIM - 1) / 2);

    auto e2 = measure_e2(E) / (lattice.volume() * NDIM);

    auto poly = measure_polyakov(U, e_t);

    out << trajectory << ' ' << plaq << ' ' << e2 << ' ' << poly << std::endl;
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
        update_E(U, E, p, p.dt);
        update_U(U, E, p.dt);
    }
    // and bring U and E to the same time value
    update_E(U, E, p, p.dt);
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
    p.dt = par.get("dt");
    p.trajlen = par.get("trajectory length");
    p.n_traj = par.get("number of trajectories");
    long seed = par.get("random seed");
    p.n_save = par.get("trajs/saved");
    p.measure_file = par.get("meas file");
    p.config_file = par.get("config name");
    if (par.get_item("polyakov range", {"off", "%f"}) == 1) {
        p.poly_range = par.get();
        p.poly_range_on = true;
    } else {
        p.poly_range = 0;
        p.poly_range_on = false;
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

    int start_traj = 0;
    if (!restore_checkpoint(U, start_traj, p)) {
        U = 1;
    }

    // open the measure file
    std::ofstream measfile;
    std::stringstream b;
    b << std::setprecision(10) << p.beta;
    p.measure_file += b.str();
    b.clear();

    // use append mode
    measfile.open(p.measure_file, std::ios::out | std::ios::app);
    measfile.precision(10); // use 10 digits precision by default

    double p_now = measure_polyakov(U, e_t).real();

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

            if (p_now > p.poly_range[0] && p_now < p.poly_range[1]) {
                poly_ok = true; // normal, nice branch
            } else if ((p_old < p.poly_range[0] && p_now > p_old && p_now < p.poly_range[1]) ||
                       (p_old > p.poly_range[1] && p_now < p_old && p_now > p.poly_range[0])) {
                poly_ok = true; // this is when we "search" for the range
            } else {
                poly_ok = false;
            }
        }

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

        measure_stuff(measfile, U, E, trajectory);

        if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
            checkpoint(U, trajectory, p);
            measfile.flush();
        }
    }
    measfile.close();

    hila::finishrun();
}
