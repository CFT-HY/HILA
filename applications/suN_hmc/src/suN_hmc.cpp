#include "hila.h"
#include "gauge/staples.h"
#include "gauge/polyakov.h"

using mygroup=SU<NCOLOR,double>;

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    double beta;
    double deltab;
    double dt;
    int trajlen;
    int n_traj;
    int n_save;
    std::string config_file;
    double time_offset;
};

///////////////////////////////////////////////////////////////////////////////////
// general functions

template <typename group>
double measure_plaq(const GaugeField<group>& U) {
    // measure the Wilson plaquette action
    Reduction<double> plaq;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            plaq+=1.0-real(trace(U[dir1][X]*U[dir2][X+dir1]*
                U[dir1][X+dir2].dagger()*U[dir2][X].dagger()))/group::size();
        }
    }
    return plaq.value();
}

template <typename group>
double measure_e2(const VectorField<Algebra<group>>& E) {
    // compute gauge kinetic energy from momentum field E
    Reduction<double> e2=0;
    e2.allreduce(false).delayed(true);
    foralldir(d) {
        onsites(ALL) e2+=E[d][X].squarenorm();
    }
    return e2.value()/2;
}

template <typename group>
void update_U(GaugeField<group>& U,const VectorField<Algebra<group>>& E,double delta) {
    // evolve U with momentum E over time stepe delta 
    foralldir(d) {
        onsites(ALL) U[d][X]=chexp(E[d][X]*delta)*U[d][X];
    }
}

template <typename group>
void regroup_gauge(GaugeField<group>& U) {
    // re-unitarize gauge filed
    foralldir(d) {
        onsites(ALL) U[d][X].reunitarize();
    }
}

// end general functions
///////////////////////////////////////////////////////////////////////////////////
// non-bulk-prevention functions

template <typename group>
void update_E(const GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p,double delta) {
    // compute the force for the plaquette action and use it to evolve E
    Field<group> staple;
    auto eps=delta*p.beta/group::size();
    foralldir(d) {
        staplesum(U,staple,d);
        onsites(ALL) {
            E[d][X]-=eps*(U[d][X]*staple[X].dagger()).project_to_algebra();
        }
    }
}

template <typename group>
double measure_action(const GaugeField<group>& U,const VectorField<Algebra<group>>& E,const parameters& p) {
    // measure the total action, consisting of plaquette and momentum term
    auto plaq=measure_plaq(U);
    auto e2=measure_e2(E);
    return p.beta*plaq+e2/2;
}

template <typename group>
void do_trajectory(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    // leap frog integration for normal action
    // start trajectory: advance U by half a time step
    update_U(U,E,p.dt/2);
    // main trajectory integration:
    for(int n=0; n<p.trajlen-1; n++) {
        update_E(U,E,p,p.dt);
        update_U(U,E,p.dt);
    }
    // end trajectory: bring U and E to the same time value
    update_E(U,E,p,p.dt);
    update_U(U,E,p.dt/2);
    regroup_gauge(U);
}

// end non-bulk-prevention functions
///////////////////////////////////////////////////////////////////////////////////
// bulk-prevention functions, cf. arXiv:2306.14319 (with n=2)

template <typename T>
T get_ch_inv(const T& U) {
    // compute inverse of the square matrix U, using the 
    // Cayley-Hamilton theorem (Faddeev–LeVerrier algorithm)
    T tB[2];
    int ip=0,iip=1;
    tB[ip]=1.;
    auto tc=trace(U);
    tB[iip]=U;
    for(int k=2; k<=T::size(); ++k) {
        tB[iip]-=tc;
        tB[ip]=U*tB[iip];
        tc=trace(tB[ip])/k;
        ip=iip;
        iip=(iip+1)%2;
    }
    return tB[ip]/tc;
}

template <typename T>
T get_bp_UAmat(const T& U) {
    // compute U*A(U)-matrix from Eq. (B3) of arXiv:2306.14319 for n=2
    T tA1=0.5*(1.+U);
    T tA2=get_ch_inv(tA1); 
    tA1=tA2*tA2.dagger();
    return U*tA1*tA1*tA2;
}

template <typename T>
T get_bp_iOsqmat(const T& U) {
    // compute matrix inside the trace on r.h.s. of Eq. (B1) of arXiv:2306.14319 for n=2
    T tA1=0.5*(1.+U);
    T tA2=tA1.dagger()*tA1;
    tA1=get_ch_inv(tA2);
    return tA1*tA1-1.;
}

template <typename group>
void update_E_bp(const GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p,double delta) {
    // compute force for BP action for n=2 according to Eq. (B5) of arXiv:2306.14319 and use it to evolve E
    Field<group> fmatp;
    Field<group> fmatm;
    auto eps=delta*2.0*p.beta/group::size();
    foralldir(d1) {
        foralldir(d2) if(d2!=d1) {
            onsites(ALL) {
                fmatm[X]=U[d1][X]*U[d2][X+d1]*(U[d2][X]*U[d1][X+d2]).dagger();
                fmatp[X]=get_bp_UAmat(fmatm[X]);
                fmatm[X]=(fmatp[X]*U[d2][X]).dagger()*U[d2][X];
            }
            onsites(ALL) {
                E[d1][X]-=eps*(fmatp[X]+fmatm[X-d2]).project_to_algebra();
            }
        }
    }
}

template <typename group>
double measure_plaq_bp(const GaugeField<group>& U) {
    // measure the BP plaquette action
    Reduction<double> plaq;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            plaq+=real(trace(get_bp_iOsqmat(U[dir1][X]*U[dir2][X+dir1]*
                (U[dir2][X]*U[dir1][X+dir2]).dagger())))/group::size();
        }
    }
    return plaq.value();
}

template <typename group>
double measure_action_bp(const GaugeField<group>& U,const VectorField<Algebra<group>>& E,const parameters& p) {
    // measure the total BP action, consisting of BP-plaquette and momentum term
    auto plaq=measure_plaq_bp(U);
    auto e2=measure_e2(E);
    return p.beta*plaq+e2/2;
}

template <typename group>
void do_trajectory_bp(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    // leap frog integration for BP action
    // start trajectory: advance U by half a time step
    update_U(U,E,p.dt/2);
    // main trajectory integration:
    for(int n=0; n<p.trajlen-1; n++) {
        update_E_bp(U,E,p,p.dt);
        update_U(U,E,p.dt);
    }
    // end trajectory: bring U and E to the same time value
    update_E_bp(U,E,p,p.dt);
    update_U(U,E,p.dt/2);
    regroup_gauge(U);
}

// end bulk-prevention functions
///////////////////////////////////////////////////////////////////////////////////
// measurement functions

template <typename group>
void measure_stuff(const GaugeField<group>& U,const VectorField<Algebra<group>>& E) {
    static bool first=true;
    if(first) {
        hila::out0<<"Legend MEAS: BP-plaq  plaq  E^2  P.real  P.imag\n";
        first=false;
    }
    auto plaqbp=measure_plaq_bp(U)/(lattice.volume()*NDIM*(NDIM-1)/2);
    auto plaq=measure_plaq(U)/(lattice.volume()*NDIM*(NDIM-1)/2);
    auto e2=measure_e2(E)/(lattice.volume()*NDIM);
    auto poly=measure_polyakov(U,e_t);
    hila::out0<<"MEAS "<<std::setprecision(8)<<plaqbp<<' '<<plaq<<' '<<e2<<' '<<poly<<'\n';
}

// end measurement functions
///////////////////////////////////////////////////////////////////////////////////
// load/save config functions

template <typename group>
void checkpoint(const GaugeField<group>& U,int trajectory,const parameters& p) {
    double t=hila::gettime();
    // save config
    U.config_write(p.config_file);
    // write run_status file
    if(hila::myrank()==0) {
        std::ofstream outf;
        outf.open("run_status",std::ios::out|std::ios::trunc);
        outf<<"trajectory  "<<trajectory+1<<'\n';
        outf<<"seed        "<<static_cast<uint64_t>(hila::random()*(1UL<<61))<<'\n';
        outf<<"time        "<<hila::gettime()<<'\n';
        outf.close();
    }
    std::stringstream msg;
    msg<<"Checkpointing, time "<<hila::gettime()-t;
    hila::timestamp(msg.str().c_str());
}

template <typename group>
bool restore_checkpoint(GaugeField<group>& U,int& trajectory,parameters& p) {
    uint64_t seed;
    bool ok=true;
    p.time_offset=0;
    hila::input status;
    if(status.open("run_status",false,false)) {
        hila::out0<<"RESTORING FROM CHECKPOINT:\n";
        trajectory=status.get("trajectory");
        seed=status.get("seed");
        p.time_offset=status.get("time");
        status.close();
        hila::seed_random(seed);
        U.config_read(p.config_file);
        ok=true;
    } else {
        std::ifstream in;
        in.open(p.config_file,std::ios::in|std::ios::binary);
        if(in.is_open()) {
            in.close();
            hila::out0<<"READING initial config\n";
            U.config_read(p.config_file);
            ok=true;
        } else {
            ok=false;
        }
    }
    return ok;
}

// end load/save config functions
///////////////////////////////////////////////////////////////////////////////////



int main(int argc,char** argv) {

    // hila::initialize should be called as early as possible
    hila::initialize(argc,argv);

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

    hila::out0<<"SU("<<mygroup::size()<<") HMC with bulk-prevention\n";

    hila::input par("parameters");

    CoordinateVector lsize;
    // reads NDIM numbers
    lsize=par.get("lattice size");
    // inverse gauge coupling
    p.beta=par.get("beta");
    // HMC step size
    p.dt=par.get("dt");
    // trajectory length in steps
    p.trajlen=par.get("trajectory length");
    // number of trajectories
    p.n_traj=par.get("number of trajectories");
    // random seed = 0 -> get seed from time
    long seed=par.get("random seed");
    // save config and checkpoint
    p.n_save=par.get("trajs/saved");
    // measure surface properties and print "profile"
    p.config_file=par.get("config name");

    /*
    if(hila::myrank()==0) {
        // test matrix logarithm
        mygroup tU;
        tU.random();
        std::cout<<"tU:"<<std::endl;
        std::cout<<tU<<std::endl;
        Algebra<mygroup> ltU=log(tU);
        std::cout<<"ltU:"<<std::endl;
        std::cout<<ltU<<std::endl;
        ltU*=0.5;
        tU=chexp(ltU);
        std::cout<<"sqrt(tU)=exp(ltU/2):"<<std::endl;
        std::cout<<tU<<std::endl;
        tU=tU*tU;
        std::cout<<"sqrt(tU)^2:"<<std::endl;
        std::cout<<tU<<std::endl;
    }
    */

    par.close(); // file is closed also when par goes out of scope

    // setting up the lattice is convenient to do after reading
    // the parameter
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);

    // Alloc gauge field and momenta (E)
    GaugeField<mygroup> U;
    VectorField<Algebra<mygroup>> E;

    int start_traj=0;
    if(!restore_checkpoint(U,start_traj,p)) {
        U=1;
    }

    bool searching=true;

    GaugeField<mygroup> U_old;
    for(int trajectory=start_traj; trajectory<p.n_traj; trajectory++) {

        U_old=U;
        double ttime=hila::gettime();

        foralldir(d) onsites(ALL) E[d][X].gaussian_random(sqrt(2.0));

        double act_old=measure_action_bp(U,E,p);

        do_trajectory_bp(U,E,p);

        double act_new=measure_action_bp(U,E,p);


        bool reject=hila::broadcast(exp(act_old-act_new)<hila::random());

        hila::out0<<std::setprecision(12)<<"HMC "<<trajectory
            <<(reject?" REJECT":" ACCEPT")<<" start "<<act_old<<" ds "
            <<std::setprecision(6)<<act_new-act_old;

        hila::out0<<" time "<<std::setprecision(3)<<hila::gettime()-ttime<<'\n';

        if(reject) {
            U=U_old;
        }

        hila::out0<<"Measure_start "<<trajectory<<'\n';

        measure_stuff(U,E);

        hila::out0<<"Measure_end "<<trajectory<<'\n';

        if(p.n_save>0&&(trajectory+1)%p.n_save==0) {
            checkpoint(U,trajectory,p);
        }
    }

    hila::finishrun();
}
