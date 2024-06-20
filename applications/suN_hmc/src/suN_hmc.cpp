#include "hila.h"
#include "gauge/polyakov.h"
#include "wilson_plaquette_action.h"
#include "bulk_prevention_action.h"
#include "wilson_line_and_force.h"
#include "improved_action.h"
#include "wilson_flow.h"
#include "energy_and_topo_charge_clover.h"
#include "energy_and_topo_charge_log.h"
#include "clover_action.h"
#include <string>

using ftype=double;
using mygroup=SU<NCOLOR,ftype>;

template<typename ... Args>
std::string string_format(const std::string& format,Args ... args) {
    // wrapper for std::snprintf which sets up buffer of required size

    // determine required buffer size :
    int size_s=std::snprintf(nullptr,0,format.c_str(),args ...)+1;
    if(size_s<=0) { 
        throw std::runtime_error("Error during formatting."); 
    }

    // allocate buffer :
    auto size=static_cast<size_t>(size_s);
    std::unique_ptr<char[]> buf(new char[size]);

    // write formatted string to buffer :
    std::snprintf(buf.get(),size,format.c_str(),args ...);

    return std::string(buf.get(),buf.get()+size-1); 
}

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    ftype beta;          // inverse gauge coupling
    ftype dt;            // HMC time step
    int trajlen;         // number of HMC time steps per trajectory
    int n_traj;          // number of trajectories to generate
    int n_therm;         // number of thermalization trajectories (counts only accepted traj.)
    int wflow_freq;      // number of trajectories between wflow measurements
    ftype wflow_max_l;   // flow scale at which wilson flow stops
    ftype wflow_l_step;  // flow scale interval between flow measurements
    ftype wflow_a_accu;  // desired absolute accuracy of wilson flow integration steps
    ftype wflow_r_accu;  // desired relative accuracy of wilson flow integration steps
    int n_save;          // number of trajectories between config. check point 
    std::string config_file;
    ftype time_offset;
};

// DBW2 action parameters:
//ftype c12=-1.4088;     // rectangle weight
//ftype c11=1.0-8.0*c12; // plaquette weight

// Iwasaki action parameters:
ftype c12=-0.331;     // rectangle weight
ftype c11=1.0-8.0*c12; // plaquette weight

// LW action parameters:
//ftype c12=-1.0/12.0;     // rectangle weight
//ftype c11=1.0-8.0*c12; // plaquette weight

// Wilson plaquette action parameters:
//ftype c12=0;
//ftype c11=1.0;


///////////////////////////////////////////////////////////////////////////////////
// general functions

template <typename group,typename atype=hila::arithmetic_type<group>>
void update_U(GaugeField<group>& U,const VectorField<Algebra<group>>& E,atype delta) {
    // evolve U with momentum E over time step delta 
    foralldir(d) {
        onsites(ALL) U[d][X]=chexp(E[d][X]*delta)*U[d][X];
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_e2(const VectorField<Algebra<group>>& E) {
    // compute gauge kinetic energy from momentum field E
    Reduction<atype> e2=0;
    e2.allreduce(false).delayed(true);
    foralldir(d) {
        onsites(ALL) e2+=E[d][X].squarenorm();
    }
    return e2.value()/2;
}

// end general functions
///////////////////////////////////////////////////////////////////////////////////
// Wilson plaquette action functions

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_local_wplaq(const GaugeField<group>& U) {
    return measure_s_wplaq(U);
}

template <typename group>
void do_trajectory_wplaq(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    // leap frog integration for normal action

    // start trajectory: advance U by half a time step
    update_U(U,E,p.dt/2);
    // main trajectory integration:
    for(int n=0; n<p.trajlen-1; n++) {
        update_E_wplaq(U,E,p.beta*p.dt);
        update_U(U,E,p.dt);
    }
    // end trajectory: bring U and E to the same time
    update_E_wplaq(U,E,p.beta*p.dt);
    update_U(U,E,p.dt/2);

    U.reunitarize_gauge();
}

// end Wilson plaquette action functions
///////////////////////////////////////////////////////////////////////////////////
// improved action functions

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_local_impr(const GaugeField<group>& U) {
    return measure_s_impr_f2(U,c11,c12);
}

template <typename group>
void do_trajectory_impr(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    // leap frog integration for normal action

    // start trajectory: advance U by half a time step
    update_U(U,E,p.dt/2);
    // main trajectory integration:
    for(int n=0; n<p.trajlen-1; n++) {
        update_E_impr_f2(U,E,c11,c12,p.beta*p.dt);
        update_U(U,E,p.dt);
    }
    // end trajectory: bring U and E to the same time
    update_E_impr_f2(U,E,c11,c12,p.beta*p.dt);
    update_U(U,E,p.dt/2);
    
    U.reunitarize_gauge();
}

// end improved action functions
///////////////////////////////////////////////////////////////////////////////////
// bulk-prevention functions

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_local_bp(const GaugeField<group>& U) {
    return measure_s_bp(U);
}

template <typename group>
void do_trajectory_bp(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    // leap frog integration for BP action

    // start trajectory: advance U by half a time step
    update_U(U,E,p.dt/2);
    // main trajectory integration:
    for(int n=0; n<p.trajlen-1; n++) {
        update_E_bp(U,E,p.beta*p.dt);
        update_U(U,E,p.dt);
    }
    // end trajectory: bring U and E to the same time
    update_E_bp(U,E,p.beta*p.dt);
    update_U(U,E,p.dt/2);

    U.reunitarize_gauge();
}

// end bulk-prevention functions
///////////////////////////////////////////////////////////////////////////////////
// clover action functions

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s_local_clover(const GaugeField<group>& U) {
    return measure_s_clover_action(U);
}

template <typename group>
void do_trajectory_clover(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    // leap frog integration for BP action

    // start trajectory: advance U by half a time step
    update_U(U,E,p.dt/2);
    // main trajectory integration:
    for(int n=0; n<p.trajlen-1; n++) {
        update_E_clover_action(U,E,p.beta*p.dt);
        update_U(U,E,p.dt);
    }
    // end trajectory: bring U and E to the same time
    update_E_clover_action(U,E,p.beta*p.dt);
    update_U(U,E,p.dt/2);

    U.reunitarize_gauge();
}

// end clover action functions
///////////////////////////////////////////////////////////////////////////////////
// HMC action and force wrappers

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_s(const GaugeField<group>& U) {
    //atype res=measure_s_local_wplaq(U);
    atype res=measure_s_local_bp(U);
    //atype res=measure_s_local_impr(U);
    //atype res=measure_s_local_clover(U);
    return res;
}

template <typename group>
void do_trajectory(GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p) {
    //do_trajectory_wplaq(U,E,p);
    do_trajectory_bp(U,E,p);
    //do_trajectory_impr(U,E,p);
    //do_trajectory_clover(U,E,p);
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_action(const GaugeField<group>& U,const VectorField<Algebra<group>>& E,const parameters& p) {
    // measure the total action, consisting of plaquette and momentum term
    atype plaq=measure_s(U);
    atype e2=measure_e2(E);
    return p.beta*plaq+e2/2;
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_action(const GaugeField<group>& U,const VectorField<Algebra<group>>& E,const parameters& p, atype& plaq) {
    // measure the total action, consisting of plaquette and momentum term
    plaq=p.beta*measure_s(U);
    atype e2=measure_e2(E);
    return plaq+e2/2;
}

// end HMC action and force wrappers
///////////////////////////////////////////////////////////////////////////////////
// measurement functions

template <typename group>
void measure_stuff(const GaugeField<group>& U,const VectorField<Algebra<group>>& E) {
    // perform measurements on current gauge and momentum pair (U, E) and
    // print results in formatted form to standard output
    static bool first=true;
    if(first) {
        // print legend for measurement output
        hila::out0<<"LMEAS:     s-local          plaq           E^2        P.real        P.imag\n";
        first=false;
    }
    auto slocal=measure_s(U)/(lattice.volume()*NDIM*(NDIM-1)/2);
    auto plaq=measure_s_wplaq(U)/(lattice.volume()*NDIM*(NDIM-1)/2);
    auto e2=measure_e2(E)/(lattice.volume()*NDIM);
    auto poly=measure_polyakov(U,e_t);
    hila::out0<<string_format("MEAS % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e",slocal,plaq,e2,poly.real(),poly.imag())<<'\n';
}

// end measurement functions
///////////////////////////////////////////////////////////////////////////////////
// Wilson flow functions

template <typename group,typename atype=hila::arithmetic_type<group>>
void get_wf_force(const GaugeField<group>& U,VectorField<Algebra<group>>& E) {
    atype eps=2.0; // factor 2.0 to switch from unoriented to oriented plaquettes
                   // (factor is usually absorbed in \beta, but gradient flow force is
                   // computed from action term with \beta-factor stripped off) 

    // force computation routine to be used for wilson flow:

    //get_force_wplaq(U,E,eps); // force for Wilson plaquette action

    //get_force_bp(U,E,eps);    // force for bulk-preventing (BP) action (should only be used
                                // for gradient flow if BP is also used for HMC evolution)

    get_force_impr_f2(U,E,eps*c11,eps*c12); // force for improved action:  
                                            // 2.0 * ( c11*\sum_{unoriented plaquettes P} ReTr(P) 
                                            //        + c12*\sum_{unoriented 1x2-rectangles R}  ReTr(R) )
    //get_force_clover_action(U,E,eps); // force for clover action (not useful by itself)

    // 2/9*BP+7/9*Wilson (cancels O(a^4) terms)
    //get_force_bp(U,E,eps*2.0/9.0);
    //get_force_wplaq_add(U,E,eps*7.0/9.0);
};

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_dE_wplaq_dt(const GaugeField<group>& U) {
    Reduction<atype> declov=0;
    declov.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K,Kc;
    get_wf_force(U,K);
    get_force_wplaq(U,Kc,2.0);
    foralldir(d) {
        onsites(ALL) {
            declov+=Kc[d][X].dot(K[d][X]);
        }
    }
    return declov.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_dE_clov_dt(const GaugeField<group>& U) {
    Reduction<atype> declov=0;
    declov.allreduce(false).delayed(true);
    VectorField<Algebra<group>> K,Kc;
    get_wf_force(U,K);
    get_force_clover_action(U,Kc,2.0);
    foralldir(d) {
        onsites(ALL) {
            declov+=Kc[d][X].dot(K[d][X]);
        }
    }
    return declov.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void measure_wflow_stuff(const GaugeField<group>& V, atype flow_l, atype t_step) {
    // perform measurements on flowed gauge configuration V at flow scale flow_l
    // [t_step is the flow time integration step size used in last gradient flow step]
    // and print results in formatted form to standard output
    static bool first=true;
    if(first) {
        // print legend for flow measurement output
        hila::out0<<"LWFLMEAS  l(ambda)        S-flow        S-plaq        E_plaq    dE_plaq/dl         E_clv     dE_clv/dl     Qtopo_clv         E_log     Qtopo_log   [t step size]\n";
        first=false;
    }
    atype slocal=measure_s(V)/(lattice.volume()*NDIM*(NDIM-1)/2); // average action per plaquette 
    atype plaq=measure_s_wplaq(V)/(lattice.volume()*NDIM*(NDIM-1)/2); // average wilson plaquette action
    atype eplaq=plaq*NDIM*(NDIM-1)*group::size(); // naive energy density (based on wilson plaquette action)

    // average energy density and toplogical charge from 
    // symmetric log definition of field strength tensor :
    atype qtopolog,elog;
    measure_topo_charge_and_energy_log(V,qtopolog,elog);
    elog/=lattice.volume();

    // average energy density and toplogical charge from 
    // clover definition of field strength tensor :
    atype qtopocl,ecl;
    measure_topo_charge_and_energy_clover(V,qtopocl,ecl);
    ecl/=lattice.volume();

    // derivative of plaquette energy density w.r.t. to flow time :
    atype deplaqdt=measure_dE_wplaq_dt(V)/lattice.volume();

    // derivative of clover energy density w.r.t. to flow time :
    atype declovdt=measure_dE_clov_dt(V)/lattice.volume();
    
    // print formatted results to standard output :
    hila::out0<<string_format("WFLMEAS  % 9.3f % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e       [%0.5f]",flow_l,slocal,plaq,eplaq,0.25*flow_l*deplaqdt,ecl,0.25*flow_l*declovdt,qtopocl,elog,qtopolog,t_step)<<'\n';
}

// end Wilson flow functions
///////////////////////////////////////////////////////////////////////////////////
// load/save config functions

template <typename group>
void checkpoint(const GaugeField<group>& U,int trajectory,const parameters& p) {
    double t=hila::gettime();
    // name of config with extra suffix
    std::string config_file=p.config_file+"_"+std::to_string(((trajectory+1)/p.n_save)%2);
    // save config
    U.config_write(config_file);
    // write run_status file
    if(hila::myrank()==0) {
        std::ofstream outf;
        outf.open("run_status",std::ios::out|std::ios::trunc);
        outf<<"trajectory  "<<trajectory+1<<'\n';
        outf<<"seed        "<<static_cast<uint64_t>(hila::random()*(1UL<<61))<<'\n';
        outf<<"time        "<<hila::gettime()<<'\n';
        // write config name to status file:
        outf<<"config name  "<<config_file<<'\n';
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
        // get config name with suffix from status file:
        std::string config_file=status.get("config name");
        status.close();
        hila::seed_random(seed);
        U.config_read(config_file);
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
    // number of thermalization trajectories
    p.n_therm=par.get("thermalization trajs");
    // wilson flow frequency (number of traj. between w. flow measurement)
    p.wflow_freq=par.get("wflow freq");
    // wilson flow max. flow-distance
    p.wflow_max_l=par.get("wflow max lambda");
    // wilson flow flow-distance step size
    p.wflow_l_step=par.get("wflow lambda step");
    // wilson flow absolute accuracy (per integration step)
    p.wflow_a_accu=par.get("wflow abs. accuracy");
    // wilson flow relative accuracy (per integration step)
    p.wflow_r_accu=par.get("wflow rel. accuracy");
    // random seed = 0 -> get seed from time
    long seed=par.get("random seed");
    // save config and checkpoint
    p.n_save=par.get("trajs/saved");
    // measure surface properties and print "profile"
    p.config_file=par.get("config name");

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

    auto orig_dt=p.dt;
    auto orig_trajlen=p.trajlen;
    GaugeField<mygroup> U_old;
    int nreject=0;
    for(int trajectory=start_traj; trajectory<p.n_traj; ++trajectory) {
        if(trajectory<p.n_therm) {
            // during thermalization: start with 10% of normal step size (and trajectory length)
            // and increse linearly with number of accepted thermalization trajectories. normal
            // step size is reached after 3/4*p.n_therm accepted thermalization trajectories.
            if(trajectory<p.n_therm*3.0/4.0) {
                p.dt=orig_dt*(0.1+0.9*4.0/3.0*trajectory/p.n_therm);
                p.trajlen=orig_trajlen;
            } else {
                p.dt=orig_dt;
                p.trajlen=orig_trajlen;
            }
            if(nreject>1) {
                // if two consecutive thermalization trajectories are rejected, decrese the
                // step size by factor 0.5 (but keeping trajectory length constant, since
                // thermalization becomes inefficient if trajectory length decreases)
                for(int irej=0; irej<nreject-1; ++irej) {
                    p.dt*=0.5;
                    p.trajlen*=2;
                }
                hila::out0<<" thermalization step size(reduzed due to multiple reject) dt="<<std::setprecision(8)<<p.dt<<'\n';
            } else {
                hila::out0<<" thermalization step size dt="<<std::setprecision(8)<<p.dt<<'\n';
            }
        } else if(trajectory==p.n_therm) {
            p.dt=orig_dt;
            p.trajlen=orig_trajlen;
            hila::out0<<" normal stepsize dt="<<std::setprecision(8)<<p.dt<<'\n';
        }

        U_old=U;
        double ttime=hila::gettime();

        foralldir(d) onsites(ALL) E[d][X].gaussian_random();

        ftype g_act_old;
        ftype act_old=measure_action(U,E,p,g_act_old);

        do_trajectory(U,E,p);

        ftype g_act_new;
        ftype act_new=measure_action(U,E,p,g_act_new);


        bool reject=hila::broadcast(exp(act_old-act_new)<hila::random());
        
        if(trajectory<p.n_therm) {
            // during thermalization: keep track of number of rejected trajectories
            if(reject) {
                ++nreject;
                --trajectory;
            } else {
                if(nreject>0) {
                    --nreject;
                }
            }
        }

        hila::out0<<std::setprecision(12)<<"HMC "<<trajectory<<" S_TOT_start "<<act_old<<" dS_TOT "
            <<std::setprecision(6)<<act_new-act_old<<std::setprecision(12);
        if(reject) {
            hila::out0<<" REJECT"<<" --> S_GAUGE "<<g_act_old;
        } else {
            hila::out0<<" ACCEPT"<<" --> S_GAUGE "<<g_act_new;
        }

        hila::out0<<"  time "<<std::setprecision(3)<<hila::gettime()-ttime<<'\n';
        
        if(reject) {
            U=U_old;
        }

        hila::out0<<"Measure_start "<<trajectory<<'\n';

        measure_stuff(U,E);

        hila::out0<<"Measure_end "<<trajectory<<'\n';

        if(trajectory>=p.n_therm) {

            if(p.wflow_freq>0&&trajectory%p.wflow_freq==0) {
                int wtrajectory=trajectory/p.wflow_freq;
                if(p.wflow_l_step>0) {
                    int nflow_steps=(int)(p.wflow_max_l/p.wflow_l_step);

                    double wtime=hila::gettime();
                    hila::out0<<"Wflow_start "<<wtrajectory<<'\n';

                    GaugeField<mygroup> V=U;
                    ftype t_step=0.001;
                    for(int i=0; i<nflow_steps; ++i) {

                        t_step=do_wilson_flow_adapt(V,i*p.wflow_l_step,(i+1)*p.wflow_l_step,get_wf_force,p.wflow_a_accu,p.wflow_r_accu,t_step);

                        measure_wflow_stuff(V,(i+1)*p.wflow_l_step,t_step);

                    }

                    hila::out0<<"Wflow_end "<<wtrajectory<<"    time "<<std::setprecision(3)<<hila::gettime()-wtime<<'\n';
                }
            }

        }

        if(p.n_save>0&&(trajectory+1)%p.n_save==0) {
            checkpoint(U,trajectory,p);
        }
    }

    hila::finishrun();
}
