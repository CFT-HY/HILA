#include "hila.h"
#include "gauge/staples.h"
#include "gauge/polyakov.h"
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
    ftype beta;
    ftype dt;
    int trajlen;
    int n_traj;
    int n_therm;
    int wflow_freq;
    ftype wflow_max_l;
    ftype wflow_l_step;
    int n_save;
    std::string config_file;
    ftype time_offset;
};


///////////////////////////////////////////////////////////////////////////////////
// general functions

template <typename group>
ftype measure_plaq(const GaugeField<group>& U) {
    // measure the Wilson plaquette action
    Reduction<ftype> plaq=0;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            plaq+=1.0-real(trace(U[dir1][X]*U[dir2][X+dir1]*
                U[dir1][X+dir2].dagger()*U[dir2][X].dagger()))/group::size();
        }
    }
    return plaq.value();
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

template <typename group,typename atype=hila::arithmetic_type<group>>
void update_U(GaugeField<group>& U,const VectorField<Algebra<group>>& E,atype delta) {
    // evolve U with momentum E over time stepe delta 
    foralldir(d) {
        onsites(ALL) U[d][X]=chexp(E[d][X]*delta)*U[d][X];
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void measure_topo_charge_and_energy_log(const GaugeField<group>& U, atype& qtopo_out, atype& energy_out) {
    // measure topological charge and field strength energy of the gauge field, using the
    // matrix logarithms of the plaquettes as components of the field strength tensor

    Reduction<atype> qtopo=0;
    Reduction<atype> energy=0;
    qtopo.allreduce(false).delayed(true);

#if NDIM == 4
    Field<group> F[6];
    // F[0]: F[0][1], F[1]: F[0][2], F[2]: F[0][3],
    // F[3]: F[1][2], F[4]: F[1][3], F[5]: F[2][3]
    Field<group> tF0,tF1;

    int k=0;
    // (Note: by replacing log() by project_to_algebra() in the following, one would get the clover F[dir1][dir2])
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            // log of dir1-dir2-plaquette that starts and ends at X; corresponds to F[dir1][dir2]
            // at center location of plaquette X+dir1/2+dir2/2 :
            tF0[X]=log((U[dir1][X]*U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger())).expand();
            // parallel transport to X+dir1
            tF1[X]=U[dir1][X].dagger()*tF0[X]*U[dir1][X];
        }
        onsites(ALL) {
            // get F[dir1][dir2] at X from average of parallel transported F[dir1][dir2] from the 
            // centers of the dir1-dir2-plaquettes that surround X :
            F[k][X]=( tF0[X]+tF1[X-dir1] + U[dir2][X-dir2].dagger()*(tF0[X-dir2]+tF1[X-dir1-dir2])*U[dir2][X-dir2] )*0.25;
        }
        ++k;
    }
    onsites(ALL) {
        qtopo+=real(trace(F[0][X].dagger()*F[5][X]));
        qtopo+=-real(trace(F[1][X].dagger()*F[4][X]));
        qtopo+=real(trace(F[2][X].dagger()*F[3][X]));

        energy+=F[0][X].squarenorm();
        energy+=F[1][X].squarenorm();
        energy+=F[2][X].squarenorm();
        energy+=F[3][X].squarenorm();
        energy+=F[4][X].squarenorm();
        energy+=F[5][X].squarenorm();
    }
#endif
    qtopo_out=qtopo.value()/(4.0*M_PI*M_PI);
    energy_out=energy.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void measure_topo_charge_and_energy_clover(const GaugeField<group>& U,atype& qtopo_out,atype& energy_out) {
    // measure topological charge and field strength energy of the gauge field, using the
    // clover matrices as components of the field strength tensor

    Reduction<atype> qtopo=0;
    Reduction<atype> energy=0;
    qtopo.allreduce(false).delayed(true);

#if NDIM == 4
    Field<Algebra<group>> F[6];
    // F[0]: F[0][1], F[1]: F[0][2], F[2]: F[0][3],
    // F[3]: F[1][2], F[4]: F[1][3], F[5]: F[2][3]
    int k=0;
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            F[k][X]=0.25*(U[dir1][X]*U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger()
                -(U[dir1][X-dir1-dir2]*U[dir2][X-dir2]).dagger()*U[dir2][X-dir1-dir2]*U[dir1][X-dir1]
                +U[dir2][X]*(U[dir2][X-dir1]*U[dir1][X-dir1+dir2]).dagger()*U[dir1][X-dir1]
                -U[dir1][X]*(U[dir1][X-dir2]*U[dir2][X+dir1-dir2]).dagger()*U[dir2][X-dir2]).project_to_algebra();
        }
        ++k;
    }
    onsites(ALL) {
        qtopo+=real(trace(F[0][X].dagger()*F[5][X]));
        qtopo+=-real(trace(F[1][X].dagger()*F[4][X]));
        qtopo+=real(trace(F[2][X].dagger()*F[3][X]));

        energy+=F[0][X].squarenorm();
        energy+=F[1][X].squarenorm();
        energy+=F[2][X].squarenorm();
        energy+=F[3][X].squarenorm();
        energy+=F[4][X].squarenorm();
        energy+=F[5][X].squarenorm();
    }
#endif
    qtopo_out=qtopo.value()/(8.0*M_PI*M_PI);
    energy_out=energy.value()*0.5;
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void measure_topo_charge_and_energy_log_clover(const GaugeField<group>& U,atype& qtopo_out,atype& energy_out) {
    // measure topological charge and field strength energy of the gauge field, using the
    // logs of the "clover leave matrices" as components of the field strength tensor
    // (this is here just for cross check of log-definition above; should be the same, 
    // but more costly since log is taken four times more)

    Reduction<atype> qtopo=0;
    Reduction<atype> energy=0;
    qtopo.allreduce(false).delayed(true);

#if NDIM == 4
    Field<Algebra<group>> F[6];
    // F[0]: F[0][1], F[1]: F[0][2], F[2]: F[0][3],
    // F[3]: F[1][2], F[4]: F[1][3], F[5]: F[2][3]
    int k=0;
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            F[k][X]=log(U[dir1][X]*U[dir2][X+dir1]*(U[dir2][X]*U[dir1][X+dir2]).dagger());
            F[k][X]-=log((U[dir1][X-dir1-dir2]*U[dir2][X-dir2]).dagger()*U[dir2][X-dir1-dir2]*U[dir1][X-dir1]);
            F[k][X]+=log(U[dir2][X]*(U[dir2][X-dir1]*U[dir1][X-dir1+dir2]).dagger()*U[dir1][X-dir1]);
            F[k][X]-=log(U[dir1][X]*(U[dir1][X-dir2]*U[dir2][X+dir1-dir2]).dagger()*U[dir2][X-dir2]);
            F[k][X]*=0.25;
        }
        ++k;
    }
    onsites(ALL) {
        qtopo+=real(trace(F[0][X].dagger()*F[5][X]));
        qtopo+=-real(trace(F[1][X].dagger()*F[4][X]));
        qtopo+=real(trace(F[2][X].dagger()*F[3][X]));

        energy+=F[0][X].squarenorm();
        energy+=F[1][X].squarenorm();
        energy+=F[2][X].squarenorm();
        energy+=F[3][X].squarenorm();
        energy+=F[4][X].squarenorm();
        energy+=F[5][X].squarenorm();
    }
#endif
    qtopo_out=qtopo.value()/(8.0*M_PI*M_PI);
    energy_out=energy.value()*0.5;
}


// end general functions
///////////////////////////////////////////////////////////////////////////////////
// non-bulk-prevention functions

template <typename group>
void get_force(const GaugeField<group>& U,VectorField<Algebra<group>>& K) {
    // compute the force for the plaquette action and write it to K
    Field<group> staple;
    foralldir(d) {
        staplesum(U,staple,d);
        onsites(ALL) {
            K[d][X]=-(U[d][X]*staple[X].dagger()).project_to_algebra();
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void update_E(const GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p,atype delta) {
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

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_action(const GaugeField<group>& U,const VectorField<Algebra<group>>& E,const parameters& p) {
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
    // end trajectory: bring U and E to the same time
    update_E(U,E,p,p.dt);
    update_U(U,E,p.dt/2);
    
    U.reunitarize_gauge();
}

// end non-bulk-prevention functions
///////////////////////////////////////////////////////////////////////////////////
// bulk-prevention functions, cf. arXiv:2306.14319 (with n=2)

template <typename T>
T get_ch_inv(const T& U) {
    // compute inverse of the square matrix U, using the 
    // Cayley-Hamilton theorem (Faddeev–LeVerrier algorithm)
    T tB[2];
    int ip=0;
    tB[ip]=1.;
    auto tc=trace(U);
    tB[1-ip]=U;
    for(int k=2; k<=T::size(); ++k) {
        tB[1-ip]-=tc;
        tB[ip]=U*tB[1-ip];
        tc=trace(tB[ip])/k;
        ip=1-ip;
    }
    return tB[ip]/tc;
}

template <typename T>
T get_bp_UAmat(const T& U) {
    // compute U*A(U) with the A-matrix from Eq. (B3) of arXiv:2306.14319 for n=2
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
void get_force_bp(const GaugeField<group>& U,VectorField<Algebra<group>>& K) {
    // compute force for BP action for n=2 according to Eq. (B5) of arXiv:2306.14319 and write it to K
    Field<group> fmatp;
    Field<group> fmatmd1;
    Field<group> fmatmd2;
    auto eps=2.0;
    bool first=true;
    foralldir(d1) {
        onsites(ALL) K[d1][X]=0;
    }
    foralldir(d1) {
        foralldir(d2) if(d2>d1) {
            onsites(ALL) {
                fmatp[X]=get_bp_UAmat(U[d1][X]*U[d2][X+d1]*(U[d2][X]*U[d1][X+d2]).dagger());
                fmatmd1[X]=(fmatp[X]*U[d2][X]).dagger()*U[d2][X]; // parallel transport fmatp[X].dagger() to X+d2
                fmatmd2[X]=U[d1][X].dagger()*fmatp[X]*U[d1][X]; // parallel transport fmatp[X] to X+d1
            }
            onsites(ALL) {
                K[d1][X]-=eps*(fmatmd1[X-d2]+fmatp[X]).project_to_algebra();
                K[d2][X]-=eps*(fmatmd2[X-d1]-fmatp[X]).project_to_algebra();
            }
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void update_E_bp(const GaugeField<group>& U,VectorField<Algebra<group>>& E,const parameters& p,atype delta) {
    // compute force for BP action for n=2 according to Eq. (B5) of arXiv:2306.14319 and use it to evolve E
    Field<group> fmatp;
    Field<group> fmatmd1;
    Field<group> fmatmd2;
    auto eps=delta*2.0*p.beta/group::size();
    foralldir(d1) {
        foralldir(d2) if(d2>d1) {
            onsites(ALL) {
                fmatp[X]=get_bp_UAmat(U[d1][X]*U[d2][X+d1]*(U[d2][X]*U[d1][X+d2]).dagger());
                fmatmd1[X]=(fmatp[X]*U[d2][X]).dagger()*U[d2][X]; // parallel transport fmatp[X].dagger() to X+d2
                fmatmd2[X]=U[d1][X].dagger()*fmatp[X]*U[d1][X]; // parallel transport fmatp[X] to X+d1
            }
            onsites(ALL) {
                E[d1][X]-=eps*(fmatmd1[X-d2]+fmatp[X]).project_to_algebra();
                E[d2][X]-=eps*(fmatmd2[X-d1]-fmatp[X]).project_to_algebra();
            }
        }
    }
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_plaq_bp(const GaugeField<group>& U) {
    // measure the BP plaquette action
    Reduction<hila::arithmetic_type<group>> plaq=0;
    plaq.allreduce(false).delayed(true);
    foralldir(dir1) foralldir(dir2) if(dir1<dir2) {
        onsites(ALL) {
            plaq+=real(trace(get_bp_iOsqmat(U[dir1][X]*U[dir2][X+dir1]*
                (U[dir2][X]*U[dir1][X+dir2]).dagger())))/group::size();
        }
    }
    return plaq.value();
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype measure_action_bp(const GaugeField<group>& U,const VectorField<Algebra<group>>& E,const parameters& p) {
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
    // end trajectory: bring U and E to the same time
    update_E_bp(U,E,p,p.dt);
    update_U(U,E,p.dt/2);

    U.reunitarize_gauge();
}

// end bulk-prevention functions
///////////////////////////////////////////////////////////////////////////////////
// measurement functions

template <typename group>
void measure_stuff(const GaugeField<group>& U,const VectorField<Algebra<group>>& E) {
    // perform measurements on current gauge and momentum pair (U, E) and
    // print results in formatted form to standard output
    static bool first=true;
    if(first) {
        // print legend for measurement output
        hila::out0<<"Legend MEAS: BP-plaq  plaq  E^2  P.real  P.imag\n";
        first=false;
    }
    auto plaqbp=measure_plaq_bp(U)/(lattice.volume()*NDIM*(NDIM-1)/2);
    auto plaq=measure_plaq(U)/(lattice.volume()*NDIM*(NDIM-1)/2);
    auto e2=measure_e2(E)/(lattice.volume()*NDIM);
    auto poly=measure_polyakov(U,e_t);
    hila::out0<<string_format("MEAS % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e",plaqbp,plaq,e2,poly.real(),poly.imag())<<'\n';
}

template <typename group,typename atype=hila::arithmetic_type<group>>
void measure_wflow_stuff(const GaugeField<group>& V, atype flow_l, atype t_step) {
    // perform measurements on flowed gauge configuration V at flow scale flow_l
    // [t_step is the flow time integration step size used in last gradient flow step]
    // and print results in formatted form to standard output
    static bool first=true;
    if(first) {
        // print legend for flow measurement output
        hila::out0<<"Legend  WFLOW MEAS  f.-scale  BP-S-plaq  S-plaq  t^2*E_plaq  t^2*E_log  Qtopo_log  t^2*E_clov  Qtopo_clov  [t step size]\n";
        first=false;
    }
    auto plaqbp=measure_plaq_bp(V)/(lattice.volume()*NDIM*(NDIM-1)/2); // average bulk-prevention plaquette action
    auto eplaq=measure_plaq(V)*2.0/lattice.volume();
    auto plaq=eplaq/(NDIM*(NDIM-1)); // average wilson plaquette action
    eplaq*=group::size(); // average naive energy density (based on wilson plaquette action)

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

    // print formatted results to standard output :
    hila::out0<<string_format("WFLOW MEAS % 6.2f % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e % 0.6e   [%0.5f]",flow_l,plaqbp,plaq,pow(flow_l,4)/64.0*eplaq,pow(flow_l,4)/64.0*elog,qtopolog,pow(flow_l,4)/64.0*ecl,qtopocl,t_step)<<'\n';
}

template <typename group>
void get_wf_force(const GaugeField<group>& U,VectorField<Algebra<group>>& E) {
    // force computation routine to be used to integrate wilson flow
    get_force(U,E);
    //get_force_bp(U,E);
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype do_wilson_flow(GaugeField<group>& V,atype l_start,atype l_end) {
    // wilson flow integration from flow scale l_start to l_end
    // using 3rd order Runge-Kutta from arXiv:1006.4518 (cf. appendix C of
    // arXiv:2101.05320 for derivation of this Runge-Kutta method)
    atype tmin=l_start*l_start/8.0;
    atype tmax=l_end*l_end/8.0;
    int nsteps=10;
    atype step=(tmax-tmin)/(atype)nsteps;
    while(step>0.08) {
        //make sure the step size is not too large
        ++nsteps;
        step=(tmax-tmin)/(atype)nsteps;
    }

    VectorField<Algebra<group>> k1,k2;

    // 3rd order Runge-Kutta coefficients from arXiv:1006.4518 :
    atype a11=0.25;
    atype a21=-17.0/36.0,a22=8.0/9.0;
    atype a33=0.75;

    for(int i=0; i<nsteps; ++i) {
        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            V[d][X]=chexp(k1[d][X]*(step*a11))*V[d][X];
        }

        get_wf_force(V,k2);
        foralldir(d) onsites(ALL) {
            k2[d][X]=k2[d][X]*(step*a22)+k1[d][X]*(step*a21);
            V[d][X]=chexp(k2[d][X])*V[d][X];
        }

        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            k1[d][X]=k1[d][X]*(step*a33)-k2[d][X];
            V[d][X]=chexp(k1[d][X])*V[d][X];
        }
    }

    V.reunitarize_gauge();   
    return step;
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype do_wilson_flow_adapt_b(GaugeField<group>& V,atype l_start,atype l_end,atype tstep=0.08, atype atol=1.0e-7, atype rtol=1.0e-7) {
    // wilson flow integration from flow scale l_start to l_end
    // using 3rd order Runge-Kutta from arXiv:1006.4518 (cf. appendix C of
    // arXiv:2101.05320 for derivation of this Runge-Kutta method)
    atype t=l_start*l_start/8.0;
    atype tmax=l_end*l_end/8.0;
    int nsteps=10;
    atype step=(tmax-t)/(atype)nsteps;
    while(step>0.08) {
        //make sure the step size is not too large
        ++nsteps;
        step=(tmax-t)/(atype)nsteps;
    }

    VectorField<Algebra<group>> k1,k2;
    GaugeField<group> V2,V0;
    Field<atype> reldiff;
    atype maxreldiff,maxstep;

    // 3rd order Runge-Kutta coefficients from arXiv:1006.4518 :
    atype a11=0.25;
    atype a21=-17.0/36.0,a22=8.0/9.0;
    atype a33=0.75;

    // 2nd order step coefficients :
    atype b21=-1.25,b22=2.0;
    int i=0;
    V0=V;
    while(i<nsteps&&t<tmax) {

        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            V[d][X]=chexp(k1[d][X]*(step*a11))*V[d][X];
        }

        get_wf_force(V,k2);
        foralldir(d) onsites(ALL) {
            V2[d][X]=chexp(k2[d][X]*(step*b22)+k1[d][X]*(step*b21))*V[d][X];

            k2[d][X]=k2[d][X]*(step*a22)+k1[d][X]*(step*a21);
            V[d][X]=chexp(k2[d][X])*V[d][X];
        }

        onsites(ALL) {
            reldiff[X]=0;
        }

        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            k1[d][X]=k1[d][X]*(step*a33)-k2[d][X];
            V[d][X]=chexp(k1[d][X])*V[d][X];

            reldiff[X]+=(V[d][X]-V2[d][X]).norm()/(atol+rtol*V[d][X].norm());
        }

        maxreldiff=reldiff.max()/(NDIM*group::size()*group::size());

        maxstep=min(step/pow(maxreldiff,1.0/3.0),0.8);
        if(step>maxstep) {
            while(step>maxstep) {
                //make sure the step size is not too large
                ++nsteps;
                step=(tmax-t)/(atype)(nsteps-i);
            }
            V=V0;
        } else {
            t+=step;
            ++i;
            while(step<0.5*maxstep&&nsteps>i+1) {
                //make sure the step size is not too small
                --nsteps;
                step=(tmax-t)/(atype)(nsteps-i);
            }
            V0=V;
        }
        /*
        if(hila::myrank()==0) {
            std::cout<<"wflow  i: "<<i<<"  nsteps: "<<nsteps<<"  step: "<<step<<"  maxstep: "<<maxstep<<"  maxreldiff: "<<maxreldiff<<std::endl;
        }
        */
    }

    V.reunitarize_gauge();
    return step;
}

template <typename group,typename atype=hila::arithmetic_type<group>>
atype do_wilson_flow_adapt(GaugeField<group>& V,atype l_start,atype l_end,atype tstep=0.001,atype atol=1.0e-7,atype rtol=1.0e-7) {
    // wilson flow integration from flow scale l_start to l_end
    // using 3rd order Runge-Kutta from arXiv:1006.4518 (cf. appendix C of
    // arXiv:2101.05320 for derivation of this Runge-Kutta method)

    // translate flow scale interval [l_start,l_end] to corresponding
    // flow time interval [t,tmax] :
    atype t=l_start*l_start/8.0;
    atype tmax=l_end*l_end/8.0;
    atype step=min(min(tstep,0.51*(tmax-t)),0.08);  //initial step size

    VectorField<Algebra<group>> k1,k2;
    GaugeField<group> V2,V0;
    Field<atype> reldiff;
    atype maxreldiff,maxstep;

    // 3rd order Runge-Kutta coefficients from arXiv:1006.4518 :
    atype a11=0.25;
    atype a21=-17.0/36.0,a22=8.0/9.0;
    atype a33=0.75;

    // 2nd order step coefficients :
    atype b21=-1.25,b22=2.0;
    V0=V;
    bool stop=false;

    while(t<tmax&&!stop) {

        if(t+step>=tmax) {
            step=tmax-t;
            stop=true;
        } else {
            if(t+2.0*step>=tmax) {
                step=0.51*(tmax-t);
            }
        }

        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            // first steps of 3rd and 2nd order RK are the same :
            V[d][X]=chexp(k1[d][X]*(step*a11))*V[d][X];
        }

        get_wf_force(V,k2);
        foralldir(d) onsites(ALL) {
            // second step for 2nd order RK :
            V2[d][X]=chexp(k2[d][X]*(step*b22)+k1[d][X]*(step*b21))*V[d][X];

            // second step for 3rd order RK :
            k2[d][X]=k2[d][X]*(step*a22)+k1[d][X]*(step*a21);
            V[d][X]=chexp(k2[d][X])*V[d][X];
        }

        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            // third step of 3rd order RK :
            k1[d][X]=k1[d][X]*(step*a33)-k2[d][X];
            V[d][X]=chexp(k1[d][X])*V[d][X];
        }

        // determine maximum difference between 2nd and 3rd order RK, 
        // relative to desired accuracy :
        onsites(ALL) {
            reldiff[X]=0;
        }
        foralldir(d) onsites(ALL) {
            reldiff[X]+=(V[d][X]-V2[d][X]).norm()/(atol+rtol*V[d][X].norm());
        }
        maxreldiff=reldiff.max()/(NDIM*group::size()*group::size());

        // max. allowed step size to achieve desired accuracy :
        maxstep=min(step/pow(maxreldiff,1.0/3.0),1.0);
        if(step>maxstep) {
            // repeat current iteration if step size was larger than maxstep
            V=V0;
        } else {
            // proceed to next iteration
            t+=step;
            V0=V;
        }
        // adjust step size for next iteration to better match accuracy goal :
        step=min(0.9*maxstep,0.08); 
    }

    V.reunitarize_gauge();
    return step;
}


template <typename group,typename atype=hila::arithmetic_type<group>>
atype do_wilson_flow_o2(GaugeField<group>& V,atype l_start,atype l_end) {
    // wilson flow integration from flow scale l_start to l_end
    // using 2rd order Runge-Kutta (cf. arXiv:2101.05320 for description
    // of method to derivation coefficients)
    atype tmin=l_start*l_start/8.0;
    atype tmax=l_end*l_end/8.0;
    int nsteps=10;
    atype step=(tmax-tmin)/(atype)nsteps;
    while(step>0.08) {
        //make sure the step size is not too large
        ++nsteps;
        step=(tmax-tmin)/(atype)nsteps;
    }

    VectorField<Algebra<group>> k1,k2;

    // 2nd order Runge-Kutta coefficients :
    atype b11=0.25;
    atype b21=-1.25,b22=2.0;

    for(int i=0; i<nsteps; ++i) {
        get_wf_force(V,k1);
        foralldir(d) onsites(ALL) {
            V[d][X]=chexp(k1[d][X]*(step*b11))*V[d][X];
        }

        get_wf_force(V,k2);
        foralldir(d) onsites(ALL) {
            V[d][X]=chexp(k2[d][X]*(step*b22)+k1[d][X]*(step*b21))*V[d][X];
        }

    }

    V.reunitarize_gauge();
    return step;
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
    // number of thermalization trajectories
    p.n_therm=par.get("thermalization trajs");
    // wilson flow frequency (number of traj. between w. flow measurement)
    p.wflow_freq=par.get("wflow freq");
    // wilson flow max. flow distance
    p.wflow_max_l=par.get("wflow max lambda");
    // wilson flow flow distance step size
    p.wflow_l_step=par.get("wflow lambda step");
    // random seed = 0 -> get seed from time
    long seed=par.get("random seed");
    // save config and checkpoint
    p.n_save=par.get("trajs/saved");
    // measure surface properties and print "profile"
    p.config_file=par.get("config name");

    par.close(); // file is closed also when par goes out of scope

    /*
    if(hila::myrank()==0) {
        // test matrix logarithm:
        mygroup tU;
        tU.random();
        std::cout<<"tU:"<<std::endl;
        std::cout<<tU<<std::endl;
        Algebra<mygroup> ltU=log(tU);
        mygroup ltUg=ltU.expand();
        std::cout<<"ltU:"<<std::endl;
        std::cout<<ltUg<<std::endl;

        Algebra<mygroup> tltU=ltU;
        std::cout<<tltU.dagger()*ltU/2<<" -- "<<real(trace(ltUg.dagger()*ltUg))<<std::endl;
        ltUg*=0.5;
        tU=chexp(ltUg);
        std::cout<<"sqrt(tU)=exp(ltU/2):"<<std::endl;
        std::cout<<tU<<std::endl;
        tU=tU*tU;
        std::cout<<"sqrt(tU)^2:"<<std::endl;
        std::cout<<tU<<std::endl;
    }
    */

    /*
    levi_civita<NDIM> LCS;
    if(hila::myrank()==0) {
        // check levi_citiva output:
        for(int i=0; i<LCS.nfac; ++i) {
            std::cout<<i<<"  :  ";
            for(int j=0; j<LCS.n; ++j) {
                std::cout<<LCS.a[i][j]<<" ";
            }
            std::cout<<"  :  "<<LCS.a[i][LCS.n]<<std::endl;
        }
    }
    */
    
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
    for(int trajectory=start_traj; trajectory<p.n_traj; trajectory++) {
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

        auto act_old=measure_action_bp(U,E,p);

        do_trajectory_bp(U,E,p);

        auto act_new=measure_action_bp(U,E,p);


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

                        t_step=do_wilson_flow_adapt(V,i*p.wflow_l_step,(i+1)*p.wflow_l_step,t_step);

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
