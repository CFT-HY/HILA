#include "hila.h"
#include "plaquettefield.h"
#include "tools/string_format.h"
#include "tools/floating_point_epsilon.h"


using ftype = double;

using mygroup = int;

#ifndef PLAQ_SHIFT
#define PLAQ_SHIFT 1
#endif

#ifndef PARITY
#define PARITY 1
#endif

::Parity uparity(const CoordinateVector &x) {
    int s = 0;
#if PARITY == 0
    // parity of x
    for (Direction d = e_x; d < NDIM; ++d) {
        s += x.e(d);
    }
#else
    // parity of spatial part of x
    for (Direction d = e_x; d < NDIM - 1; ++d) {
        s += x.e(d);
    }
#endif
    if (s % 2 == 0) {
        return Parity::even;
    } else {
        return Parity::odd;
    }
}

template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
T umod(T val, T per) {
    int nrval = abs(val) / per;
    if(val < 0) {
        return val + nrval * per;
    } else {
        return val - nrval * per;
    }
}

// define a struct to hold the input parameters: this
// makes it simpler to pass the values around
struct parameters {
    ftype beta;         // inverse gauge coupling
    int n_traj;         // number of trajectories to generate
    int n_therm;        // number of thermalization trajectories (counts only accepted traj.)
    int n_update;       // number of heat-bath sweeps per "trajectory"
    int n_or_update;    // number of overrelaxation updates per "trajectory"
    int n_ps_update; // number of plaquette shift (s_{x,\mu\nu}-vars) updates per "trajectory"
    int n_save; // number of trajectories between config. check point
    std::string config_file;
    ftype time_offset;
};


///////////////////////////////////////////////////////////////////////////////////
// Metropolis update functions

//template <typename T>
//using sw_t = std::array<std::array<Field<T>, NDIM>, NDIM>;
template <typename T>
using sw_t = PlaquetteField<T>;
/**
 * @brief Sum the staples of link variables to direction dir taking into account plaquette
 * orientations and shift weights
 *
 * action of all plaquettes containing h_{x,\mu}:
 * \sum_{\nu!=\mu} {(h_{x,\mu} + h_{x+\hat{\mu},\nu} - h_{x+\hat{\nu},\mu} - h_{x,\nu} + s_{x,\mu\nu})^2
 * + (h_{x-\hat{\nu},\mu} + h_{x-\hat{\nu}+\hat{\mu},\nu} - h_{x,\mu} - h_{x-\hat{\nu},\nu} +
 * s_{x-\hat{\nu},\mu\nu})^2}
 * the staple sum for h_{x,\mu} is defined by the sum of terms in the above expression which are
 * linear in h_{x,\mu}:
 * 2 (d-1) h_{x,\mu}^2
 *   + 2 h_{x,\mu} \sum_{\nu} {
 *              (h_{x+\hat{\mu},\nu} - h_{x+\hat{\nu},\mu} - h_{x,\nu}
 *               + s_{x,\mu\nu})
 *            + (-h_{x-\hat{\nu}+\hat{\mu},\nu} - h_{x-\hat{\nu},mu} + h_{x-\hat{\nu},\nu}
 *               - s_{x-\hat{\nu},\mu\nu})
 *     }
 *   + terms independent of h_{x,\mu}
 *
 * @tparam T Z-link group type
 * @tparam fT plaquette shift and staplesum type
 * @param H GaugeField to compute staples for
 * @param staplesum Filed to compute staplesum into at each lattice point
 * @param d1 Direction to compute staplesum for
 * @param sw plaquette shift
 * @param par Parity to compute staplesum for
 */
template <typename T, typename fT>
void staplesum(const GaugeField<T> &H, Field<fT> &staples, Direction d1,
               const sw_t<fT> &sw, Parity par = ALL) {

    Field<fT> lower;

    bool first = true;
    foralldir(d2) if (d2 != d1) {

        // anticipate that these are needed
        // not really necessary, but may be faster
        H[d2].start_gather(d1, ALL);
        H[d1].start_gather(d2, par);

        // calculate first lower 'U' of the staple sum
        // do it on opp parity
        onsites(opp_parity(par)) {
            lower[X] = ((fT)(-H[d2][X + d1] - H[d1][X] + H[d2][X]) - sw[d1][d2][X]);
        }

        // calculate then the upper 'n', and add the lower
        // lower could also be added on a separate loop
        if (first) {
            onsites(par) {
                staples[X] =
                    2.0 * ((fT)(H[d2][X + d1] - H[d1][X + d2] - H[d2][X]) + sw[d1][d2][X] + lower[X - d2]);
            }
            first = false;
        } else {
            onsites(par) {
                staples[X] +=
                    2.0 * ((fT)(H[d2][X + d1] - H[d1][X + d2] - H[d2][X]) + sw[d1][d2][X] + lower[X - d2]);
            }
        }
    }
}


/**
 * @brief Z-link theory metropolis update
 * @details --
 * @tparam T Group element type such as long or int
 * @tparam fT staple type such as double or float
 * @param h link variable to be updated
 * @param stapsum staplesum of plaquettes containing h
 * @param beta coupling constant
 * @return double change in plaquette action
 */
template <typename T, typename fT>
fT z_metropolis(T &h, const fT &stapsum, double beta) {
    fT nstap = (fT)(NDIM - 1) * 2.0; 
    fT si = nstap * h * h + (fT)h * stapsum;
    int he = h + (1 - 2 * (int)(hila::random() * 2.0));
    fT nds = nstap * he * he + (fT)he * stapsum - si;
    if (nds < 0 || hila::random() < exp(-0.5 * beta * nds)) {
        h = he;
        return -nds;
    } else {
        return 0;
    }
}

/**
 * @brief Z-link theory quasi-overrelax update
 * @details --
 * @tparam T Group element type such as long or int
 * @tparam fT staple type such as double or float
 * @param h link variable to be updated
 * @param stapsum staplesum of plaquettes containing h
 * @param beta coupling constant
 * @return double change in plaquette action
 */
template <typename T, typename fT>
fT z_overrelax(T &h, const fT &stapsum, double beta) {
    fT nstap = (fT)(NDIM - 1) * 2.0;
    fT si = nstap * h * h + (fT)h * stapsum;
    fT se = si;
    fT tsqrt = sqrt(4.0 * nstap * si + stapsum * stapsum);
    int he = (int)lround((tsqrt - stapsum) / (2.0 * nstap));
    if (he != h) {
        se = nstap * he * he + (fT)he * stapsum;
        tsqrt = sqrt(4.0 * nstap * se + stapsum * stapsum);
        if((int)lround((-tsqrt - stapsum) / (2.0 * nstap)) != h) {
            // perform only invertible overrelaxation updates
            return 0;
        }
    } else {
        he = (int)lround((-tsqrt - stapsum) / (2.0 * nstap));
        se = nstap * he * he + (fT)he * stapsum;
        tsqrt = sqrt(4.0 * nstap * se + stapsum * stapsum);
        if ((int)lround((tsqrt - stapsum) / (2.0 * nstap)) != h) {
            // perform only invertible overrelaxation updates
            return 0;
        }
    }
    fT nds = se - si;
    if (nds < 0 || hila::random() < exp(-0.5 * beta * nds)) {
        h = he;
        return -nds;
    } else {
        return 0;
    }
}

/**
 * @brief plaquette shift metropolis update
 * @details --
 * @tparam fT shift type such as double or float
 * @tparam T Group element type such as long or int
 * @param sw plaquette shift variable
 * @param tplaq plaquette variable
 * @param beta coupling constant
 * @return double change in plaquette action
 */
template <typename fT, typename T>
fT sw_metropolis(fT &sw,const T &tplaq, double beta) {
    fT si = ((fT)tplaq + sw) * ((fT)tplaq + sw);
    fT nds = ((fT)tplaq - sw) * ((fT)tplaq - sw) - si;
    if (nds < 0 || hila::random() < exp(-0.5 * beta * nds)) {
        sw = -sw;
        return -nds;
    } else {
        return 0;
    }
}

/**
 * @brief Wrapper function to update GaugeField per direction
 * @details Computes first staplesum, then uses computed result to evolve GaugeField with Metropolis
 * updates
 *
 * @tparam T Z-link group type
 * @tparam fT plaquette shift type
 * @param H GaugeField to evolve
 * @param p parameter struct
 * @param par Parity
 * @param d Direction to evolve
 * @param sw plaquette shifs
 */
template <typename T, typename fT>
void update_parity_dir(GaugeField<T> &H, const parameters &p, Parity par, Direction d, const sw_t<fT> &sw) {

    static hila::timer me_timer("Metropolis (z)");
    static hila::timer staples_timer("Staplesum");

    Field<fT> staples;

    staples_timer.start();

    staplesum(H, staples, d, sw, par);

    staples_timer.stop();

    me_timer.start();
    onsites(par) {
        z_metropolis(H[d][X], staples[X], p.beta);
    }
    me_timer.stop();

}

/**
 * @brief Wrapper function to update GaugeField per direction with overrelaxation
 * @details Computes first staplesum, then uses computed result to evolve GaugeField with
 * overrelaxation updates
 *
 * @tparam T Z-link group type
 * @tparam fT plaquette shift type
 * @param H GaugeField to evolve
 * @param p parameter struct
 * @param par Parity
 * @param d Direction to evolve
 * @param sw plaquette shifs
 */
template <typename T, typename fT>
void update_or_parity_dir(GaugeField<T> &H, const parameters &p, Parity par, Direction d,
                       const sw_t<fT> &sw) {

    static hila::timer or_timer("Overrelax (z)");
    static hila::timer staples_timer("Staplesum");

    Field<fT> staples;

    staples_timer.start();

    staplesum(H, staples, d, sw, par);

    staples_timer.stop();

    or_timer.start();
    onsites(par) {
        z_overrelax(H[d][X], staples[X], p.beta);
    }
    or_timer.stop();
}

/**
 * @brief Wrapper function to update plaquette shift variables
 * @details --
 *
 * @tparam fT shift type such as double or float
 * @tparam T Group element type such as long or int
 * @param sw plaquette shift variable
 * @param p paramaters
 * @param H GaugeField from which to compute the plaquettes
 * @return double change in plaquette action
 */
template <typename fT, typename T>
void update_sw_parity_dir(sw_t<fT> &sw, const parameters &p, Parity par, Direction d1,
                          const GaugeField<T> &H) {
    static hila::timer sw_timer("Metropolis (ps)");

    sw_timer.start();

    foralldir(d2) if (d2 != d1) {
        H[d2].start_gather(d1, par);
        H[d1].start_gather(d2, par);
        onsites(par) {
            T tplaq = H[d1][X] + H[d2][X + d1] - H[d1][X + d2] - H[d2][X];
            sw_metropolis(sw[d1][d2][X], tplaq, p.beta);
            sw[d2][d1][X] = -sw[d1][d2][X];
        }
    }

    sw_timer.stop();
}

/**
 * @brief Wrapper update function
 * @details Gauge Field update sweep with randomly chosen parities and directions
 *
 * @tparam T Z-gauge group type
 * @tparam fT plaquette shift type
 * @param H GaugeField to update
 * @param p Parameter struct
 * @param sw plaquette shifts
 */
template <typename T, typename fT>
void update(GaugeField<T> &H, sw_t<fT> &sw, const parameters &p) {

    for (int i = 0; i < 2 * NDIM; ++i) {
        int ud_type =
            hila::broadcast((int)(hila::random() * (p.n_update + p.n_or_update + p.n_ps_update)));

        int tdp = hila::broadcast((int)(hila::random() * 2 * NDIM));
        int tdir = tdp / 2;
        int tpar = 1 + (tdp % 2);
        // hila::out0 << "   " << Parity(tpar) << " -- " << Direction(tdir);
        if(ud_type < p.n_update) {
            // update Z-link variable on links (x,\mu) = (x_{tpar},tdir) with metropolis
            update_parity_dir(H, p, Parity(tpar), Direction(tdir), sw);
        } else if (ud_type < p.n_update + p.n_or_update) {
            // update Z-link variable on links (x,\mu) = (x_{tpar},tdir) with overrelaxation
            update_or_parity_dir(H, p, Parity(tpar), Direction(tdir), sw);
        } else {
            // update plaquette shift variables for plaquettes spanned by (x,\mu\nu) =
            // (x_{tpar},tdir,d2) for all d2!=tdir
            update_sw_parity_dir(sw, p, Parity(tpar), Direction(tdir), H);
        }
    }
    // hila::out0 << "\n";
}

/**
 * @brief apply site shift
 * @details Gauge Field is moved by one site
 *
 * @tparam T Z-gauge group type
 * @param H GaugeField to update
 */
template <typename T>
void site_shift(GaugeField<T> &H) {
#if PLAQ_SHIFT==2 && PARITY==1
    int tdp = hila::broadcast((int)(hila::random() * 2 * (NDIM - 1)));
    int tdir = tdp / 2;
    int tsgn = tdp % 2;
    Direction d = Direction(tdir);
    if (tsgn == 1) {
        d = opp_dir(d);
    }
    GaugeField<T> tH = H;
    foralldir(d1) {
        onsites(ALL) {
            H[d1][X] = tH[d1][X + d];
        }
    }
    onsites(ALL) {
        if (tsgn == 0) {
            if (uparity(X.coordinates()) == EVEN) {
                H[e_t][X] += 1;
            }
        } else {
            if (uparity(X.coordinates()) == ODD) {
                H[e_t][X] -= 1;
            }
        }
    }
#endif
}

/**
 * @brief Evolve Z-link field (and plaquette shift field, if set to be dynamic)
 * @details Evolution happens by means of metropolis updates.
 *
 * @tparam T Z-link group type
 * @tparam fT plaquette shift type
 * @param H GaugeField to evolve
 * @param sw plaquette shifs
 * @param p parameter struct
 */
template <typename T, typename fT>
void do_trajectory(GaugeField<T> &H, sw_t<fT> &sw, const parameters &p) {
    for (int n = 0; n < p.n_update + p.n_or_update + p.n_ps_update; n++) {
        update(H, sw, p);
    }
    
    site_shift(H);

    // value of the action is invariant under direction-dependent global Z-shifts.
    // to prevent the values from diverging, we shift them so that the link-variables
    // emerging in positive directions from the site at the origin are zero:
    CoordinateVector v = {0, 0, 0, 0};
    foralldir(d1) {
        auto dHpd = -H[d1][v];
        if (dHpd != 0) {
            onsites(ALL) {
                H[d1][X] += dHpd;
            }
        }
    }
}

// Metropolis update functions
///////////////////////////////////////////////////////////////////////////////////
// measurement functions


template <typename T, typename fT>
double measure_s_plaq(const GaugeField<T> &H, const sw_t<fT> &sw) {
    // measure the total plaquette action for the link field H (taking into account the plaquette
    // shifts)
    Reduction<double> plaq = 0;
    plaq.allreduce(false).delayed(true);
    foralldir(d1) foralldir(d2) if (d1 < d2) {
        H[d2].start_gather(d1, ALL);
        H[d1].start_gather(d2, ALL);
        onsites(ALL) {
            plaq += pow(((double)(H[d1][X] + H[d2][X + d1] - H[d1][X + d2] - H[d2][X]) +
                              (double)sw[d1][d2][X]),
                             2.0);
        }
    }
    return plaq.value();
}

template <typename T>
double measure_ns_plaq(const GaugeField<T> &H) {
    // measure the total plaquette action for the link field H (without plaquette shifts)
    Reduction<double> plaq = 0;
    plaq.allreduce(false).delayed(true);
    foralldir(d1) foralldir(d2) if (d1 < d2) {
        H[d2].start_gather(d1, ALL);
        H[d1].start_gather(d2, ALL);
        onsites(ALL) {
            plaq += pow((double)(H[d1][X] + H[d2][X + d1] - H[d1][X + d2] - H[d2][X]),
                        2.0);
        }
    }
    return plaq.value();
}




template <typename T>
double measure_hsq(const GaugeField<T> &H, double(out_only &hsq_per_par_dir)[2][NDIM],
                   double(out_only &h_per_par_dir)[2][NDIM]) {
    // measure the average magnitude squared of H and its average value per direction and parity
    Reduction<double> thsq = 0;
    thsq.allreduce(false).delayed(true);
    ReductionVector<double> h_per_p_d(2 * NDIM);
    h_per_p_d = 0.0;
    h_per_p_d.allreduce(false).delayed(true);
    ReductionVector<double> hsq_per_p_d(2 * NDIM);
    hsq_per_p_d = 0.0;
    hsq_per_p_d.allreduce(false).delayed(true);
    foralldir(d1) {
        onsites(ALL) {
            double tthsq = pow(H[d1][X], 2.0);
            thsq += tthsq;
            int tpar = (int)uparity(X.coordinates()) - 1;
            h_per_p_d[tpar * NDIM + d1] += H[d1][X];
            hsq_per_p_d[tpar * NDIM + d1] += tthsq;
        }
    }
    h_per_p_d.reduce();
    for (int par = 0; par < 2; ++par) {
        for (int dir = 0; dir < NDIM; ++dir) {
            h_per_par_dir[par][dir] = h_per_p_d[par * NDIM + dir] * 2.0 / lattice.volume();
        }
    }
    hsq_per_p_d.reduce();
    for (int par = 0; par < 2; ++par) {
        for (int dir = 0; dir < NDIM; ++dir) {
            hsq_per_par_dir[par][dir] = hsq_per_p_d[par * NDIM + dir] * 2.0 / lattice.volume();
        }
    }
    return thsq.value() / (lattice.volume() * NDIM);
}

template <typename T>
void measure_plaq_per_par_and_plane(const sw_t<T> &plaq, double(out_only &plaq_per_par_pl)[2][NDIM][NDIM]) {
    // measure the average plaquette value per plane and parity
    ReductionVector<double> h_per_p_pl(2 * NDIM * NDIM);
    h_per_p_pl = 0.0;
    h_per_p_pl.allreduce(false).delayed(true);
    foralldir(d1) foralldir(d2) if(d1 < d2) {
        onsites(ALL) {
            int tpar = (int)uparity(X.coordinates()) - 1;
            h_per_p_pl[(tpar * NDIM + d1) * NDIM + d2] +=
                (double)plaq[d1][d2][X];
        }
    }
    h_per_p_pl.reduce();
    for (int par = 0; par < 2; ++par) {
        for (int d1 = 0; d1 < NDIM; ++d1) {
            plaq_per_par_pl[par][d1][d1] = 0;
            for (int d2 = d1 + 1; d2 < NDIM; ++d2) {
                plaq_per_par_pl[par][d1][d2] =
                    h_per_p_pl[(par * NDIM + d1) * NDIM + d2] * 2.0 / lattice.volume();
                plaq_per_par_pl[par][d2][d1] = -plaq_per_par_pl[par][d1][d2];
            }
        }
    }
}

template <typename T>
void measure_plaq_sign_per_par_and_plane(const sw_t<T> &plaq,
                                    double(out_only &plaq_per_par_pl)[2][NDIM][NDIM]) {
    // measure the average plaquette value per plane and parity
    ReductionVector<double> h_per_p_pl(2 * NDIM * NDIM);
    h_per_p_pl = 0.0;
    h_per_p_pl.allreduce(false).delayed(true);
    foralldir(d1) foralldir(d2) if (d1 < d2) {
        onsites(ALL) {
            int tpar = (int)uparity(X.coordinates()) - 1;
            h_per_p_pl[(tpar * NDIM + d1) * NDIM + d2] += (double)sgn(plaq[d1][d2][X]);
        }
    }
    h_per_p_pl.reduce();
    for (int par = 0; par < 2; ++par) {
        for (int d1 = 0; d1 < NDIM; ++d1) {
            plaq_per_par_pl[par][d1][d1] = 0;
            for (int d2 = d1 + 1; d2 < NDIM; ++d2) {
                plaq_per_par_pl[par][d1][d2] =
                    h_per_p_pl[(par * NDIM + d1) * NDIM + d2] * 2.0 / lattice.volume();
                plaq_per_par_pl[par][d2][d1] = -plaq_per_par_pl[par][d1][d2];
            }
        }
    }
}

template <typename T, typename fT>
void measure_splaq_per_par_and_plane(const sw_t<T> &plaq, const sw_t<fT> &sw,
                                     double(out_only &plaq_per_par_pl)[2][NDIM][NDIM]) {
    // measure the average plaquette action (with and without shift) value per plane and parity
    ReductionVector<double> h_per_p_pl(2 * NDIM * NDIM);
    h_per_p_pl = 0.0;
    h_per_p_pl.allreduce(false).delayed(true);
    foralldir(d1) foralldir(d2) if (d1 < d2) {
        onsites(ALL) {
            int tpar = (int)uparity(X.coordinates()) - 1;
            h_per_p_pl[(tpar * NDIM + d1) * NDIM + d2] +=
                pow((double)plaq[d1][d2][X] + sw[d1][d2][X], 2.0);
            h_per_p_pl[(tpar * NDIM + d2) * NDIM + d1] += pow((double)plaq[d1][d2][X], 2.0);
        }
    }
    h_per_p_pl.reduce();
    for (int par = 0; par < 2; ++par) {
        for (int d1 = 0; d1 < NDIM; ++d1) {
            plaq_per_par_pl[par][d1][d1] = 0;
            for (int d2 = 0; d2 < NDIM; ++d2) {
                plaq_per_par_pl[par][d1][d2] =
                    h_per_p_pl[(par * NDIM + d1) * NDIM + d2] * 2.0 / lattice.volume();
            }
        }
    }
}

/**
 * @brief compute the H-plaquette field from the H-field
 * @tparam T Z-link group type
 * @param H Z-link field
 * @param plaq GaugeField[NDIM][NDIM] plaquette field
 */
template <typename T>
void plaq_field(const GaugeField<T> &H, out_only sw_t<T> &plaq) {
    foralldir(d1) {
        onsites(ALL) plaq[d1][d1][X] = 0;
        foralldir(d2) if (d1 < d2) {
            H[d2].start_gather(d1, ALL);
            H[d1].start_gather(d2, ALL);
            onsites(ALL) {
                plaq[d1][d2][X] = (H[d1][X] + H[d2][X + d1] - H[d1][X + d2] - H[d2][X]);
                plaq[d2][d1][X] = -plaq[d1][d2][X];
            }
        }
    }
}

/**
 * @brief compute the hodge-dual of a plaquette field
 * @tparam T Z-link group type
 * @param plaqin GaugeField[NDIM][NDIM] input plaquette field
 * @param plaqout GaugeField[NDIM][NDIM] output hodge-dual plaquette field
 */
template <typename T>
void dual_plaq(const sw_t<T> &plaqin, out_only sw_t<T> &plaqout) {
    if(NDIM==4) {
        foralldir(d) {
            onsites(ALL) plaqout[d][d][X] = 0;
        }
        int sign = -1;
        for (int i = 0; i < NDIM; ++i) {
            Direction d0 = Direction((0 + i) % NDIM);
            Direction d1 = Direction((1 + i) % NDIM);
            Direction d2 = Direction((2 + i) % NDIM);
            Direction d3 = Direction((3 + i) % NDIM);
            onsites(ALL) {
                plaqout[d0][d1][X] = sign * plaqin[d2][d3][X];
                plaqout[d0][d2][X] = -sign * plaqin[d1][d3][X];
                plaqout[d0][d3][X] = -sign * plaqin[d2][d1][X];
            }
            sign = -sign;
        }
    }
}

/**
 * @brief measure the absolue value of the monopole density per "time-direction" d0,
 *  as well as the monopole density per site-parity and "time-direction" d0.
 * @tparam T plaquette field data type
 * @param PlaquetteField plaq
 * @param monop_dens_per_d[NDIM] output absolute value of plaquette staple sum density per direction
 * @param monop_dens_per_p_d[2][NDIM] output plaquette staple sum density per parity and direction
 */
template <typename T>
void measure_stap_dens(const sw_t<T>& plaq, double(out_only &monop_dens_per_d)[NDIM],
                             double(out_only &monop_dens_per_p_d)[2][NDIM]) {
    if (NDIM == 4) {
        Field<double> stap;
        ReductionVector<double> mpdens_per_d(NDIM);
        mpdens_per_d = 0.0;
        mpdens_per_d.allreduce(false).delayed(true);
        ReductionVector<double> mpdens_per_p_d(NDIM * 2);
        mpdens_per_p_d = 0.0;
        mpdens_per_p_d.allreduce(false).delayed(true);
        foralldir(d1) {
            onsites(ALL) stap[X] = 0;
            foralldir(d2) if (d2 != d1) {
                onsites(ALL) {
                    stap[X] += (double)(plaq[d1][d2][X] - plaq[d1][d2][X - d2]);
                }
            }
            onsites(ALL) {
                mpdens_per_d[d1] += abs(stap[X]);
                int tpar = (int)uparity(X.coordinates()) - 1;
                mpdens_per_p_d[tpar * NDIM + d1] += stap[X];
            }
        }
        mpdens_per_d.reduce();
        mpdens_per_p_d.reduce();
        for (int i = 0; i < NDIM; ++i) {
            monop_dens_per_d[i] = mpdens_per_d[i] / lattice.volume();
            monop_dens_per_p_d[0][i] = mpdens_per_p_d[i] * 2.0 / lattice.volume();
            monop_dens_per_p_d[1][i] = mpdens_per_p_d[NDIM + i] * 2.0 / lattice.volume();
        }
    }
}

/**
 * @brief measure the absolue value of the monopole density per "time-direction" d0,
 *  as well as the monopole density per site-parity and "time-direction" d0.
 * @tparam T plaquette field data type
 * @param PlaquetteField plaq
 * @param monop_dens_per_d[NDIM] output absolute value of plaquette staple sum density per direction
 * @param monop_dens_per_p_d[2][NDIM] output plaquette staple sum density per parity and direction
 */
template <typename T>
void measure_staps_dens(const sw_t<T> &plaq, double(out_only &monop_dens_per_p_d)[2][NDIM]) {
    if (NDIM == 4) {
        Field<double> avplaq;
        Field<double> stap;
        ReductionVector<double> mpdens_per_p_d(NDIM * 2);
        mpdens_per_p_d = 0.0;
        mpdens_per_p_d.allreduce(false).delayed(true);
        foralldir(d1) {
            onsites(ALL) {
                avplaq[X] = 0;
                stap[X] = 0;
            }
            foralldir(d2) if (d2 != d1) {
                onsites(ALL) {
                    avplaq[X] += (double)(plaq[d1][d2][X] - plaq[d1][d2][X - d2]) / 6.0;
                }
            }
            foralldir(d2) if (d2 != d1) {
                onsites(ALL) {
                    stap[X] += (double)pow(plaq[d1][d2][X] - avplaq[X], 2.0) +
                               (double)pow(-plaq[d1][d2][X - d2] - avplaq[X], 2.0);
                }
            }
            onsites(ALL) {
                int tpar = (int)uparity(X.coordinates()) - 1;
                mpdens_per_p_d[tpar * NDIM + d1] += stap[X] / 6.0;
            }
        }
        mpdens_per_p_d.reduce();
        for (int i = 0; i < NDIM; ++i) {
            monop_dens_per_p_d[0][i] = mpdens_per_p_d[i] * 2.0 / lattice.volume();
            monop_dens_per_p_d[1][i] = mpdens_per_p_d[NDIM + i] * 2.0 / lattice.volume();
        }
    }
}


/**
 * @brief measure the absolue value of the monopole density per "time-direction" d0,
 *  as well as the monopole density per site-parity and "time-direction" d0.
 * @tparam T plaquette data type
 * @param PlaquetteField plaq
 * @param monop_dens_per_d[NDIM] output absolute value of monopole density per direction
 * @param monop_dens_per_p_d[2][NDIM] output monopole density per parity and direction
 */
template <typename T>
void measure_monop_dens(const sw_t<T> &plaq, double(out_only &monop_dens_per_d)[NDIM],
                           double(out_only &monop_dens_per_p_d)[2][NDIM]) {
    if(NDIM==4) {
        VectorField<double> M;
        ReductionVector<double> mpdens_per_d(NDIM);
        mpdens_per_d = 0.0;
        mpdens_per_d.allreduce(false).delayed(true);
        ReductionVector<double> mpdens_per_p_d(NDIM * 2);
        mpdens_per_p_d = 0.0;
        mpdens_per_p_d.allreduce(false).delayed(true);
        int sign = 1;
        for (int i = 0; i < NDIM; ++i) {
            Direction d0 = Direction((0 + i) % NDIM);
            Direction d1 = Direction((1 + i) % NDIM);
            Direction d2 = Direction((2 + i) % NDIM);
            Direction d3 = Direction((3 + i) % NDIM);
            onsites(ALL) {
                M[d0][X] = (double)(plaq[d1][d2][X + d3] - plaq[d1][d2][X]);
                M[d0][X] += (double)(plaq[d2][d3][X + d1] - plaq[d2][d3][X]);
                M[d0][X] += (double)(plaq[d3][d1][X + d2] - plaq[d3][d1][X]);
                M[d0][X] *= (double)sign;
                mpdens_per_d[d0] += abs(M[d0][X]);
                int tpar = (int)uparity(X.coordinates()) - 1;
                mpdens_per_p_d[tpar * NDIM + d0] += M[d0][X];
            }
            sign = -sign;
        }
        mpdens_per_d.reduce();
        mpdens_per_p_d.reduce();
        for (int i = 0; i < NDIM; ++i) {
            monop_dens_per_d[i] = mpdens_per_d[i] / lattice.volume();
            monop_dens_per_p_d[0][i] = mpdens_per_p_d[i] * 2.0 / lattice.volume();
            monop_dens_per_p_d[1][i] = mpdens_per_p_d[NDIM + i] * 2.0 / lattice.volume();
        }
    }
}

/**
 * @brief measure the average sign or sign density of H-loops winding around each of
 *  the NDIM periodic lattice directions.
 * @tparam T Z-link fied data type
 * @param hls_dens_per_p_d[2][NDIM] output H-loop sign density per parity and direction
 */
template <typename T>
void measure_nonplanar_curve(const GaugeField<T> &H, double(out_only &npl_dens_per_p_d)[2]) {
    
    ReductionVector<double> npldens_per_p_d(2);
    npldens_per_p_d = 0.0;
    npldens_per_p_d.allreduce(false).delayed(true);
    Field<T> Ff, Fb;
    PlaquetteField<T> Rf, Rb;
    foralldir(d1) foralldir(d2) if(d1 < d2) {
        onsites(ALL) {
            Rf[d1][d2][X] = H[d1][X] + H[d2][X + d1];
            Rf[d2][d1][X] = H[d2][X] + H[d1][X + d2];
            Rb[d1][d2][X] = H[d1][X - d1] + H[d2][X];
            Rb[d2][d1][X] = H[d2][X - d2] + H[d1][X];
        }
    }
    Direction d1 = e_x;
    Direction d2 = e_y;
    Direction d3 = e_z;
    Direction d4 = e_t;

    onsites(ALL) {
        Ff[X] = Rb[d1][d2][X] + Rf[d3][d4][X + d2]; //  {(x-d1),(x),(x+d2),(x+d2+d3),(x+d2+d3+d4)}
        Fb[X] = Rb[d4][d3][X] + Rf[d2][d1][X + d3]; //  {(x-d4),(x),(x+d3),(x+d2+d3),(x+d2+d3+d1)}
    }

    onsites(ALL) {
        int tpar = (int)uparity(X.coordinates()) - 1;
        npldens_per_p_d[tpar] += pow(Ff[X + d1] - Fb[X + d4], 2.0);
    }

    npldens_per_p_d.reduce();
    npl_dens_per_p_d[0] = npldens_per_p_d[0] * 2.0 / lattice.volume();
    npl_dens_per_p_d[1] = npldens_per_p_d[1] * 2.0 / lattice.volume();
}

/**
 * @brief measure the average sign or sign density of H-loops winding around each of
 *  the NDIM periodic lattice directions.
 * @tparam T Z-link fied data type
 * @param hls_dens_per_p_d[2][NDIM] output H-loop sign density per parity and direction
 */
template <typename T>
void measure_nonplanar_curve(const GaugeField<T> &H, double(out_only &npl_dens_per_p_d)[2][NDIM]) {

    ReductionVector<double> npldens_per_p_d(2 * NDIM);
    npldens_per_p_d = 0.0;
    npldens_per_p_d.allreduce(false).delayed(true);
    Field<T> Ff, Fb;
    PlaquetteField<T> Rf, Rb;
    foralldir(d1) foralldir(d2) if (d1 < d2) {
        onsites(ALL) {
            Rf[d1][d2][X] = H[d1][X] + H[d2][X + d1];
            Rf[d2][d1][X] = H[d2][X] + H[d1][X + d2];
            Rb[d1][d2][X] = H[d1][X - d1] + H[d2][X];
            Rb[d2][d1][X] = H[d2][X - d2] + H[d1][X];
        }
    }
    for (int i = 0; i < NDIM; ++i) {
        Direction d1 = Direction((i + 0) % NDIM);
        Direction d2 = Direction((i + 1) % NDIM);
        Direction d3 = Direction((i + 2) % NDIM);
        Direction d4 = Direction((i + 3) % NDIM);

        onsites(ALL) {
            Ff[X] = Rb[d1][d2][X] + H[d3][X + d2]; //  {(x-d1),(x),(x+d2),(x+d2+d3)}
            Fb[X] = H[d3][X - d3] + Rf[d2][d1][X]; //  {(x-d3),(x),(x+d2),(x+d2+d1)}
        }

        onsites(ALL) {
            int tpar = (int)uparity(X.coordinates()) - 1;
            npldens_per_p_d[tpar * NDIM + d4] += pow(Ff[X + d1] - Fb[X + d3], 2.0); // {(x),(x+d1),(x+d1+d2),(x+d1+d2+d3),(x+d2+d3),(x+d3),(x)}
        }
    }

    npldens_per_p_d.reduce();
    for (int i = 0; i < NDIM; ++i) {
        npl_dens_per_p_d[0][i] = npldens_per_p_d[0 * NDIM + i] * 2.0 / lattice.volume();
        npl_dens_per_p_d[1][i] = npldens_per_p_d[1 * NDIM + i] * 2.0 / lattice.volume();
    }
}

/**
 * @brief measure the average sign or sign density of H-loops winding around each of
 *  the NDIM periodic lattice directions.
 * @tparam T Z-link fied data type
 * @param hls_dens_per_p_d[2][NDIM] output H-loop sign density per parity and direction
 */
template <typename T>
void measure_hloop_sign_dens(const GaugeField<T> &H, double(out_only &hls_dens_per_p_d)[2][NDIM]) {
    Field<long> hloop;
    ReductionVector<double> hlsdens_per_p_d(NDIM * 2);
    hlsdens_per_p_d = 0.0;
    hlsdens_per_p_d.allreduce(false).delayed(true);

    foralldir(d) {
        onsites(ALL) hloop[X] = 0;
        // mult links so that hloop[X.dir == 0] contains the polyakov loop
        for (int i = lattice.size(d) - 2; i >= 0; --i) {

            // safe_access(polyakov) pragma allows the expression below, otherwise
            // hilapp would reject it because X and X+dir can refer to the same
            // site on different "iterations" of the loop.  However, here this
            // is restricted on single dir-plane so it works but we must tell it to hilapp.

    #pragma hila safe_access(hloop)
            onsites(ALL) {
                if (X.coordinate(d) == i) {
                    hloop[X] = H[d][X] + hloop[X + d];
                }
            }
        }

        onsites(ALL) if (X.coordinate(d) == 0) {
            int tpar = (int)uparity(X.coordinates()) - 1;
            hlsdens_per_p_d[tpar * NDIM + d] += (double)sgn(hloop[X]);
        }
    }

    hlsdens_per_p_d.reduce();
    for (int i = 0; i < NDIM; ++i) {
        hls_dens_per_p_d[0][i] = 
            hlsdens_per_p_d[i] * 2.0 * lattice.size(i) / lattice.volume();
        hls_dens_per_p_d[1][i] =
            hlsdens_per_p_d[NDIM + i] * 2.0 * lattice.size(i) / lattice.volume();
    }
}


template <typename T, typename fT>
void measure_stuff(const GaugeField<T> &H, const sw_t<fT> &sw, parameters& p) {
    // perform measurements on current link field H and plaquette shift field sw
    // and print results in formatted form to standard output
    static bool first = true;
    if (first) {
        // print legend for measurement output
        hila::out0 << "LSPLAQ    :         splaq        nsplaq\n";

        hila::out0 << "LSPLAQPPP :";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                    hila::out0 << "         p" << par << "d" << dir1 << dir2;
                }
            }
        }
        hila::out0 << "\n";

        if(0) {
            hila::out0 << "LNSPLAQPPP:";
            for (int par = 0; par < 2; ++par) {
                for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                    for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                        hila::out0 << "         p" << par << "d" << dir1 << dir2;
                    }
                }
            }
            hila::out0 << "\n";
        }

        hila::out0 << "LPLAQPPP  :";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                    hila::out0 << "         p" << par << "d" << dir1 << dir2;
                }
            }
        }
        hila::out0 << "\n";

        if(0) {
            hila::out0 << "LPLAQSPPP :";
            for (int par = 0; par < 2; ++par) {
                for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                    for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                        hila::out0 << "         p" << par << "d" << dir1 << dir2;
                    }
                }
            }
            hila::out0 << "\n";
        }

        if(0) {
            hila::out0 << "LNPLAQSPPP:";
            for (int par = 0; par < 2; ++par) {
                for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                    for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                        hila::out0 << "         p" << par << "d" << dir1 << dir2;
                    }
                }
            }
            hila::out0 << "\n";
        }

        if(0) {
            hila::out0 << "LHSQ      :           hsq\n";

            hila::out0 << "LHSQPPD   :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";

            hila::out0 << "LHPPD     :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }

        if (0) {
            hila::out0 << "LHLOOPSPPD:";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }
        if (0) {
            hila::out0 << "LNONPLCPPD:";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }
        if (1) {
            hila::out0 << "LSTAPPD   :";
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << "            d" << dir;
            }
            hila::out0 << "\n";

            hila::out0 << "LSTAPPPD  :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";

            hila::out0 << "LSTAPSPPD :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }
        if (0) {
            hila::out0 << "LMONPD    :";
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << "            d" << dir;
            }
            hila::out0 << "\n";

            hila::out0 << "LMONPPD   :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }
        if(0) {
            hila::out0 << "LDMONPD   :";
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << "            d" << dir;
            }
            hila::out0 << "\n";

            hila::out0 << "LDMONPPD  :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }

        if (p.n_ps_update > 0) {
            hila::out0 << "LSWPLAQPPP:";
            for (int par = 0; par < 2; ++par) {
                for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                    for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                        hila::out0 << "         p" << par << "d" << dir1 << dir2;
                    }
                }
            }
            hila::out0 << "\n";

            hila::out0 << "LSWMONPD  :";
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << "            d" << dir;
            }
            hila::out0 << "\n";

            hila::out0 << "LSWMONPPD :";
            for (int par = 0; par < 2; ++par) {
                for (int dir = 0; dir < NDIM; ++dir) {
                    hila::out0 << "          p" << par << "d" << dir;
                }
            }
            hila::out0 << "\n";
        }


        first = false;
    }
    auto splaq = measure_s_plaq(H, sw) / (lattice.volume() * NDIM * (NDIM - 1) / 2);
    auto nsplaq = measure_ns_plaq(H) / (lattice.volume() * NDIM * (NDIM - 1) / 2);

    hila::out0 << string_format("SPLAQ       % 0.6e % 0.6e\n", splaq, nsplaq);

    sw_t<T> plaq;
    plaq_field(H, plaq);

    double plaq_per_par_pl[2][NDIM][NDIM];
    measure_splaq_per_par_and_plane(plaq, sw, plaq_per_par_pl);
    hila::out0 << "SPLAQPPP   ";
    for (int par = 0; par < 2; ++par) {
        for (int dir1 = 0; dir1 < NDIM; ++dir1) {
            for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                hila::out0 << string_format(" % 0.6e", plaq_per_par_pl[par][dir1][dir2]);
            }
        }
    }
    hila::out0 << '\n';

    if(0) {
        hila::out0 << "NSPLAQPPP  ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                    hila::out0 << string_format(" % 0.6e", plaq_per_par_pl[par][dir2][dir1]);
                }
            }
        }
        hila::out0 << '\n';
    }

    sw_t<fT> totplaq;
    foralldir(d1) foralldir(d2) {
        onsites(ALL) totplaq[d1][d2][X] = (fT)plaq[d1][d2][X] + sw[d1][d2][X];
    }

    measure_plaq_per_par_and_plane(totplaq, plaq_per_par_pl);
    hila::out0 << "PLAQPPP    ";
    for (int par = 0; par < 2; ++par) {
        for (int dir1 = 0; dir1 < NDIM; ++dir1) {
            for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                hila::out0 << string_format(" % 0.6e", plaq_per_par_pl[par][dir1][dir2]);
            }
        }
    }
    hila::out0 << '\n';

    if(0) {
        measure_plaq_sign_per_par_and_plane(totplaq, plaq_per_par_pl);
        hila::out0 << "PLAQSPPP   ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                    hila::out0 << string_format(" % 0.6e", plaq_per_par_pl[par][dir1][dir2]);
                }
            }
        }
        hila::out0 << '\n';
    }

    if(0) {
        measure_plaq_sign_per_par_and_plane(plaq, plaq_per_par_pl);
        hila::out0 << "NPLAQSPPP  ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                    hila::out0 << string_format(" % 0.6e", plaq_per_par_pl[par][dir1][dir2]);
                }
            }
        }
        hila::out0 << '\n';
    }

    if (0) {
        double h_per_par_dir[2][NDIM];
        double hsq_per_par_dir[2][NDIM];
        auto hsq = measure_hsq(H, hsq_per_par_dir, h_per_par_dir);
        hila::out0 << string_format("HSQ         % 0.6e\n", hsq);

        hila::out0 << "HSQPPD     ";
        for (int par = 0; par < 2; ++par) {
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << string_format(" % 0.6e", hsq_per_par_dir[par][dir]);
            }
        }
        hila::out0 << '\n';

        hila::out0 << "HPPD       ";
        for (int par = 0; par < 2; ++par) {
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << string_format(" % 0.6e", h_per_par_dir[par][dir]);
            }
        }
        hila::out0 << '\n';
    }

    if (0) {
        double hl_per_par_dir[2][NDIM];
        measure_hloop_sign_dens(H, hl_per_par_dir);
        hila::out0 << "HLOOPSPPD  ";
        for (int par = 0; par < 2; ++par) {
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << string_format(" % 0.6e", hl_per_par_dir[par][dir]);
            }
        }
        hila::out0 << '\n';
    }

    if (0) {
        double nplc_per_par[2][NDIM];
        measure_nonplanar_curve(H, nplc_per_par);
        hila::out0 << "NONPLCPPD  ";
        for (int par = 0; par < 2; ++par) {
            for (int dir = 0; dir < NDIM; ++dir) {
                hila::out0 << string_format(" % 0.6e", nplc_per_par[par][dir]);
            }
        }
        hila::out0 << '\n';
    }

    if (1) {
        double m_per_dir[NDIM];
        double m_per_par_dir[2][NDIM];
        measure_stap_dens(totplaq, m_per_dir, m_per_par_dir);

        hila::out0 << "STAPPD     ";
        for (int dir1 = 0; dir1 < NDIM; ++dir1) {
            hila::out0 << string_format(" % 0.6e", m_per_dir[dir1]);
        }
        hila::out0 << '\n';

        hila::out0 << "STAPPPD    ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                hila::out0 << string_format(" % 0.6e", m_per_par_dir[par][dir1]);
            }
        }
        hila::out0 << '\n';

        measure_staps_dens(totplaq, m_per_par_dir);
        hila::out0 << "STAPSPPD   ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                hila::out0 << string_format(" % 0.6e", m_per_par_dir[par][dir1]);
            }
        }
        hila::out0 << '\n';
    }

    if (0) {
        double m_per_dir[NDIM];
        double m_per_par_dir[2][NDIM];
        measure_monop_dens(totplaq, m_per_dir, m_per_par_dir);

        hila::out0 << "MONPD      ";
        for (int dir1 = 0; dir1 < NDIM; ++dir1) {
            hila::out0 << string_format(" % 0.6e", m_per_dir[dir1]);
        }
        hila::out0 << '\n';

        hila::out0 << "MONPPD     ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                hila::out0 << string_format(" % 0.6e", m_per_par_dir[par][dir1]);
            }
        }
        hila::out0 << '\n';
    }

    if(0) {
        sw_t<fT> dualplaq;
        dual_plaq(totplaq, dualplaq);
        double m_per_dir[NDIM];
        double m_per_par_dir[2][NDIM];
        measure_monop_dens(dualplaq, m_per_dir, m_per_par_dir);

        hila::out0 << "DMONPD     ";
        for (int dir1 = 0; dir1 < NDIM; ++dir1) {
            hila::out0 << string_format(" % 0.6e", m_per_dir[dir1]);
        }
        hila::out0 << '\n';

        hila::out0 << "DMONPPD    ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                hila::out0 << string_format(" % 0.6e", m_per_par_dir[par][dir1]);
            }
        }
        hila::out0 << '\n';

    }

    if (p.n_ps_update > 0) {
        measure_plaq_per_par_and_plane(sw, plaq_per_par_pl);
        hila::out0 << "SWPLAQPPP  ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                for (int dir2 = dir1 + 1; dir2 < NDIM; ++dir2) {
                    hila::out0 << string_format(" % 0.6e", plaq_per_par_pl[par][dir1][dir2]);
                }
            }
        }
        hila::out0 << '\n';

        double m_per_dir[NDIM];
        double m_per_par_dir[2][NDIM];
        measure_monop_dens(sw, m_per_dir, m_per_par_dir);
        hila::out0 << "SWMONPD    ";
        for (int dir1 = 0; dir1 < NDIM; ++dir1) {
            hila::out0 << string_format(" % 0.6e", m_per_dir[dir1]);
        }
        hila::out0 << '\n';

        hila::out0 << "SWMONPPD   ";
        for (int par = 0; par < 2; ++par) {
            for (int dir1 = 0; dir1 < NDIM; ++dir1) {
                hila::out0 << string_format(" % 0.6e", m_per_par_dir[par][dir1]);
            }
        }
        hila::out0 << '\n';
    }
}

// end measurement functions
///////////////////////////////////////////////////////////////////////////////////
// load/save config functions

template <typename T, typename fT>
void checkpoint(const GaugeField<T> &U, const PlaquetteField<fT> &sw, int trajectory, const parameters &p) {
    double t = hila::gettime();
    // name of config with extra suffix
    std::string config_file =
        p.config_file + "_" + std::to_string(abs((trajectory + 1) / p.n_save) % 2);
    // save config
    U.config_write(config_file);
    sw.config_write(config_file + "_sw");
    // write run_status file
    if (hila::myrank() == 0) {
        std::ofstream outf;
        outf.open("run_status", std::ios::out | std::ios::trunc);
        outf << "trajectory  " << trajectory + 1 << '\n';
        outf << "seed        " << static_cast<uint64_t>(hila::random() * (1UL << 61)) << '\n';
        outf << "time        " << hila::gettime() << '\n';
        // write config name to status file:
        outf << "config name  " << config_file << '\n';
        outf.close();
    }
    std::stringstream msg;
    msg << "Checkpointing, time " << hila::gettime() - t;
    hila::timestamp(msg.str().c_str());
}

template <typename T, typename fT>
bool restore_checkpoint(GaugeField<T> &U, PlaquetteField<fT> &sw, int &trajectory, parameters &p) {
    uint64_t seed;
    bool ok = true;
    p.time_offset = 0;
    hila::input status;
    if (status.open("run_status", false, false)) {
        hila::out0 << "RESTORING FROM CHECKPOINT:\n";
        trajectory = status.get("trajectory");
        seed = status.get("seed");
        p.time_offset = status.get("time");
        // get config name with suffix from status file:
        std::string config_file = status.get("config name");
        status.close();
        hila::seed_random(seed);
        U.config_read(config_file);
        sw.config_read(config_file + "_sw");
        ok = true;
    } else {
        std::ifstream in;
        in.open(p.config_file, std::ios::in | std::ios::binary);
        if (in.is_open()) {
            in.close();
            hila::out0 << "READING initial config\n";
            U.config_read(p.config_file);
            std::ifstream in_sw;
            in_sw.open(p.config_file, std::ios::in | std::ios::binary);
            if (in_sw.is_open()) {
                in_sw.close();
                sw.config_read(p.config_file + "_sw");
                ok = true;
            } else {
                ok = false;
            }
        } else {
            ok = false;
        }
    }
    return ok;
}

// end load/save config functions
///////////////////////////////////////////////////////////////////////////////////


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

    hila::out0 << "Z-link theory with plaquette shift field\n";

    hila::out0 << "Using floating point epsilon: " << fp<ftype>::epsilon << "\n";

#if PARITY==0
    hila::out0 << "Using full 4D coordiantes for staggering\n";
#else
    hila::out0 << "Using spatial 3D coordiantes for staggering\n";
#endif

#if PLAQ_SHIFT == 1
    hila::out0 << "Using site-parity-dependent plaquette shifts on spatial plaquettes\n";
#elif PLAQ_SHIFT == 2
    hila::out0 << "Using site-parity-dependent plaquette shifts on temporal plaquettes\n";
#elif PLAQ_SHIFT == 3
    hila::out0 << "Using site-parity-dependent plaquette shifts on temporal and spatial plaquettes\n";
#else
    hila::out0 << "All plaquette shifts set to zero\n";
#endif

    parameters p;

    hila::input par("parameters");

    CoordinateVector lsize;
    // reads NDIM numbers
    lsize = par.get("lattice size");
    // gauge coupling
    p.beta = par.get("beta");
    // number of trajectories
    p.n_traj = par.get("number of trajectories");
    // number of Metropolis sweeps per trajectory
    p.n_update = par.get("metropolis updates");
    // number of overrelaxation sweeps per Metropolis sweep
    p.n_or_update = par.get("overrelax updates");
    // number of s_{x,\mu\nu} relaxation sweeps per Metropolis sweep
    p.n_ps_update = par.get("plaq shift updates");
    // number of thermalization trajectories
    p.n_therm = par.get("thermalization trajs");
    // random seed = 0 -> get seed from time
    long seed = par.get("random seed");
    // save config and checkpoint
    p.n_save = par.get("trajs/saved");
    // measure surface properties and print "profile"
    p.config_file = par.get("config name");

    par.close(); // file is closed also when par goes out of scope

    // set up the lattice
    lattice.setup(lsize);

    // We need random number here
    hila::seed_random(seed);

    // instantiate the gauge field
    GaugeField<mygroup> H;

 


    // define the plaquette shifts 
    sw_t<ftype> sw;
    foralldir(d1) foralldir(d2) {
        onsites(ALL) sw[d1][d2][X] = 0;
    }
#if PLAQ_SHIFT == 1 || PLAQ_SHIFT == 3
    for (int i = 0; i < NDIM - 1; ++i) {
        Direction d1 = Direction((1 + i) % (NDIM - 1));
        Direction d2 = Direction((2 + i) % (NDIM - 1));
        onsites(ALL) {
            if (uparity(X.coordinates()) == Parity::even) {
                sw[d1][d2][X] = 0.5;
            } else {
                sw[d1][d2][X] = -0.5;
            }
            sw[d2][d1][X] = -sw[d1][d2][X];
        }
    }
#endif
#if PLAQ_SHIFT == 2 || PLAQ_SHIFT == 3
    Direction d4 = Direction::e_t;
    for (int i = 0; i < NDIM - 1; ++i) {
        Direction d1 = Direction(i);
        onsites(ALL) {
            if (uparity(X.coordinates()) == Parity::even) {
                sw[d1][d4][X] = 0.5;
            } else {
                sw[d1][d4][X] = -0.5;
            }
            sw[d4][d1][X] = -sw[d1][d4][X];
        }
    }
#endif


    // use negative trajectory numbers for thermalisation
    int start_traj = -p.n_therm;

    if (!restore_checkpoint(H, sw, start_traj, p)) {
        foralldir(d) {
            onsites(ALL) {
                H[d][X] = 0;
            }
        }
    }


    hila::timer update_timer("Updates");
    hila::timer measure_timer("Measurements");

    ftype t_step0 = 0;
    for (int trajectory = start_traj; trajectory < p.n_traj; ++trajectory) {

        ftype ttime = hila::gettime();

        update_timer.start();

        do_trajectory(H, sw, p);

        // put sync here in order to get approx gpu timing
        hila::synchronize_threads();
        update_timer.stop();

        measure_timer.start();

        hila::out0 << "Measure_start " << trajectory << '\n';

        measure_stuff(H, sw, p);

        hila::out0 << "Measure_end " << trajectory << '\n';

        measure_timer.stop();

        if (p.n_save > 0 && (trajectory + 1) % p.n_save == 0) {
            checkpoint(H, sw, trajectory, p);
        }
    }

    hila::finishrun();
}
