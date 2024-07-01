/** @file wilson_line_and_force.h */

#ifndef WILSON_LINE_AND_FORCE_H_
#define WILSON_LINE_AND_FORCE_H_

#include "hila.h"
#include <vector>

// functions to compute general Wilson lines and gauge force for closed Wilson lines

template <typename group>
void get_wilson_line(const GaugeField<group> &U, const std::vector<Direction> &path,
                     Field<group> &R) {
    // compute the Wilson line defined by the list of directions "path"
    int i, ip;
    int L = path.size();
    int udirs[NDIRS] = {0};
    Direction dir;
    for (i = 0; i < L; ++i) {
        dir = path[i];
        if (is_up_dir(dir)) {
            if ((udirs[(int)dir]++) == 0) {
                U[dir].start_gather(-dir, ALL);
            }
        }
    }

    Field<group> R0[2];

    i = 0;
    ip = 0;
    dir = path[i];
    // initialize R0[0] with first link variable of the Wilson line:
    if (is_up_dir(dir)) {
        // link points in positive direction
        onsites(ALL) R0[ip][X] = U[dir][X - dir];
    } else {
        // link points in negative direction
        onsites(ALL) R0[ip][X] = U[-dir][X].dagger();
    }

    // multiply R0[ip] successively with the L-2 intermediate link variables of the Wilson line
    // and store the result in R0[1-ip]
    for (i = 1; i < L - 1; ++i) {
        dir = path[i];
        R0[ip].start_gather(-dir, ALL);
        if (is_up_dir(dir)) {
            // link points in positive direction
            onsites(ALL) mult(R0[ip][X - dir], U[dir][X - dir], R0[1 - ip][X]);
        } else {
            // link points in negative direction
            onsites(ALL) mult(R0[ip][X - dir], U[-dir][X].dagger(), R0[1 - ip][X]);
        }
        ip = 1 - ip;
    }

    // multiply R0[ip] by the last link variable of the Wilson line and store the result in R
    dir = path[i];
    R0[ip].start_gather(-dir, ALL);
    if (is_up_dir(dir)) {
        // link points in positive direction
        onsites(ALL) mult(R0[ip][X - dir], U[dir][X - dir], R[X]);
    } else {
        // link points in negative direction
        onsites(ALL) mult(R0[ip][X - dir], U[-dir][X].dagger(), R[X]);
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_wloop_force_from_wl_add(const GaugeField<group> &U, const std::vector<Direction> &path,
                                 const Field<group> &W, atype eps, VectorField<Algebra<group>> &K) {
    // compute gauge force of Wilson loop "W", corresponding to "path" and add result to vector
    // field "K"
    Field<group> R = W;
    Field<group> R0;
    int L = path.size();
    int udirs[NDIRS] = {0};
    Direction dir;
    for (int i = 0; i < L; ++i) {
        dir = path[i];
        if (is_up_dir(dir)) {
            if ((udirs[(int)dir]++) == 0) {
                U[dir].start_gather(-dir, ALL);
            }
        }
    }

    for (int i = 0; i < L; ++i) {
        dir = path[i];
        R.start_gather(-dir, ALL);
        if (is_up_dir(dir)) {
            // link points in positive direction
            onsites(ALL) K[dir][X] -= R[X].project_to_algebra_scaled(eps);

            onsites(ALL) mult(U[dir][X - dir].dagger(), R[X - dir], R0[X]);
            onsites(ALL) mult(R0[X], U[dir][X - dir], R[X]);
        } else {
            // link points in negative direction
            onsites(ALL) mult(U[-dir][X], R[X - dir], R0[X]);
            onsites(ALL) mult(R0[X], U[-dir][X].dagger(), R[X]);

            onsites(ALL) K[-dir][X] += R[X].project_to_algebra_scaled(eps);
        }
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_wloop_force_from_wl(const GaugeField<group> &U, const std::vector<Direction> &path,
                             const Field<group> &W, atype eps, VectorField<Algebra<group>> &K) {
    // compute gauge force of Wilson loop "W", corresponding to "path" and store result in vector
    // field "K"
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_wloop_froce_add(U, path, W, eps, K);
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_wloop_force_add(const GaugeField<group> &U, const std::vector<Direction> &path, atype eps,
                         VectorField<Algebra<group>> &K) {
    // compute gauge force of Wilson loop along "path" and add result to vector field "K"
    Field<group> R, R0;
    int L = path.size();
    get_wilson_line(U, path, R);

    Direction dir;
    int udirs[NDIRS] = {0};
    for (int i = 0; i < L; ++i) {
        dir = path[i];
        if (is_up_dir(dir)) {
            if ((udirs[(int)dir]++) == 0) {
                U[dir].start_gather(-dir, ALL);
            }
        }
    }

    for (int i = 0; i < L; ++i) {
        dir = path[i];
        R.start_gather(-dir, ALL);
        if (is_up_dir(dir)) {
            // link points in positive direction
            onsites(ALL) K[dir][X] -= R[X].project_to_algebra_scaled(eps);

            onsites(ALL) mult(U[dir][X - dir].dagger(), R[X - dir], R0[X]);
            onsites(ALL) mult(R0[X], U[dir][X - dir], R[X]);
        } else {
            // link points in negative direction
            onsites(ALL) mult(U[-dir][X], R[X - dir], R0[X]);
            onsites(ALL) mult(R0[X], U[-dir][X].dagger(), R[X]);

            onsites(ALL) K[-dir][X] += R[X].project_to_algebra_scaled(eps);
        }
    }
}

template <typename group, typename atype = hila::arithmetic_type<group>>
void get_wloop_force(const GaugeField<group> &U, const std::vector<Direction> &path, atype eps,
                     VectorField<Algebra<group>> &K) {
    // compute gauge force of Wilson loop defined by "path" and store result in vector field "K"
    foralldir(d1) {
        K[d1][ALL] = 0;
    }
    get_wloop_froce_add(U, path, eps, K);
}


/*
* old, slow implementation:
template <typename group,int L>
void get_wilson_line_b(const GaugeField<group>& U,const Direction(&path)[L],Field<group>& R) {
    // compute the Wilson line defined by the list of directions "path"
    CoordinateVector v=0;
    // initialize R with first link of the Wilson line:
    if(path[0]<NDIM) {
        // link points in positive direction
        onsites(ALL) R[X]=U[path[0]][X];
        v+=path[0];
    } else {
        // link points in negative direction
        v+=path[0];
        onsites(ALL) R[X]=U[-path[0]][X+v].dagger();
    }

    // multiply R successively with the remaining links of the Wilson line
    for(int i=1; i<L; ++i) {
        if(path[i]<NDIM) {
            // link points in positive direction
            onsites(ALL) R[X]*=U[path[i]][X+v];
            v+=path[i];
        } else {
            // link points in negative direction
            v+=path[i];
            onsites(ALL) R[X]*=U[-path[i]][X+v].dagger();
        }
    }
}

template <typename group,int L>
void get_wloop_force_add_b(const GaugeField<group>& U,const Direction(&path)[L],ftype
eps,VectorField<Algebra<group>>& K) {
    // compute gauge force of Wilson loop defined by "path" and add result to vector field "K"
    Field<group> R;
    get_wilson_line_b(U,path,R);

    CoordinateVector v=0;
    for(int i=0; i<L; ++i) {
        if(path[i]<NDIM) {
            // link points in positive direction
            onsites(ALL) K[path[i]][X]-=R[X-v].project_to_algebra_scaled(eps);

            onsites(ALL) R[X]=U[path[i]][X+v].dagger()*R[X]*U[path[i]][X+v];

            v+=path[i];
        } else {
            // link points in negative direction
            v+=path[i];

            onsites(ALL) R[X]=U[-path[i]][X+v]*R[X]*U[-path[i]][X+v].dagger();

            onsites(ALL) K[-path[i]][X]+=R[X-v].project_to_algebra_scaled(eps);
        }
    }
}
*/


#endif