#ifndef PARAMETERS_H
#define PARAMETERS_H

using mygroup = SU<NCOLOR, double>;

enum class poly_limit { OFF, RANGE, PARABOLIC };


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
    int twist_coeff;
};

#endif