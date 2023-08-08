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
    poly_limit polyakov_pot;
    double poly_min, poly_max, poly_m2;
    std::vector<int> n_smear;
    double smear_coeff;
    std::vector<int> z_smear;
    int n_surface;
    int n_dump_polyakov;
};

#endif