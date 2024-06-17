#ifndef POLYAKOV_SURFACE_H_
#define POLYAKOV_SURFACE_H_

#include "hila.h"
#include "utility.h"
/// Helper function to get valid z-index coordinate

void spectraldensity_surface(std::vector<float> &surf, std::vector<double> &npow,
                             std::vector<int> &hits) {

    // do fft for the surface
    static bool first = true;
    static Complex<double> *buf;
    static fftw_plan fftwplan;

    int area = lattice.size(e_x) * lattice.size(e_y);

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

    int pow_size = npow.size();

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
}

template <typename T>
void measure_polyakov_field(const Field<T> &Ut, Field<Complex<float>> &polyakov_field) {
    Field<T> polyakov = Ut;

    // mult links so that polyakov[X.dir == 0] contains the polyakov loop
    for (int plane = lattice.size(e_t) - 2; plane >= 0; plane--) {

        // safe_access(polyakov) pragma allows the expression below, otherwise
        // hilapp would reject it because X and X+dir can refer to the same
        // site on different "iterations" of the loop.  However, here this
        // is restricted on single dir-plane so it works but we must tell it to hilapp.

#pragma hila safe_access(polyakov)
        onsites(ALL) {
            if (X.coordinate(e_t) == plane) {
                polyakov[X] = Ut[X] * polyakov[X + e_t];
            }
        }
    }

    onsites(ALL) if (X.coordinate(e_t) == 0) {
        polyakov_field[X] = trace(polyakov[X]);
    }
}

////////////////////////////////////////////////////////////////////
template <typename T>
void smear_polyakov_field(Field<T> &polyakov_field, int nsmear, float smear_coeff, double twist) {
    if (nsmear > 0) {
        Field<T> pl2 = 0;

        for (int i = 0; i < nsmear; i++) {
            onsites(ALL) if (X.coordinate(e_t) == 0) {
                
                pl2[X] = polyakov_field[X] +
                         smear_coeff * (polyakov_field[X + e_x] + polyakov_field[X - e_x] +
                                        polyakov_field[X + e_y] + polyakov_field[X - e_y]);
                if (X.coordinate(e_z) == 0)
                    pl2[X] = polyakov_field[X] +
                        smear_coeff * (expi(-2 * M_PI * twist/NCOLOR)*polyakov_field[X + e_z] + polyakov_field[X - e_z]);
                else if (X.coordinate(e_z) == 1)
                    pl2[X] = polyakov_field[X] +
                        smear_coeff *(polyakov_field[X + e_z] + expi(2 * M_PI * twist/NCOLOR)*polyakov_field[X - e_z]);
                else
                    pl2[X] = polyakov_field[X] +
                        smear_coeff *(polyakov_field[X + e_z] + polyakov_field[X - e_z]);
                
            }

            onsites(ALL) if (X.coordinate(e_t) == 0) {
                polyakov_field[X] = pl2[X] / (1 + 6 * smear_coeff);
            }
        }
    }
}


/**
 * @brief Polyakov average and value over z-index
 *
 * @tparam T
 * @param polyakov_field
 * @param surface_origin_profile z-index of Polyakov field at \f$x=0\f$ and \f$y=0\f$
 * @param surface_average_profile z-index Polyakov field sum
 */
template <typename T>
void measure_polyakov_profile(Field<T> &polyakov_field, std::vector<T> &surface_origin_profile,
                              std::vector<T> &surface_average_profile) {
    ReductionVector<T> p_surface_average(lattice.size(e_z)), p_origin(lattice.size(e_z));
    p_surface_average.allreduce(false);
    p_origin.allreduce(false);
    onsites(ALL) if (X.coordinate(e_t) == 0) {
        p_surface_average[X.z()] += polyakov_field[X];
        if (X.x() == 0 && X.y() == 0)
            p_origin[X.z()] += polyakov_field[X];
    }
    double inverse_surface_are = 1.0 / (lattice.size(e_x) * lattice.size(e_y));
    surface_origin_profile = p_origin.vector();
    surface_average_profile= p_surface_average.vector();
    std::transform(surface_average_profile.begin(), surface_average_profile.end(), surface_average_profile.begin(),
                [inverse_surface_are](T value) { return value * inverse_surface_are; });
}

/**
 * @brief Measrue surface tension through Polyakov loop fluctuations
 *
 * @tparam group gauge field group
 * @param U gauge field
 * @param p parameter list
 * @param traj current update step
 */
template <typename group>
void measure_polyakov_surface(GaugeField<group> &U, const parameters &p, int traj) {
    // Keep track of printing label
    static bool first = true;
    using MyType = Complex<float>;
    Field<MyType> polyakov_field;
    if (first) {
        std::vector<MyType> dummy(lattice.size(e_z));
        for (int sl = 0; sl < p.n_smear.size(); sl++) {
            print_formatted_numbers(dummy, "PRO " + std::to_string(sl) +" Origin polyakov smeared", true);
            print_formatted_numbers(dummy,"PRO " + std::to_string(sl) + " Average polyakov smeared", true);
            first = false;
        }
    }
    // 1. First we measure Polyakov field
    if (0) {
        // this section does local sums of poly lines
        Field<group> staples, Ut;

        Ut = U[e_t];

        for (int s = 0; s < 20; s++) {
            staplesum(U, staples, e_t);
            U[e_t][ALL] = (U[e_t][X] + 0.5 * staples[X]) / (1 + 6 * 0.5);
        }

        measure_polyakov_field(U[e_t], polyakov_field);

        U[e_t] = Ut;

    } else {
        // here standard non-link integrated field
        measure_polyakov_field(U[e_t], polyakov_field);
    }

    hila::out0 << std::setprecision(5);

    std::vector<MyType> surface_origin_profile, surface_average_profile;

    // 2. Smearing subprocess
    int prev_smear = 0;
    for (int sl = 0; sl < p.n_smear.size(); sl++) {

         // 2.1 Smear polyakov field
        int smear = p.n_smear.at(sl);
        smear_polyakov_field(polyakov_field, smear - prev_smear, p.smear_coeff, p.twist_coeff);
        prev_smear = smear;

        Field<MyType> polyakov_field_z = polyakov_field;
        if (p.z_smear.at(sl) > 0) {
            Field<MyType> sub_polyakov_field;
            for (int j = 0; j < p.z_smear.at(sl); j++) {
                double twist = p.twist_coeff;
                twist /= NCOLOR;
                hila::out0 << twist << '\n';
                onsites(ALL) if (X.coordinate(e_t) == 0) {
 
                    if (X.z() == 0)
                        sub_polyakov_field[X] =
                            polyakov_field_z[X] + p.smear_coeff * (expi(-2 * M_PI * twist)*polyakov_field_z[X + e_z] + polyakov_field_z[X - e_z]);
                    else if (X.z() == 1)
                        sub_polyakov_field[X] =
                            polyakov_field_z[X] + p.smear_coeff * (polyakov_field_z[X + e_z] + expi(2 * M_PI * twist)*polyakov_field_z[X - e_z]);
                    else
                        sub_polyakov_field[X] =
                            polyakov_field_z[X] + p.smear_coeff * (polyakov_field_z[X + e_z] + polyakov_field_z[X - e_z]);
                }
                onsites(ALL) if (X.coordinate(e_t) == 0) {
                    polyakov_field_z[X] = sub_polyakov_field[X] / (1 + 2 * p.smear_coeff);
                }
            }
        }

        measure_polyakov_profile(polyakov_field_z, surface_origin_profile, surface_average_profile);

        // double inverse_surface_are = 1.0 / (lattice.size(e_x) * lattice.size(e_y));
        // for (int i = 0; i < surface_average_profile.size(); i++) {
        //     surface_average_profile[i] *= inverse_surface_are;
        //     hila::out0 << "PRO " << sl << " " << i << " (" << surface_average_profile[i].re << ","
        //                << surface_average_profile[i].im << ") (" << surface_origin_profile[i].re
        //                << "," << surface_origin_profile[i].im << ")\n";
        // }
        // for (int i = 0; i < surface_average_profile.size(); i++) {
        //     surface_average_profile[i] *= inverse_surface_are;
        //     hila::out0 << "Polyakov_smeared: " << surface_average_profile[i].re << " " << surface_average_profile[i].im << ", ";
        // }
        hila::out0 << "PRO " << sl << "\n";
        
        print_formatted_numbers(surface_origin_profile,"PRO " + std::to_string(sl) + " Origin polyakov smeared", false);
        print_formatted_numbers(surface_average_profile,"PRO " + std::to_string(sl) + " Average polyakov smeared", false);
        float min = 1e8, max = 0;
        int minloc, maxloc;
        for (int i = 0; i < surface_average_profile.size(); i++) {
            if (min > surface_average_profile[i].abs()) {
                min = surface_average_profile[i].abs();
                minloc = i;
            }
            if (max < surface_average_profile[i].abs()) {
                max = surface_average_profile[i].abs();
                maxloc = i;
            }
        }
        hila::out0 << "PRO " << sl << " Min: " << surface_average_profile[minloc].abs() << " " << minloc << 
                   "\nPRO " << sl << " Max: " << surface_average_profile[maxloc].abs()<< " " << maxloc << std::endl;
        hila::synchronize_threads();
        //find the surface between minloc and maxloc
        float surface_level = max * 0.5; // assume min is really 0
        int area = lattice.size(e_x) * lattice.size(e_y);

        hila::out0 << "Surface_level" << sl << ' ' << surface_level << '\n';

        int startloc;
        if (maxloc > minloc)
            startloc = (maxloc + minloc) / 2;
        else
            startloc = ((maxloc + minloc + lattice.size(e_z)) / 2) % lattice.size(e_z);

        //starting positio for the other surface
        //startloc2 = z_ind(startloc + lattice.size(e_z) / 2);

        hila::out0 << "Start location: " << startloc << std::endl;

        std::vector<float> surf_interpolated, surf2;
        // Only allocate on first rank
        if (hila::myrank() == 0) {
            surf_interpolated.resize(area);
            //surf2.resize(area);
        }

        hila::out0 << std::setprecision(6);

        std::vector<MyType> polyakov_3D_volume;
        std::vector<MyType> line(lattice.size(e_z));

        //get full xyz-volume t=0 slice to main node
        polyakov_3D_volume = polyakov_field_z.get_slice({-1, -1, -1, 0});

        for (int y = 0; y < lattice.size(e_y); y++) {
            //get now full xz-plane polyakov line to main node
            //reduces MPI calls compared with doing line-by-line
            polyakov_3D_volume = polyakov_field_z.get_slice({-1, y, -1, 0});
            if (hila::myrank() == 0) {
                for (int x = 0; x < lattice.size(e_x); x++) {
                    line = polyakov_field_z.get_slice({x, y, -1, 0});

                    //copy ploop data to line - x runs fastest
                    for (int z = 0; z < lattice.size(e_z); z++) {
                        line[z] = polyakov_3D_volume[x + lattice.size(e_x) * (z)];
                    }

                    if (hila::myrank() == 0) {
                    //start search of the surface from the center between min and max
                    int z = startloc;

                    while (line[z_ind(z)].abs() > surface_level && startloc - z < lattice.size(e_z) *
                    0.4)
                        z--;

                    while (line[z_ind(z + 1)].abs() <= surface_level &&
                            z - startloc < lattice.size(e_z) * 0.4)
                        z++;


                    //do linear interpolation
                    //surf[x + y * lattice.size(e_x)] = z;
                    // surf_interpolated[x + y * lattice.size(e_x)] =
                    //     z +
                    //     (surface_level - line[z_ind(z)]) / (line[z_ind(z + 1)] - line[z_ind(z)]);

                    // if (p.n_surface > 0 && (traj + 1) % p.n_surface == 0) {
                    //     hila::out0 << "SURF" << sl << ' ' << x << ' ' << y << ' '
                    //                << surf_interpolated[x + y * lattice.size(e_x)] << '\n';
                    // }

                    //and locate the other surface - start from Lz/2 offset

                    // z = startloc2;

                    // while (line[z_ind(z)] <= surface_level &&
                    //        startloc2 - z < lattice.size(e_z) * 0.4)
                    //     z--;

                    // while (line[z_ind(z + 1)] > surface_level &&
                    //        z - startloc2 < lattice.size(e_z) * 0.4)
                    //     z++;

                    // //do linear interpolation
                    // surf[x + y * lattice.size(e_x)] = z;
                    // surf2[x + y * lattice.size(e_x)] =
                    //     z +
                    //     (surface_level - line[z_ind(z)]) / (line[z_ind(z + 1)] - line[z_ind(z)]);
                    }
                }
            }

            // if (hila::myrank() == 0) {
            //     constexpr int pow_size = 200;
            //     std::vector<double> npow(pow_size);
            //     std::vector<int> hits(pow_size);

            //     spectraldensity_surface(surf_interpolated, npow, hits);
            //     //spectraldensity_surface(surf2, npow, hits);

            //     for (int i = 0; i < pow_size; i++) {
            //         if (hits[i] > 0)
            //             hila::out0 << "POW" << sl << ' ' << i << ' ' << npow[i] / hits[i] << ' '
            //                     << hits[i] << '\n';
            //     }
            // }
        }
    }
}

#endif