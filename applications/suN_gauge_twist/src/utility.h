#ifndef UTILITY_H_
#define UTILITY_H_

#include "hila.h"
#include <algorithm>
int z_ind(int z) {
    return (z + lattice.size(e_z)) % lattice.size(e_z);
}

template <typename T>
void print_formatted_numbers(const std::vector<Complex<T>> &numbers, std::string label,
                             bool indices, bool sum = false) {

    if (indices) {
        hila::out0 << label << ": ";
        for (int i = 0; i < numbers.size() - 1; i++) {
            hila::out0 << i << ", ";
        }
        if (sum == true) {
            hila::out0 << "sum\n";
        } else {
            hila::out0 << numbers.size() - 1 << "\n";
        }
    } else {
        hila::out0 << label << ": ";
        for (int i = 0; i < numbers.size() - 1; i++) {
            hila::out0 << numbers[i].re << " " << numbers[i].im << ", ";
        }
        hila::out0 << numbers[numbers.size() - 1].re << " " << numbers[numbers.size() - 1].im
                   << "\n";
    }
}

template <typename T>
void write_surface(std::vector<T> surf) {
    std::ofstream MFile;
    MFile.open("surface", std::ios_base::app);
    hila::out0 << "OPENED surface file" << std::endl;
    for (int y = 0; y < lattice.size(e_y); y++) {
        for (int x = 0; x < lattice.size(e_x); x++) {
            if (hila::myrank() == 0)
                MFile << x << ' ' << y << ' ' << surf[x + y * lattice.size(e_x)] << '\n';
        }
    }
#endif