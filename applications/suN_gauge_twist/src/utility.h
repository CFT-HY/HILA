#ifndef UTILITY_H_
#define UTILITY_H_

#include "hila.h"
#include <algorithm>

enum class APPEND_FILE { TRUE, FALSE };
enum class CLOSE_FILE { TRUE, FALSE };

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

// template <typename T>
// void write_surface(std::vector<T> surf, std::string file_name, APPEND_FILE
// append_true,
//                    CLOSE_FILE close_file) {

//     static std::ofstream MFile; // Declare static ofstream object

//     // Open file if not already open, based on append_true parameter
//     if (!MFile.is_open()) {
//         if (append_true == APPEND_FILE::TRUE)
//             MFile.open(file_name, std::ios::app); // Append mode
//         else
//             MFile.open(file_name); // Default (overwrite) mode

//         if (!MFile.is_open()) {
//             std::cerr << "Error: Failed to open file " << file_name <<
//             std::endl; return;
//         }

//         std::cout << "Opened file " << file_name << std::endl;
//     }
//     // MFile.open("surface", std::ios_base::app);
//     MFile << "volume: " << lattice.size(e_x) << " " << lattice.size(e_y) << "
//     " << lattice.size(e_z)
//           << " " << lattice.size(e_t) << std::endl;
//     MFile << "x y z" << std::endl;
//     for (int y = 0; y < lattice.size(e_y); y++) {
//         for (int x = 0; x < lattice.size(e_x); x++) {
//             if (hila::myrank() == 0)
//                 MFile << x << ' ' << y << ' ' << surf[x + y *
//                 lattice.size(e_x)] << std::endl;
//         }
//     }
//     // Close file if specified by close_file parameter
//     if (close_file == CLOSE_FILE::TRUE) {
//         MFile.close();
//         std::cout << "Closed file " << file_name << std::endl;
//     }
// }
template <typename T>
void write_surface(std::vector<T> surf, std::string file_name, APPEND_FILE append_true,
                   CLOSE_FILE close_file) {

    std::ofstream MFile; // Declare non-static ofstream object

    // Open file based on append_true parameter
    if (append_true == APPEND_FILE::TRUE)
        MFile.open(file_name, std::ios::app); // Append mode
    else
        MFile.open(file_name); // Default (overwrite) mode

    if (!MFile.is_open()) {
        std::cerr << "Error: Failed to open file " << file_name << std::endl;
        return;
    }

    // Write lattice size and surface data
    MFile << "volume: " << lattice.size(e_x) << " " << lattice.size(e_y) << " " << lattice.size(e_z)
          << " " << lattice.size(e_t) << std::endl;
    MFile << "x y z" << std::endl;

    for (int y = 0; y < lattice.size(e_y); y++) {
        for (int x = 0; x < lattice.size(e_x); x++) {
            if (hila::myrank() == 0) {
                MFile << x << ' ' << y << ' ' << surf[x + y * lattice.size(e_x)] << std::endl;
            }
        }
    }

    // Close file if specified by close_file parameter
    if (close_file == CLOSE_FILE::TRUE) {
        MFile.close();
        std::cout << "Closed file " << file_name << std::endl;
    }
}

template <typename T>
void write_surface_binary(std::vector<T> surf, std::string file_name, APPEND_FILE append_true,
                          CLOSE_FILE close_file) {

    std::ofstream MFile; // Declare ofstream object

    // Open file based on append_true parameter
    if (append_true == APPEND_FILE::TRUE)
        MFile.open(file_name, std::ios::app | std::ios::binary); // Append mode
    else
        MFile.open(file_name, std::ios::binary); // Default (overwrite) mode

    if (!MFile.is_open()) {
        std::cerr << "Error: Failed to open file " << file_name << std::endl;
        return;
    }

    // Write lattice size
    int lattice_sizes[] = {lattice.size(e_x), lattice.size(e_y), lattice.size(e_z),
                           lattice.size(e_t)};
    MFile.write(reinterpret_cast<char *>(lattice_sizes), sizeof(lattice_sizes));

    // Write surface data
    for (int y = 0; y < lattice.size(e_y); y++) {
        for (int x = 0; x < lattice.size(e_x); x++) {
            if (hila::myrank() == 0) {
                int coords[] = {x, y};
                T value = surf[x + y * lattice.size(e_x)];
                MFile.write(reinterpret_cast<char *>(coords), sizeof(coords));
                MFile.write(reinterpret_cast<char *>(&value), sizeof(value));
            }
        }
    }

    // Close file if specified by close_file parameter
    if (close_file == CLOSE_FILE::TRUE) {
        MFile.close();
        std::cout << "Closed file " << file_name << std::endl;
    }
}

// template <typename T>
// void write_fourier(std::vector<T> npow, std::vector<int> hits, int pow_size,
// std::string file_name,
//                    APPEND_FILE append_true, CLOSE_FILE close_file) {

//     static std::ofstream MFile; // Declare static ofstream object

//     // Open file if not already open, based on append_true parameter

//     if (!MFile.is_open()) {
//         if (append_true == APPEND_FILE::TRUE)
//             MFile.open(file_name, std::ios::app); // Append mode
//         else
//             MFile.open(file_name); // Default (overwrite) mode

//         if (!MFile.is_open()) {
//             std::cerr << "Error: Failed to open file " << file_name <<
//             std::endl; return;
//         }

//         std::cout << "Opened file " << file_name << std::endl;
//     }
//     // MFile.open("surface", std::ios_base::app);
//     MFile << "volume: " << lattice.size(e_x) << " " << lattice.size(e_y) << "
//     " << lattice.size(e_z)
//           << " " << lattice.size(e_t) << std::endl;
//     MFile << "i"
//           << " n/h"
//           << " h" << std::endl;
//     for (int i = 0; i < pow_size; i++) {
//         if (hits[i] > 0)
//             MFile << i << ' ' << npow[i] / hits[i] << ' ' << hits[i] << '\n';
//     }
//     // Close file if specified by close_file parameter
//     if (close_file == CLOSE_FILE::TRUE) {
//         MFile.close();
//         std::cout << "Closed file " << file_name << std::endl;
//     }
// }

template <typename T>
void write_fourier(std::vector<T> npow, std::vector<int> hits, int pow_size, std::string file_name,
                   APPEND_FILE append_true, CLOSE_FILE close_file) {

    std::ofstream MFile; // Declare ofstream object

    // Open file based on append_true parameter
    if (append_true == APPEND_FILE::TRUE)
        MFile.open(file_name, std::ios::app); // Append mode
    else
        MFile.open(file_name); // Default (overwrite) mode

    if (!MFile.is_open()) {
        std::cerr << "Error: Failed to open file " << file_name << std::endl;
        return;
    }

    // Write lattice size and data
    MFile << "volume: " << lattice.size(e_x) << " " << lattice.size(e_y) << " " << lattice.size(e_z)
          << " " << lattice.size(e_t) << std::endl;
    MFile << "i"
          << " n/h"
          << " h" << std::endl;
    for (int i = 0; i < pow_size; i++) {
        if (hits[i] > 0)
            MFile << i << ' ' << npow[i] / hits[i] << ' ' << hits[i] << '\n';
    }

    // Close file if specified by close_file parameter
    if (close_file == CLOSE_FILE::TRUE) {
        MFile.close();
        std::cout << "Closed file " << file_name << std::endl;
    }
}

template <typename T>
void write_fourier_binary(std::vector<T> npow, std::vector<int> hits, int pow_size,
                          std::string file_name, APPEND_FILE append_true, CLOSE_FILE close_file) {

    std::ofstream MFile; // Declare ofstream object

    // Open file based on append_true parameter
    if (append_true == APPEND_FILE::TRUE)
        MFile.open(file_name, std::ios::app | std::ios::binary); // Append mode
    else
        MFile.open(file_name, std::ios::binary); // Default (overwrite) mode

    if (!MFile.is_open()) {
        std::cerr << "Error: Failed to open file " << file_name << std::endl;
        return;
    }

    // Write lattice size and data
    int lattice_sizes[] = {lattice.size(e_x), lattice.size(e_y), lattice.size(e_z),
                           lattice.size(e_t)};
    MFile.write(reinterpret_cast<char *>(lattice_sizes), sizeof(lattice_sizes));

    for (int i = 0; i < pow_size; i++) {
        if (hits[i] > 0) {
            MFile.write(reinterpret_cast<char *>(&i), sizeof(i));
            T value = npow[i] / hits[i];
            MFile.write(reinterpret_cast<char *>(&value), sizeof(value));
            MFile.write(reinterpret_cast<char *>(&hits[i]), sizeof(hits[i]));
        }
    }

    // Close file if specified by close_file parameter
    if (close_file == CLOSE_FILE::TRUE) {
        MFile.close();
        std::cout << "Closed file " << file_name << std::endl;
    }
}


#endif