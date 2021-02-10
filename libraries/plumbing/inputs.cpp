#include <sstream>
#include <iostream>
#include <fstream>
#include <regex>
#include <type_traits>
#include "inputs.h"
#include "globals.h"

#include "com_mpi.h"

#ifdef USE_MPI
static int myrank = 0;

#endif

/// handle one line of input in parameter file
void input::handle(const std::string &line) {
    std::regex pattern("\\s*([a-zA-Z_-]+[0-9]*)\\s*=\\s*([^\\s]+)\\s*");
    std::smatch results;
    if (!std::regex_match(line, results, pattern)) {
        return;
    };
    std::string variable(results[1]);
    std::string value(results[2]);
    bool is_numeric =
        (!value.empty() && value.find_first_not_of("0123456789.-") == std::string::npos);
    if (essentials.count(variable) != 0)
        essentials[variable] = true;
    if (is_numeric) {
        values[variable] = std::stod(value);
        hila::output << "found " + variable + " = " << values[variable] << "\n";
    } else {
        names[variable] = value;
        hila::output << "found " + variable + " = " << names[variable] << "\n";
    }
    if (essentials.count(variable) == 1) {
        essentials[variable] = true;
    }
}

void input::import(const std::string &fname) {

#ifdef USE_MPI
    int dummy = 0;
    char **argvp;

    MPI_Comm_rank(lattice->mpi_comm_lat, &myrank);
    if (myrank == 0) {
        read(fname);
        check_essentials();
    }
    broadcast_values();
    broadcast_names();

#else

    read(fname);
    check_essentials();

#endif
}

void input::import(int &argc, char ***argvp, const std::string &fname) {

#ifdef USE_MPI

    MPI_Comm_rank(lattice->mpi_comm_lat, &myrank);
    if (myrank == 0) {
        read(fname);
        check_essentials();
    }

    broadcast_values();
    broadcast_names();

#else

    read(fname);
    check_essentials();

#endif
}

void input::read(const std::string &fname) {
    std::ifstream inputfile;
    inputfile.open(fname);
    if (inputfile.is_open()) {
        std::string line;
        getline(inputfile, line, '\n');
        while (!inputfile.eof()) {
            handle(line);
            getline(inputfile, line, '\n');
        }
    } else {
        hila::output << "input file couldn't be opened. Checking default params...\n";
    }
    inputfile.close();
}

void input::add_essential(const std::string &var) {
    essentials.insert(std::pair<std::string, bool>(var, false));
}

template <typename T>
void input::add_essential(std::string const &var, T const &default_value) {
    hila::output << "type of default param not recognized for: " + var + "\n";
    essentials[var] = false;
}

template <> void input::add_essential(std::string const &var, int const &default_value) {
    values[var] = default_value;
    essentials[var] = true;
}
template <>
void input::add_essential(std::string const &var, float const &default_value) {
    values[var] = default_value;
    essentials[var] = true;
}
template <>
void input::add_essential(std::string const &var, double const &default_value) {
    values[var] = default_value;
    essentials[var] = true;
}
template <>
void input::add_essential(std::string const &var, std::string const &default_value) {
    names[var] = default_value;
    essentials[var] = true;
}

void input::check_essentials() {
    bool fail = false;
    for (auto i = essentials.begin(); i != essentials.end(); ++i) {
        if (!(*i).second) {
            hila::output << "required parameter " + (*i).first + " not found\n";
            fail = true;
        }
    }
    if (fail) {
        hila::output << "exiting...\n";
#ifdef USE_MPI
        hila::finishrun();
#else
        exit(1);
#endif
    }
}

input::returntype input::get(const std::string &variable) {
    return returntype(variable, this);
}

void input::close() { this->~input(); }

#ifdef USE_MPI
void input::broadcast_values() {
    double *vals; // buffer containing values for each name
    char *names;  // buffer containing variable names
    int lengths[3];
    int counter = 0;

    if (myrank == 0) {
        lengths[0] = values.size();
        lengths[1] = 0;
        lengths[2] = 0;
        for (auto i = values.begin(); i != values.end(); ++i) {
            lengths[1] += (int)(*i).first.size();
            lengths[2] += 1;
        }
    }

    // broadcast lengths to other processes
    MPI_Bcast(&lengths, 3, MPI_INT, 0, lattice->mpi_comm_lat);

    vals = new double[lengths[0]];
    names = new char[lengths[1] + lengths[0]];
    counter = lengths[2];

    // construct name and value lists in root node
    if (myrank == 0) {
        std::string buffer;
        int index = 0;
        for (auto i = values.begin(); i != values.end(); ++i) {
            vals[index] = (*i).second;
            buffer.append((*i).first + "\t");
            index++;
        }
        snprintf(names, lengths[1] + lengths[0], "%s", buffer.c_str());
    }

    MPI_Bcast(vals, lengths[2], MPI_DOUBLE, 0, lattice->mpi_comm_lat);
    MPI_Bcast(names, lengths[1] + lengths[0], MPI_CHAR, 0, lattice->mpi_comm_lat);

    // construct map in other nodes
    if (myrank != 0) {
        std::istringstream iss(std::string(names), std::istringstream::in);
        for (int j = 0; j < counter; j++) {
            std::string temp;
            iss >> temp;
            values[temp] = vals[j];
        }
    }

    delete[] vals;
    delete[] names;
}

void input::broadcast_names() {
    char *vars;    // buffer containing variables
    char *strings; // buffer containing the strings for each variable
    int lengths[3];

    if (myrank == 0) {
        lengths[0] = names.size();
        lengths[1] = 0;
        lengths[2] = 0;
        for (auto i = names.begin(); i != names.end(); ++i) {
            lengths[1] += (int)(*i).first.size();
            lengths[2] += (int)(*i).second.size();
        }
    }

    // broadcast lengths to other processes
    MPI_Bcast(&lengths, 3, MPI_INT, 0, lattice->mpi_comm_lat);

    vars = new char[lengths[1] + lengths[0]];
    strings = new char[lengths[2] + lengths[0]];

    int counter = 0;
    if (myrank == 0) {
        std::string buffer1;
        std::string buffer2;
        for (auto i = names.begin(); i != names.end(); ++i) {
            buffer1.append((*i).first + "\t");
            buffer2.append((*i).second + "\t");
            counter++;
        }
        snprintf(vars, lengths[1] + lengths[0], "%s", buffer1.c_str());
        snprintf(strings, lengths[2] + lengths[0], "%s", buffer1.c_str());
    }

    MPI_Bcast(vars, lengths[0], MPI_DOUBLE, 0, lattice->mpi_comm_lat);
    MPI_Bcast(strings, lengths[1] + lengths[0], MPI_CHAR, 0, lattice->mpi_comm_lat);

    // construct name-string map in other nodes
    if (myrank != 0) {
        std::istringstream iss1(std::string(vars), std::istringstream::in);
        std::istringstream iss2(std::string(strings), std::istringstream::in);
        for (int j = 0; j < counter; j++) {
            std::string temp1;
            std::string temp2;
            iss1 >> temp1;
            iss2 >> temp2;
            names[temp1] = temp2;
        }
    }

    delete[] strings;
    delete[] vars;
}
#endif
