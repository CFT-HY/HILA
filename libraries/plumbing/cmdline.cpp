#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <regex>
#include <charconv>
#include "hila.h"

using strvec = std::vector<std::string>;

////////////////////////////////////////////////////////////////////////////////
/// @brief Prints the help texts associated with recognized flags and quits the
///        program
////////////////////////////////////////////////////////////////////////////////
void cmdlinearguments::quit_with_help() {
    print_help();
    hila::finishrun();
}

////////////////////////////////////////////////////////////////////////////////
/// @struct argmap_val
/// @brief Struct to hold structured information on used command line arguments
///
/// @var strvec argmap_val::val_strings
/// @brief Strings corresponding to read in cmdline parameters
///
/// @var std::string argmap_val::help_text
/// @brief A string describing the use of the flag
///
/// @var int argmap_val::number
/// @brief number of args for each flag.  < 0: unspecified
///
/// @var bool argmap_val::present
/// @brief A boolean telling whether the flag was found in argv
////////////////////////////////////////////////////////////////////////////////
struct argmap_val {
    strvec val_strings;
    std::string aux;
    std::string help_text;
    int number;
    bool present;
};

// class instance cmdline is globally available through hila.h
namespace hila {
cmdlinearguments cmdline;
}

cmdlinearguments::cmdlinearguments() {}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copies provided argc and argv into internal variables.
///
/// @param argc   length of argv
/// @param argv   array of pointers to command-line arguments
////////////////////////////////////////////////////////////////////////////////
void cmdlinearguments::initialise_args(int argc0, char **argv0) {
    argc = argc0;
    argv = (const char **)malloc(argc * sizeof(const char *));
    for (int i = 0; i < argc; i++)
        argv[i] = argv0[i];

    // Proceed to parse argv
    fill_argmap();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Returns the vector of strings associated with flag 'flag'.
///
/// @param flag   flag given as string
/// @return vector of arguments
////////////////////////////////////////////////////////////////////////////////
strvec cmdlinearguments::values(const std::string &flag) {
    strvec valvec = argmap[flag].val_strings;
    if (argmap.count(flag) == 0) {
        hila::out0 << "Flag '" << flag << "' is not recognized!\n";
        quit_with_help();
    }

    if (valvec.size() == 0) {
        hila::out0 << "\n\nFlag '" << flag << "' has no entries"
                   << " and the associated vector of strings is of length zero!\n\n";
    }
    return valvec;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Adds a new recognized flag 'flag' and associated information to the
///        argmap.
///
/// @param flag        string with the name of the flag, always of format "^-[a-z].*"
///                    (no whitespace allowed), e.g. '-my_flag123'
/// @param help_text   a string containing useful knowledge of the usage of the flag.
///                    Will be printed when error occurs. For formatting, a newline
///                    can be used to split the help_text, as it will help with
///                    correct formatting in print_help().
/// @param number      required number of args for each key appearance.  default -1: number
///                    unspecified
////////////////////////////////////////////////////////////////////////////////
void cmdlinearguments::add_flag(const std::string &flag, const std::string &help_text,
                                std::string aux, int number) {
    if (argmap.count(flag) > 0) {
        hila::out0 << "\n###################################################\n";
        hila::out0 << "# Flag " << flag << " is already set! Terminating.\n";
        hila::out0 << "###################################################\n\n";
        quit_with_help();
    } else {
        argmap[flag] = {std::vector<std::string>(), aux, help_text, number, false};
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Given flag '-flag' argv is scanned for valid entries that follow.
/// @details Once '-flag' is found the following strings are scanned. If they
///          are not in the format of a flag (described in add_flag()) they are
///          considered valid input. argv with
///          "-flag a -b_flag b -flag c d -flag -flag 4"
///          should result in {"a", "c", "d", "4"}.
///
/// @param  flag   string corresponding to the flag
/// @return Vector containing the found entries.
////////////////////////////////////////////////////////////////////////////////
strvec cmdlinearguments::read_arg_vector(const char *flag) {
    strvec uargs;

    std::vector<int> u_ind(argc);
    for (int i = 0; i < argc; i++)
        u_ind[i] = -1;
    int *p_ind = u_ind.data();

    // Anything of the format "-[a letter][whatever]" is parsed as a flag
    const std::regex forbidden_user_input("^-[a-zA-Z].*");

    for (int i = 0; i < argc; i++) {
        const char *p = argv[i];
        // check if its the flag
        int n = 0;
        bool match = false;
        if (std::strcmp(p, flag) == 0) {
            match = true;
            // Indicate that the flag has been spotted in argv
            argmap[flag].present = true;
            // Slate for removal and move onto the
            // options
            *(p_ind++) = 1;
            i++;
            if (i < argc)
                p = argv[i];
            else {
                break;
            }
            // Check if the string is a valid input
            // (Not of type ^-[a-zA-Z].*)
            // and push into vector uargs
            while (!std::regex_match(p, forbidden_user_input)) {
                *(p_ind++) = 1;
                uargs.push_back(std::string(p));
                i++;
                n++;
                if (i < argc)
                    p = argv[i];
                else
                    break;
            }

            // If we stopped on another 'flag' for whatever
            // reason, set the current index to be removed
            // and compensate for the for loop increase in
            // i and p_ind
            if (std::strcmp(p, flag) == 0) {
                *(p_ind--) = 1;
                i--;
            }
        }

        if (match) {

            if (argmap[flag].number >= 0) {
                if (argmap[flag].number != n) {
                    hila::out0 << "Error: command line flag " << flag << " requires "
                               << argmap[flag].number << " arguments\n";
                    quit_with_help();
                }
            }
        }

        p_ind++;
    }
    // Effectively remove user arguments from argv
    int j = 0;
    for (int i = 0; i < argc; i++)
        if (u_ind[i] < 0)
            argv[j++] = argv[i];
    argc = j;

    return uargs;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Loops over the known flags and fills the argmap struct from argv
////////////////////////////////////////////////////////////////////////////////
void cmdlinearguments::fill_argmap() {
    // Loop through the flags (keys of the map)
    for (auto const &p : argmap) {
        // Get corresponding input and set it to uarg[flag]
        std::vector<std::string> arg_vec = read_arg_vector(p.first.c_str());
        argmap[std::string(p.first)].val_strings = arg_vec;
    }
    if (argc > 1) {
        hila::out0 << "There remain unprocessed command-line arguments:\n";
        for (int i = 1; i < argc; i++)
            hila::out0 << "'" << argv[i] << "', ";
        hila::out0 << "\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Splits a help_text string into a vector of strings at newline.
/// @return A vector containing the lines.
////////////////////////////////////////////////////////////////////////////////
strvec cmdlinearguments::parse_help(std::string help_text) {
    strvec lines;
    std::stringstream stream(help_text);
    std::string line;

    while (std::getline(stream, line)) {
        lines.push_back(line);
    }

    return lines;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Prints the recognized flags and their help texts.
////////////////////////////////////////////////////////////////////////////////
void cmdlinearguments::print_help() {
    hila::out0 << "Recognized command-line flags and their possible arguments:\n";

    const std::string padding = "                        ";
    for (auto const &p : argmap) {
        const std::string &flag = p.first;
        const std::string &help = p.second.help_text;
        const std::string &aux = p.second.aux;
        strvec help_vec = parse_help(help);
        hila::out0 << "  " << flag << " " << aux << std::setw(22 - flag.length() - aux.length() - 1)
                   << ": " << help_vec[0] << "\n";
        for (int i = 1; i < help_vec.size(); i++) {
            hila::out0 << padding << help_vec[i] << "\n";
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Checks whether a flag has entries by returning the count. Does not
///        differentiate on whether the flag is known.
///
/// @param flag
////////////////////////////////////////////////////////////////////////////////
int cmdlinearguments::flag_set(const char *flag) {
    if (argmap.count(flag) > 0)
        return argmap[flag].val_strings.size();
    else
        return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Checks if a flag has been found in argv
///
/// @return Boolean indicating the statement.
////////////////////////////////////////////////////////////////////////////////
bool cmdlinearguments::flag_present(const char *flag) {
    return argmap[flag].present;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Attempts to return a long integer from i:th entry of flag 'flag'.
///        Failure results in terminating the program.
/// @param flag
/// @param i      index of the desired entry (default = 0);
/// @return long conversion of the desired entry
////////////////////////////////////////////////////////////////////////////////
long cmdlinearguments::get_int(const char *flag, int i) {
    // If flag is found and non-empty
    int set = flag_set(flag);
    if (set) {
        std::string opt;
        if (i < set)
            opt = argmap[flag].val_strings[i];
        else {
            hila::out0 << "Flag '" << flag << "' has only " << set + 1
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            quit_with_help();
        }

        long val = 0;
        // Check for format
        const std::regex permitted_user_input("^[+-]?[0-9]+");
        if (std::regex_match(opt, permitted_user_input)) {
            val = std::stol(opt);
        } else {
            hila::out0 << "Expected a number (integer) after command line parameter '" << flag
                       << "'\n";
            quit_with_help();
        }
        return val;
    } else {
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
        return -1;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Attempts to return a double from i:th entry of flag 'flag'.
///        Failure results in terminating the program. If std::stod
///        doesn't complain, the conversion is considered valid.
/// @param flag
/// @param i      index of the desired entry (default = 0);
/// @return double conversion of the desired entry
////////////////////////////////////////////////////////////////////////////////
double cmdlinearguments::get_double(const char *flag, int i) {
    // If flag is found and non-empty
    int set = flag_set(flag);
    if (set) {
        std::string opt;
        if (i < set)
            opt = argmap[flag].val_strings[i];
        else {
            hila::out0 << "Flag '" << flag << "' has only " << set
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            quit_with_help();
            return -1;
        }

        double val;
        // We're not going to manually check the format for this
        auto [p, ec] = std::from_chars(opt.data(), opt.data() + opt.size(), val);
        if (p == opt.data()) {
            hila::out0 << "Expected a number (double) after command line parameter '" << flag
                       << "'\n";
            quit_with_help();
        }

        // try {
        //     val = std::stod(opt);
        // } catch (std::exception &e) {
        //     hila::out0 << "Expected a number (double) after command line parameter '" << flag
        //                << "'\n";
        //     quit_with_help();
        // }

        return val;
    }
    // if not found
    else {
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
        return -1;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Attempts to return the i:th entry of flag 'flag'.
///        Failure results in terminating the program.
///
/// @param flag
/// @param i      index of the desired entry (default = 0);
/// @return the desired entry
////////////////////////////////////////////////////////////////////////////////
std::string cmdlinearguments::get_string(const char *flag, int i, bool clear) {
    int set = flag_set(flag);
    if (set) {
        if (i < set) {
            auto t = argmap[flag].val_strings[i];
            if (clear)
                argmap[flag].val_strings[i].clear();
            return t;
        } else {
            hila::out0 << "Flag '" << flag << "' has only " << set + 1
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            quit_with_help();
            return "";
        }
    } else {
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
        return "";
    }
}
