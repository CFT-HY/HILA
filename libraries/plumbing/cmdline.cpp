#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <regex>
#include <limits.h>
#include "cmdline.h"
#include "hila.h"

using strvec = std::vector<std::string>;

void cmdlinearguments::quit_with_help()
{
    GLOB_ARGS.print_help();
    hila::finishrun();
}

cmdlinearguments GLOB_ARGS;

cmdlinearguments::cmdlinearguments() {
    // Set default flags and their documentation when initializing
    add_flag("-t","cpu time limit\n");
    add_flag("-o","output filename (default: stdout)");
    add_flag("-i","input filename (overrides the 1st hila::input() name)\nuse '-i -' for standard input");
    add_flag("-device","in GPU runs using only 1 GPU, choose this GPU number (default 0)");
    add_flag("-check","check input & layout with <nodes>-nodes & exit\nonly with 1 real MPI node (without mpirun)");
    add_flag("-n","number of nodes used in layout check, only relevant with -check");
    add_flag("-partitions","number of partitioned lattice streams");
    add_flag("-sync","synchronize partition runs (on/off) (default = off)");
}

void cmdlinearguments::initialise_args(int argc0, char **argv0) {
    argc = argc0;
    argv = (const char **)malloc(argc * sizeof(const char *));
    for (int i = 0; i < argc; i++)
        argv[i] = argv0[i];
}

strvec cmdlinearguments::values(std::string flag)
{
    return argmap[flag].val_strings;
}

void cmdlinearguments::add_flag(std::string flag, std::string help_text)
{
    if (argmap.count(flag) > 0)
    {
        hila::out0 << "###################################################\n";
        hila::out0 << "# Flag " << flag << " is already set! Terminating.\n";
        hila::out0 << "###################################################\n";
        quit_with_help();
    }
    else
    {
        argmap[flag] = {std::vector<std::string>(), help_text};
    }
}

strvec cmdlinearguments::read_arg_vector(const char *flag)
{
    strvec uargs;

    int u_ind[argc];
    for (int i = 0; i < argc; i++) u_ind[i] = -1;
    int *p_ind = u_ind;

    // Anything of the format "-[a letter][whatever]" is parsed as a flag
    const std::regex forbidden_user_input("^-[a-zA-Z].*");

    for (int i = 0; i < argc; i++)
    {
        const char *p = argv[i];
        // check if its the flag
        if (std::strcmp(p, flag) == 0)
        {
            // Slate for removal and move onto the
            // options
            *(p_ind++) = 1;
            i++;
            if (i < argc)
                    p = argv[i];
            else
            {
                hila::out0 << "Flag " << flag << " is missing an entry!\n";
                hila::out0 << "Terminating.\n";
                hila::finishrun();
                break;
            }
            // Check if the string is a valid input
            // (Not of type ^-[a-zA-Z].*)
            // and push into vector uargs
            //while (p[0] != '-')
            while (!std::regex_match(p, forbidden_user_input))
            {
                *(p_ind++) = 1;
                uargs.push_back(std::string(p));
                i++;
                if (i < argc)
                    p = argv[i];
                else
                    break;
            }
            // If we stopped on another 'flag' for whatever
            // reason, set the current index to be removed
            // and compensate for the for loop increase in
            // i and p_ind
            if (std::strcmp(p, flag) == 0)
            {
                *(p_ind--) = 1;
                i--;
            }
        }
        p_ind++;
    }
    // Effectively remove user arguments from argv
    int j = 0;
    for (int i = 0; i < argc; i++) if (u_ind[i] < 0) argv[j++] = argv[i];
    argc = j;

    return uargs;
}

void cmdlinearguments::fill_argmap()
{
    // Loop through the flags (keys of the map)
    for (auto const& p : argmap)
    {
        // Get corresponding input and set it to uarg[flag]
        std::vector<std::string> arg_vec = read_arg_vector(p.first.c_str());
        argmap[std::string(p.first)].val_strings = arg_vec;
    }
}

strvec cmdlinearguments::parse_help(std::string help_text)
{
    strvec lines;
    std::stringstream stream(help_text);
    std::string line;

    while (std::getline(stream, line))
    {
        lines.push_back(line);
    }

    return lines;
}

void cmdlinearguments::print_help()
{
    hila::out0 << "Recognized command-line flags:\n";

    for (auto const& p : argmap)
    {
        std::string flag = p.first;
        std::string help = p.second.help_text;
        strvec help_vec = parse_help(help);
        hila::out0 << "    " << flag << std::setw(20 - flag.length())
                   << ": " << help_vec[0] << "\n";
        for (int i = 1; i < help_vec.size(); i++)
        {
            std::string padding = "                        ";
            hila::out0 << padding << help_vec[i] << "\n";
        }
    }

}

int cmdlinearguments::flag_set(const char *flag)
{
    if (argmap.count(flag) > 0)
        return argmap[flag].val_strings.size();
    else return 0;
}

long cmdlinearguments::get_int(const char *flag, int i)
{
    // If flag is found and non-empty
    int set = flag_set(flag);
    if (set)
    {
        std::string opt;
        if (i < set)
            opt = argmap[flag].val_strings[i];
        else {
            hila::out0 << "Flag '" << flag << "' has only " << set + 1
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            hila::finishrun();
        }

        long val = 0;
        // Check for format
        const std::regex permitted_user_input("^[+-]?[0-9]+");
        if (std::regex_match(opt, permitted_user_input)) {
            val = std::stol(opt);
        }
        else
        {
            hila::out0 << "Expected a number (integer) after command line parameter '" << flag
                       << "'\n";
            quit_with_help();
        }
        return val;
    }
    else {
        hila::out0 << "Requested flag '" << flag << "' is missing an entry!\n";
        return LONG_MAX;
    }
}

double cmdlinearguments::get_double(const char *flag, int i)
{
    // If flag is found and non-empty
    int set = flag_set(flag);
    if (set)
    {
        std::string opt;
        if (i < set)
            opt = argmap[flag].val_strings[i];
        else {
            hila::out0 << "Flag '" << flag << "' has only " << set + 1
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            hila::finishrun();
            return -1;
        }

        double val;
        // We're not going to manually check the format for this
        try
        {
            val = std::stod(opt);
        }
        catch(std::exception &e)
        {
            hila::out0 << "Expected a number (double) after command line parameter '" << flag
                       << "'\n";
            quit_with_help();
        }
        return val;
    }
    // if not found
    else return -1;
}

std::string cmdlinearguments::get_string(const char *flag, int i)
{
    int set = flag_set(flag);
    if (set)
    {
        if (i < set)
            return argmap[flag].val_strings[i];
        else {
            hila::out0 << "Flag '" << flag << "' has only " << set + 1
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            hila::finishrun();
            return "";
        }
    }
    else
    {
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        hila::finishrun();
        return "";
    }
}

