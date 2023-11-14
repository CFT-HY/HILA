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
    }

void cmdlinearguments::initialise_args(int argc0, char **argv0) {
    argc = argc0;
    argv = (const char **)malloc(argc * sizeof(const char *));
    for (int i = 0; i < argc; i++)
        argv[i] = argv0[i];
}

strvec cmdlinearguments::values(std::string flag)
{
    strvec valvec = argmap[flag].val_strings;
    if (argmap.count(flag) == 0)
    {
        hila::out0 << "Flag '" << flag << "' is not recognized!\n";
        quit_with_help();
    }

    if (valvec.size() == 0)
    {
        hila::out0 << "\n\nFlag '" << flag << "' has no entries"
                   << " and the associated vector of strings is of length zero!\n\n";
    }
    return valvec;
}

void cmdlinearguments::add_flag(std::string flag, std::string help_text)
{
    if (argmap.count(flag) > 0)
    {
        hila::out0 << "\n###################################################\n";
        hila::out0 << "# Flag " << flag << " is already set! Terminating.\n";
        hila::out0 << "###################################################\n\n";
        quit_with_help();
    }
    else
    {
        argmap[flag] = {std::vector<std::string>(), help_text, false};
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
            // Indicate that the flag has been spotted on the command-line
            argmap[flag].used = true;
            // Slate for removal and move onto the
            // options
            *(p_ind++) = 1;
            i++;
            if (i < argc)
                    p = argv[i];
            else
            {
                break;
            }
            // Check if the string is a valid input
            // (Not of type ^-[a-zA-Z].*)
            // and push into vector uargs
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
    if (argc > 1)
    {
        hila::out0 << "There remain unprocessed command-line arguments:\n";
        for (int i = 1; i < argc; i++)
            hila::out0 << "'" << argv[i] << "', ";
        hila::out0 << "\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
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

bool cmdlinearguments::flag_used(const char *flag)
{
    return argmap[flag].used;
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
            quit_with_help();
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
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
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
            hila::out0 << "Flag '" << flag << "' has only " << set
                       << " entries. Can't return index " << i << ".\n";
            hila::out0 << "Terminating.\n";
            quit_with_help();
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
    else
    {
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
        return -1;
    }
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
            quit_with_help();
            return "";
        }
    }
    else
    {
        hila::out0 << "Flag '" << flag << "' is missing an entry!\n";
        hila::out0 << "Terminating.\n";
        quit_with_help();
        return "";
    }
}

