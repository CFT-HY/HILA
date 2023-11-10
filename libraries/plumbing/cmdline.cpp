#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <vector>
#include <regex>
#include "cmdline.h"
#include "hila.h"

using strvec = std::vector<std::string>;

cmdlinearguments GLOB_ARGS;

cmdlinearguments::cmdlinearguments() {
    // Set default flags and their documentation
    add_flag("-t", "cpu time limit\n");
    add_flag("-o", "output filename (default: stdout)");
    add_flag("-i", "input filename (overrides the 1st hila::input() name)\nuse '-i -' for standard input");
    add_flag("-device","in GPU runs using only 1 GPU, choose this GPU number (default 0)");
    add_flag("-check","check input & layout with <nodes>-nodes & exit\nonly with 1 real MPI node (without mpirun)");
    add_flag("-n nodes","number of nodes used in layout check, only relevant with -check");
    add_flag("-partitions","number of partitioned lattice streams");
    add_flag("-sync","synchronize partition runs (on/off) (default=no)");
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
        hila::terminate(0);
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
        // check if user option flag
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
                hila::out0 << "Aborting run.\n";
                hila::terminate(0);
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

void cmdlinearguments::fill_argmap(std::map<std::string, argmap_val> *argmap)
{
    // Loop through the flags (keys of the map)
    for (auto const& p : *argmap)
    {
        // Get corresponding input and set it to uarg[flag]
        std::vector<std::string> arg_vec = read_arg_vector(p.first.c_str());
        (*argmap)[std::string(p.first)].val_strings = arg_vec;
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
        //<< "  -t <seconds>    : cpu time limit\n"
        //<< "  -o <name>       : output filename (default: stdout)\n"
        //<< "  -i <name>       : input filename (overrides the 1st hila::input() name)\n"
        //<< "                    use '-i -' for standard input\n"
        //<< "  -device <number>: in GPU runs using only 1 GPU, choose this GPU number (default 0)\n"
        //<< "  -check          : check input & layout with <nodes>-nodes & exit\n"
        //<< "                    only with 1 real MPI node (without mpirun)\n"
        //<< "  -n nodes        : number of nodes used in layout check, only relevant with -check\n"
        //<< "  -partitions n   : number of partitioned lattice streams\n"
        //<< "  -sync on/off    : synchronize partition runs (default=no)\n";
         
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
