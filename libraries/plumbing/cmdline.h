#ifndef CMDLINE_H
#define CMDLINE_H

struct argmap_val
{
    std::vector<std::string> val_strings;
    std::string help_text;
};

class cmdlinearguments {
  private:
    int argc = 1;
    const char **argv;
    
    // map to hold the flag-sorted command line arguments
    std::map<std::string, argmap_val> argmap;
    
    std::vector<std::string> read_arg_vector(const char *flag);
    void fill_argmap(std::map<std::string, argmap_val> *argmap);
    std::vector<std::string> parse_help(std::string);

  public:
    cmdlinearguments();
    //cmdlinearguments(int argc0, char **argv0);
    void initialise_args(int argc0, char **argv0);
    std::vector<std::string> values(std::string);
    void add_flag(std::string flag, std::string help_text);
    void print_help();

};

extern cmdlinearguments GLOB_ARGS;

#endif
