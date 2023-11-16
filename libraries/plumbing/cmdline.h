#ifndef CMDLINE_H
#define CMDLINE_H



class cmdlinearguments {
  private:
    int argc;
    const char **argv;
    
    // map to hold the flag-sorted command line arguments
    // argmap_val is defined in cmdline.cpp
    std::map<std::string, struct argmap_val> argmap;
    
    std::vector<std::string> read_arg_vector(const char *flag);
    std::vector<std::string> parse_help(std::string);
    void quit_with_help();

  public:
    cmdlinearguments();
    void fill_argmap();
    void initialise_args(int argc0, char **argv0);
    std::vector<std::string> values(std::string);
    void add_flag(std::string flag, std::string help_text);
    void print_help();
    int flag_set(const char *flag);
    bool flag_used(const char *flag);
    long get_int(const char *flag, int i = 0);
    double get_double(const char *flag, int i = 0);
    std::string get_string(const char *flag, int i = 0);

};

namespace hila{
extern cmdlinearguments cmdline;
}

#endif
