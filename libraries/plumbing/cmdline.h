#ifndef HILA_CMDLINE_H_
#define HILA_CMDLINE_H_

#include <map>
#include <vector>
#include <string>

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
    void fill_argmap();

  public:
    cmdlinearguments();
    void initialise_args(int argc0, char **argv0);
    std::vector<std::string> values(const std::string &);
    void add_flag(const std::string &flag, const std::string &help_text,
                  const std::string aux = "", int number = -1);
    void print_help();
    int flag_set(const char *flag);
    bool flag_present(const char *flag);
    long get_int(const char *flag, int i = 0);
    double get_double(const char *flag, int i = 0);
    std::string get_string(const char *flag, int i = 0, bool clear = false);
};

namespace hila {
extern cmdlinearguments cmdline;
}

#endif
