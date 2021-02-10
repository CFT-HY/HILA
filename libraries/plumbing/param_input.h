#ifndef PARAM_INPUT_H
#define PARAM_INPUT_H

#include <string>
#include "defs.h"
#include "lattice.h"
#ifdef USE_MPI
#include "plumbing/com_mpi.h"
#endif

////////////////////////////////////////////////////////////////////////
/// input - Class for parsing runtime parameter files
///
/// Input files consist normally of "key <value>" -pairs (typically 1 per line).
/// <value> can be  int, double, string.
/// Whitespace or sometimes comma (,) acts as a delimiter.
/// Comment character '#': everything after # is a comment to the end of the line.
///
/// The values are read in sequentially, thus, they must be in the file in the order
/// they are read in.
///
/// Class provides the following user functions:
///
///     input in("filename");               - opens the "filename" for input
///     input in; in.open("filename");      - alternative way to open
///
///     int i = in.get("key1");
///     int i = in.get_int("key1");  (alternative, more robust form)
///          read in "key1 <integer value>" pair and return the value.
///          The first form infers the type from the assigned to variable.
///
///          Key can contain whitespace, but it has to match exactly the form in the file:
///          if file contains "nice value  3",  i = in.get("nice value");  returns 3.
///
///          key-value pair must be the next line to be read in the file, if the key
///          or integer value do not match an error is raised.
///
///     double d = in.get("key2");       or   d = in.get_double("key2");
///
///     std::string s = in.get("key3");  or   s = in.get_string("key3");
///          matches "key3 <string value>" where string value is either
///            a) sequence of non-whitespace chars, delimited by whitespace, eol, ',' or
///            '#'. b) characters enclosed by quote marks "..".  These have to pair within
///            the same line.
///               Quote marks are removed.
///
///     std::vector<int> vi = in.get("vector-key")   or  in.get_int_vector("vector-key");
///           matches line "vector-key  <int1>, <int2>, <int3>"
///           returns a std::vector<int> containing the comma-separated integer values.
///
///     std::vector<double> vd = in.get("vector-key1")        or
///     in.get_double_vector("vector-key1"); std::vector<std::string> vs =
///     in.get("vector-key2")   or  in.get_string_vector("vector-key2");
///
///     int in.get_item("key", items);
///     int in.get_item("key", items, val);     where val is (assignable) double variable
///     (*)
///          where  "std::vector<std::string> items" contains a list of allowed values.
///          Function returns the index of found item string.
///          E.g. if the file contains a line
///               animal  cat
///          the code
///               std::vector<std::string> items = {"mouse", "bird", "cat", "cow"};
///               i = in.get_item("animal", items);
///          will return value 2.  If there is no match an error is raised.
///
///          Form (*) also matches a double precision value.  If this matches, return
///          value is -1 and val is assigned to the found value.  For example, code
///               double c_sw;
///               i = in.get_item("c_sw", {"tree","perturbative"},c_sw );
///          on lines
///               c_sw  perturbative              - returns 1
///               c_sw  1.343                     - returns -1 and variable c_sw is
///               assigned 1.343 c_sw  abcd                      - gives error and quits
///
///      in.close();   - closes the file.  Now "in.open("file")" can be used again.
///          File is also closed when variable "in" goes out of scope.
///
///      NOTE: an empty key "" always matches, and can be used to read values without
///      keys.
///            E.g.
///                 i = get_int("");  or   i = get("");   reads in next integer.
///            For functions get_int(), get_double(), get_string() the empty string can be
///            left out.
///
///      NOTE: returned values are always broadcast to all MPI nodes.  No need for the
///      user to do it.
///
////////////////////////////////////////////////////////////////////////

class input {
  public:
    input(){};
    ~input() { close(); }
    input(const std::string &fname) { open(fname); };

    // Trick to "specialize" to return type
    class returntype {
      public:
        const std::string &label;
        input *parent;

        returntype(const std::string &str, input *a) : label(str), parent(a) {}

        operator double() { return parent->get_double(label); }
        operator float() { return (float)parent->get_double(label); }

        operator int() { return parent->get_int(label); }
        operator std::string() { return parent->get_string(label); }

        operator std::vector<int>() { return parent->get_int_vector(label); }

        operator std::vector<double>() { return parent->get_double_vector(label); }

        operator std::vector<std::string>() { return parent->get_string_vector(label); }
    };

    returntype get(const std::string &label);
    returntype get();

    void open(const std::string &fname);

    double get_double(const std::string &label);
    int get_int(const std::string &label);
    std::string get_string(const std::string &label);
    double get_double();
    int get_int();
    std::string get_string();

    int get_item(const std::string &label, const std::vector<std::string> &items,
                 double *dval = nullptr);
    int get_item(const std::string &label, const std::vector<std::string> &items,
                 double &val);
    std::vector<int> get_int_vector(const std::string &label);
    std::vector<double> get_double_vector(const std::string &label);
    std::vector<std::string> get_string_vector(const std::string &label);

    void close();

  private:
    std::ifstream inputfile;
    bool is_initialized = false;
    std::string filename;

    std::string linebuffer;
    size_t lb_start = 0; // linebuffer start index
    bool is_line_printed;

    // these store the values scanned from list items
    int item_int_val = 0;
    double item_double_val = 0;
    std::string item_string_val;

    bool peek_token(std::string &tok);
    bool get_token(std::string &tok);
    bool match_token(const std::string &tok);

    bool scan_string(std::string &val);
    bool is_value(const std::string &s, double &val);
    bool is_value(const std::string &s, int &val);
    bool is_value(const std::string &s, std::string &val);
    bool contains_word_list(const std::string &list, int &end_of_key);

    std::string remove_quotes(const std::string &val);

    void print_linebuf(int eok);

    bool get_line();
    bool handle_key(const std::string &c);

    bool remove_whitespace();

    template <typename T>
    void get_type_vector(const std::string &label, std::vector<T> &res, const char *name);
};

#endif
