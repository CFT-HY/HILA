#ifndef PARAM_INPUT_H
#define PARAM_INPUT_H

#include<string>
#include "defs.h"
#include "lattice.h"
#ifdef USE_MPI
#include "plumbing/com_mpi.h"
#endif


////////////////////////////////////////////////////////////////////////
/// input - Class for parsing runtime parameter files using std c++ libraries
/// 
/// Fulfills three simple functions: 
///
/// 1. Parses a text file or command line for runtime variables
/// 2. Allows user to define which variables should be found in parameter files
/// 3. Allows the user to assign these runtime variables by name 
///    using the following syntax:
///     
///    input input1
///    input1.import("params.txt")
///    int nx = input1.get("nx")
///    string out = input1.get("outputfname")
///
/// When using MPI, the input data is read by one process and broadcasted to
/// the others.      
////////////////////////////////////////////////////////////////////////

class input {
    public:

        input(){};
        ~input() { close(); }
        input(const std::string & fname) { init(fname); };

        // Trick to "specialize" to return type
        class returntype {
            public:
            const std::string & label;
            input * parent;

            returntype(const std::string & str, input * a) : label(str), parent(a) {}

            operator double() {
                return parent->get_double(label);
            }
            operator float() {
                return (float) parent->get_double(label);
            }

            operator int() {
                return parent->get_int(label);
            }
            operator std::string() {
                return parent->get_string(label);
            }

            operator std::vector<int>() {
                return parent->get_int_list(label);
            }

            operator std::vector<double>() {
                return parent->get_double_list(label);
            }

            operator std::vector<std::string>() {
                return parent->get_string_list(label);
            }

        };

        void init(const std::string &fname);
        returntype get(const std::string &label);
        double get_double(const std::string &label);
        int    get_int(const std::string &label);
        std::string get_string(const std::string &label);

        int get_item(const std::string &label, const std::vector<std::string> &items);
        std::vector<int> get_int_list(const std::string & label);
        std::vector<double> get_double_list(const std::string & label);
        std::vector<std::string> get_string_list(const std::string & label);

        void close();


    private:

        std::ifstream inputfile;
        bool is_initialized = false;
        std::string filename;

        std::string linebuffer;
        size_t lb_start = 0;   // linebuffer start index

        // these store the values scanned from list items
        int item_int_val = 0;     
        double item_double_val = 0;
        std::string item_string_val;

        bool peek_token( std::string & tok ); 
        bool get_token( std::string & tok );
        bool match_token(const std::string & tok);
 
        bool scan_string(std::string &val);
        bool is_value(const std::string &s, double &val);
        bool is_value(const std::string &s, int &val);
        bool is_value(const std::string &s, std::string & val);

        std::string remove_quotes(const std::string &val);

        void print_linebuf();

        bool get_line();
        bool handle_label(const std::string & c);

        template <typename T>
        void get_type_list(const std::string & label, std::vector<T> & res, 
                           const char * name);

     
};

#endif
