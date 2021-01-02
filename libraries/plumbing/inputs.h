
#ifndef INPUTS
#define INPUTS

#include<string>
#include<map>

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
        ~input(){};

        //read runtime parameters from file
        void import(const std::string & fname);
        //read runtime parameters from file and cmd line 
        void import(int & argc, char *** argvp, const std::string & fname);

        //add an essential variable - triggers complaint if not in parameter file or commandline
        void add_essential(const std::string & );
        //same as above, except defines a default value to be used in case it is not found
        template< typename T >
        void add_essential(std::string const &, T const & default_value);

        /// The return type of input parameters. Avoids the problem of specializing by
        /// return type: the returned value caan be cast to double, float, int or 
        /// std::string.
        class returntype {
            public:
            returntype(const std::string & str, input * a) : name(str), parent(a) {}
            const std::string & name;
            input * parent;
            operator double(){
                return parent->values[name];
            }
            operator float(){
                return (float) parent->values[name];
            }
            operator int(){
                return (int) parent->values[name]; 
            }
            operator std::string(){
                return parent->names[name];
            }
        };

        returntype get(const std::string &);
        void close();

    protected:

        std::map<std::string, double> values;
        std::map<std::string, std::string> names;

    private:

        std::map<std::string, bool> essentials;

        void check_essentials();
        void handle(const std::string &);
        void define_essentials();
        void read(const std::string &);

        //internal state broadcast in mpi implementation
        #ifdef USE_MPI
        void broadcast_values();
        void broadcast_names();
        #endif
};

#endif
