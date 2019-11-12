
#ifndef INPUTS
#define INPUTS

#include<string>
#include<map>
#include<list>

////////////////////////////////////////////////////////////////////////
/// input - Class for parsing runtime parameter files using std c++ libraries
/// 
/// Fulfills two simple functions: 
///
/// 1. parses a text file for variables
/// 2. allows the programmer to assign these runtime variables by name 
///    using the following syntax:
///     
///    input input1("")
///    int nx = input1.get("nx")
///    string out = input1.get("outputfname")
///   
///    TODO: edit the constructor so that some root process reads the parameters,
///    and then speads the data amongst different processes. Ie. the root reads
///    and processes, then sends the input data to all processes. 
////////////////////////////////////////////////////////////////////////

class input {
    public:

        //read runtime parameters from file
        input(const std::string & fname);
        ~input(){};
    
        std::map<std::string, double> values;
        std::map<std::string, std::string> names;

        //add an essential variable - triggers complaint if not in parameter file
        void add_essential(const std::string & );
        template< typename T >
        //same as above, except defines a default value to be used in case it is not found
        void add_essential(const std::string &, const T & default_value);

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

    private:

        std::map<std::string, bool> essentials;
        void check_essentials();
        //handle one line of input
        void handle(const std::string &);
        //define necessary runtime variables
        void define_essentials();

        void read(const std::string &);
        
        //internal state broadcast in mpi implementation
        #ifdef USE_MPI
        void broadcast_values();
        void broadcast_names();
        #endif
};

#endif
