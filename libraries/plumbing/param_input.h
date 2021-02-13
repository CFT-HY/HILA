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
/// Input files consist normally of "key <value>" -pairs.
///
/// Comment character '#': everything after # is a comment to the end of the line.
///
/// The values are read in sequentially, thus, they must be in the file in the order
/// they are read in.
///
/// Class provides user functions 
///    open(), close(), get(), get_item(), get_value(), quiet()
///
///
///       input f("filename");        - opens the "filename" for input
///
///       input f; 
///       f.open("filename");         - this also opens the file 
///
/// get() is the main method for reading input:
///
///       var = f.get("key");      
///
///    reads in a key-value pair
///             key   <value(s)>
///    from the input file, and returns the value of type of variable var.
///    The value is broadcast to all MPI nodes.  The method infers
///    the type of the returned variable from the type of the assingnment.
///
///    Key is alphanumeric string, which may contain words separated by whitespace.
///    The match has to be exact, except for the amount of whitespace.
///        
///    Recognized types:
///        int, long, float, double, Cmplx<float>, Cmplx<double>, 
///        std::string, 
///        Vector<n,T>, std::vector<T>     where T is one of the above types
///
///  If there is an error, an error message is printed and the program quits.
///
///  Examples:
///
///        int i = f.get("number of cats");
///
///    matches " number of cats   5 "  from the file f, and returns 5 on
///    variable i, which is broadcast to all nodes.
///
///        CoordinateVector v;
///        v = f.get("lattice size");
///
///    matches "lattice  size     32, 32, 32, 32"  (if NDIM == 4).
///    CoordinateVector is a derived type of Vector<NDIM,int>.
///    Multiple items are separated by commas, amount of whitespace is not
///    significant.
///
///        std::vector<double> dvec = f.get("vec");
///
///    matches "vec  3, 4,5.5 , 7.8, 4"  and returns a vector of double values.
///    The numbers are read until they are not followed by a comma.  If comma is
///    a last caracter on a line, reading continues on the next line.
///
///        Cmplx<double> phase = f.get("complex phase");
///
///    matches "complex phase   (0.4, 0.5)"
///    complex values are given in pairs within ( , )
///
///        std::string s = f.get("key");
///    
///    matches "key3 <string value>" where string value is either
///     a) sequence of non-whitespace chars, delimited by whitespace, eol, ',' or '#'.
///     b) characters enclosed by quote marks "..".  These have to pair within
///        the same line.  Quote marks are removed.
///
///    If there is no key label, the value is read without requiring any key:
///
///         int i = f.get();    // read an int
///
///
///
/// get_item() is used to select one from a "menu" of alternatives:     
///
///         int i;
///         i = f.get_item("key", items);
///         i = f.get_item("key", items, val);  (*) where val is (assignable) 
///                                                 double variable
///     
///    where  "std::vector<std::string> items" contains a list of allowed items.
///    Function returns the index of found item string.
///    E.g. if the file contains a line
///          "animal  cat"
///    calling
///         std::vector<std::string> items = {"mouse", "bird", "cat", "cow"};
///         i = f.get_item("animal", items);
///    will return value 2.  If there is no match an error is raised.
///
///    Form (*) also matches a double precision value.  If this matches, return
///    value is -1 and val is assigned to the found value.  For example, code
///         double c_sw;
///         i = in.get_item("c_sw", {"tree","perturbative"},c_sw );
///    and file content: 
///         c_sw  perturbative   - returns 1
///         c_sw  1.343          - returns -1 and variable c_sw <- 1.343
///         c_sw  abcd           - error message and quit
///
///
/// close(): 
///
///    closes the input f.  Now "f.open("file")" can be used again.
///    File is also closed when variable "f" goes out of scope.
///
///    NOTE: an empty key "" always matches, get("") is equivalent to get(),
///
///    NOTE: get() sends the data to all nodes, and quits on error.  Different
///          interface is provided by get_value():
///
/// get_value():
///
///         template <typename T>
///         bool get_value(T & val,std::string key, bool broadcast=true);
///
///    val can be any value used in get()-method above.  If broadcast==false, the 
///    value is not broadcast to other nodes.  The return value is false if
///    the value could not be read successfully, true otherwise.  
///    This method does not exit on error (but an error message is printed)
///    Example:
///          int i,j;
///          bool success;
///          success = get_value(i, "key", false);   // only node 0 gets i
///          success = get_value(j,"key2");          // all nodes get j
///  
///    NOTE: if broadcast == false the return value is correct only on node 0.
///          Synchronization must be done afterwards
///
///
/// quiet()  and  quiet(false):
///
///    By default the methods print everything read to standard output, for 
///    logging.   f.quiet()  disables this, f.quiet(false) re-enables.
///
///
////////////////////////////////////////////////////////////////////////




class input {

  private:

    std::ifstream inputfile;
    bool is_initialized = false;
    std::string filename;

    std::string linebuffer;
    size_t lb_start = 0; // linebuffer start index
    bool is_line_printed;
    bool speaking = true;          // false does not print information

  public:

    input(){}
    ~input() { close(); }
    input(const std::string &fname) {  open(fname); }

    bool open(const std::string &fname, bool exit_on_error = true);
    void close();

    // make class quiet (no printouts), quiet(false) returns to normal
    void quiet(bool really=true) { speaking = !really; }
 

    // Trick to "specialize" .get("label") -method to return type
    class returntype {
      public:
        const std::string &label;
        input *parent;

        returntype(const std::string &str, input *a) : label(str), parent(a) {}

        template <typename T>
        operator T() { 
            T val;
            if (!parent->get_value(val,label,true))
                hila::finishrun();
            return val; 
        }

        // operator int() { return (int)parent->get_value<long>(label); }

        // operator long() { return parent->get_value<long>(label); }

        // operator float() { return (float)parent->get_value<double>(label); }

        // operator double() { return parent->get_value<double>(label); }

        // operator Cmplx<double>() { return parent->get_value<Cmplx<double>>(label); }

        // operator Cmplx<float>() { return (Cmplx<float>)parent->get_value<Cmplx<double>>(label); }

        // operator CoordinateVector() { return parent->get_value<CoordinateVector>(label); }

        // operator std::vector<double>() { return parent->get_vector<double>(label); }

        // operator std::vector<int>() { return parent->get_vector<int>(label); }

        // operator std::vector<Cmplx<double>>() { return parent->get_vector<Cmplx<double>>(label); }


    };

    // The main get() method is simply constructor for returntype

    inline returntype get(const std::string &key) { 
        return returntype(key, this);
    }

    inline returntype get() { 
        return returntype("", this);
    }

    /// General single-value input method, can be called with
    /// <input>.get_value<type>("key");
    /// but typically used in .get() -methods

    template <typename T>
    bool get_value(T & val, const std::string &label, bool bcast=true) {
        val = {};
        bool no_error = handle_key(label); // removes whitespace

        if (hila::myrank() == 0 && no_error) {
            std::string tok;
            if (!(get_token(tok) && is_value(tok, val))) {
                
                if (speaking)
                    hila::output << "Error: expecting " << type_id<T>() << " after '" << label << "'\n";
    
                no_error = false;
            }
        }

        if (bcast) {
            // string has to be treated separately
            if constexpr (std::is_same<T,std::string>::value) {
                broadcast(val);
                broadcast(no_error);
            } else {
                broadcast(val,no_error);
            }
        }

        return no_error;
    }

    /// Specialize the above method to Cmplx -pair:  (re,im)

    template <typename T>
    bool get_value(Cmplx<T> & val, const std::string &label, bool bcast=true) {
        val = 0;
        bool no_error = handle_key(label); // removes whitespace

        if (hila::myrank() == 0 && no_error) {
            std::string tok;
            T re, im;

            no_error = ( match_token("(") &&
                         get_token(tok) && is_value(tok,re) &&
                         match_token(",") &&
                         get_token(tok) && is_value(tok,im) &&
                         match_token(")") );
            if (!no_error && speaking) {
                hila::output << "Error: expecting complex value '(re,im)' after '" << label << "'\n";
            }

            val = Cmplx<T>(re,im);
        }

        if (bcast) {
            broadcast(val, no_error);
        }
        return no_error;
    }


    /// Specialize .get_value<Vector<n,T>>() : which includes CoordinateVector

    template <int n, typename T>
    bool get_value(Vector<n,T> & val, const std::string &label, bool bcast=true) {
        val = 0;
        bool no_error = true;

        if (hila::myrank() == 0) {
            no_error = get_value(val[0],label,false);
            for (int i=1; i<n && no_error; i++) {
                no_error = get_value(val[i],",",false);
            }

            if (!no_error && speaking) {
                hila::output << "Error: expecting " << NDIM << " comma-separated " 
                          << type_id<T>() << "s after '" << label << "'\n";
            }
        }

        if (bcast) {
            broadcast(val,no_error);
        }
        return no_error;
    }

    /// Specialize -get_value() to std::vector<> 

    template <typename T>
    bool get_value(std::vector<T> & val, const std::string &label, bool bcast=true) {
        val = {};
        bool no_error = true;

        if (hila::myrank() == 0) {
            T v;
            no_error = get_value(v,label,false);
            val.push_back(v);
            while (no_error && match_token(",")) {
                no_error = get_value(v,"",false);
                val.push_back(v);
            }

            if (!no_error && speaking) {
                hila::output << "Error: expecting a comma-separated list of "
                          << type_id<T>() << "s after '" << label << "'\n";
            }
        }

        if (bcast) {
            broadcast(no_error);
            broadcast(val);
        }

        return no_error;
    }

    // get_item selects one from a "menu" of items

    int get_item(const std::string &label, const std::vector<std::string> &items,
                 bool bcast=true);

  private:

    /// a helper method to give type name
    template <typename T>  
    inline const char * type_id() { return nullptr; }


    bool peek_token(std::string &tok);
    bool get_token(std::string &tok);
    bool match_token(const std::string &tok);

    bool scan_string(std::string &val);

    bool is_value(const std::string &s, int &val);
    bool is_value(const std::string &s, long &val);
    bool is_value(const std::string &s, std::string &val);
    bool is_value(const std::string &s, double &val);
    bool is_value(const std::string &s, float &val);


    bool contains_word_list(const std::string &list, int &end_of_key);

    std::string remove_quotes(const std::string &val);

    void print_linebuf(int eok);

    bool get_line();
    bool handle_key(const std::string &c);

    bool remove_whitespace();

};

/// give here specializations of the type_id helper

template <> inline const char * input::type_id<int>() { return "int"; }
template <> inline const char * input::type_id<long>() { return "long"; }
template <> inline const char * input::type_id<float>() { return "float"; }
template <> inline const char * input::type_id<double>() { return "double"; }
template <> inline const char * input::type_id<std::string>() { return "string"; }
template <> inline const char * input::type_id<Cmplx<float>>() { return "complex value"; }
template <> inline const char * input::type_id<Cmplx<double>>() { return "complex value"; }




#endif
