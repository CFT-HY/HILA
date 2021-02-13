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
/// Comment character '#': everything after # is a comment to the end of the
/// line.
///
/// The values are read in sequentially, thus, they must be in the file in the
/// order they are read in.
///
/// Class provides user functions
///    open(), close(), get(), get_item(), get_value(), quiet()
///
///
///    input f("filename");   - initialization with filename opens
///                             the file for input
///    input f;
///
/// open(): open file for reading
///
///      bool open(std::string filename, bool exit_on_error=true )
///
///    If exit_on_error == true, quit the
///    program on error.
///
///      f.open("filename");
///      bool success = f.open("filename",false);
///
///
/// get(std::string key) - read input values
///
///       var = f.get("key");
///
///    reads in a key-value pair
///             key   <value(s)>
///    from the input file, and returns the value of type of variable var.
///    The value is broadcast to all MPI nodes.  The method infers
///    the type of the returned variable from the type of the assignment.
///
///    Key is an alphanumeric string, which may contain words separated by
///    whitespace.
///
///    Recognized types:
///        int, long, float, double, Cmplx<float>, Cmplx<double>,
///        std::string, CoordinateVector,
///        Vector<n,T>, std::vector<T>     where T is one of the above types
///
///    If there is an error, an error message is printed and the program quits.
///
///    Examples:
///
///        int i = f.get("number of cats");
///
///    matches " number of cats   5 "  from the file f, and returns 5 on
///    variable i, which is broadcast to all nodes.
///
///        CoordinateVector v;
///        v = f.get("lattice size");
///
///    matches "lattice  size  32, 32, 32, 32"  (if NDIM == 4).
///
///    Multiple items are separated by commas, whitespace is not
///    significant.
///
///        std::vector<double> dvec = f.get("vec");
///
///    matches "vec  3,4, 5.5, 7.8, 4"  and returns a vector of double values.
///    The numbers are read until they are not followed by a comma.  If comma is
///    the last caracter on a line, reading continues to the next line.
///
///        Cmplx<double> phase = f.get("complex phase");
///
///    matches "complex phase   (0.4, 0.5)"
///    complex values are given in pairs within ( , )
///
///        std::string s = f.get("key");
///
///    matches "key3 <string value>" where string value is either
///     a) sequence of non-whitespace chars, delimited by whitespace, eol, ','
///        or '#'.
///     b) characters enclosed by quotes "..".  These have to pair
///        within the same line.  Quote marks are removed.
///
///    If there is no key label, the value is read without requiring any key:
///
///         int i = f.get();    // read an int
///
/// get_value(): read input (alternative to get())
///
///         template <typename T>
///         bool get_value(T & val,std::string key, bool broadcast=true);
///
///    Val can be any value used in get()-method above.  If broadcast==false,
///    the value is not broadcast to other nodes.  The return value is false if
///    the value could not be read successfully, true otherwise.
///    This method does not exit on error (but an error message may be printed)
///    Example:
///          int i,j;
///          bool success;
///          success = get_value(i, "key", false);   // only node 0 gets i
///          success = get_value(j,"key2");          // all nodes get j
///
///    NOTE: if broadcast == false the return value is correct only on node 0.
///
///
///
/// get_item(): select one item from a "menu":
///
///         int get_item(std::string key, std::vector<std::string> items,
///                      bool broadcast = true);
///
///    "items" contains the allowed entries. Return value is the
///    index of the item found in input file.
///
///    If the value of the optional bool parameter broadcast is:
///    - true (default): the result is broadcast to all nodes
///                      and the program exits if no matches found.
///    - false: result is not broadcast, and if no match found returns -1.
///
///    Special item values:
///      "%f"  - matches a float or double value
///      "%i"  - matches an int or long
///      "%s"  - matches any string value
///
///    If one of these is matched, it has to be read again with corresponding
///    get() or get_value() -method.
///
///    Examples:
///
///         i = f.get_item("colour", {"red","green","blue"});
///
///    will return value 1 if f contains "colour  green", and quits the
///    program if none of the 3 alternatives are found.
///
///         double clover;
///         int i = f.get_item("c_sw", {"tree","perturbative","%f"} );
///         if (i == 2)
///              clover = f.get();
///         else { ...
///
///    If file contains:
///         c_sw  perturbative   - get_item() returns 1
///         c_sw  1.343          - get_item() returns 2 and subsequent
///                                get() sets c_sw = 1.343
///         c_sw  abcd           - error message and quit
///
///    NOTE: "%s" matches anything. It should be the last item in the list.
///         (The items are tested in order and first to match is returned.)
///
///
/// close():
///
///    closes the input f.  Now "f.open("file")" can be used again.
///    File is also closed when variable "f" goes out of scope.
///
///
/// quiet()  and  quiet(false):
///
///    By default the methods print everything read to standard output, for
///    logging.   f.quiet()  disables this, f.quiet(false) re-enables.
///
///
/// NOTE: methods which broadcast to all nodes (default) must be called
///       from all nodes synchronously.  Thus,
///
///        if (hila::myrank() == 0) {
///            ...
///            double v = f.get("a value");
///            ...
///        }
///
///    fails.  get_value() with broadcast=false can be used in this case.
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
    bool speaking = true; // false does not print information

  public:
    input() {}
    ~input() { close(); }
    input(const std::string &fname) { open(fname); }

    bool open(const std::string &fname, bool exit_on_error = true);
    void close();

    // make class quiet (no printouts), quiet(false) returns to normal
    void quiet(bool really = true) { speaking = !really; }

    // Trick to "specialize" .get("label") -method to return type
    class returntype {
      public:
        const std::string &label;
        input *parent;

        returntype(const std::string &str, input *a) : label(str), parent(a) {}

        template <typename T> operator T() {
            T val;
            if (!parent->get_value(val, label, true))
                hila::finishrun();
            return val;
        }

        // operator int() { return (int)parent->get_value<long>(label); }

        // operator long() { return parent->get_value<long>(label); }

        // operator float() { return (float)parent->get_value<double>(label); }

        // operator double() { return parent->get_value<double>(label); }

        // operator Cmplx<double>() { return
        // parent->get_value<Cmplx<double>>(label); }

        // operator Cmplx<float>() { return
        // (Cmplx<float>)parent->get_value<Cmplx<double>>(label); }

        // operator CoordinateVector() { return
        // parent->get_value<CoordinateVector>(label); }

        // operator std::vector<double>() { return
        // parent->get_vector<double>(label); }

        // operator std::vector<int>() { return parent->get_vector<int>(label);
        // }

        // operator std::vector<Cmplx<double>>() { return
        // parent->get_vector<Cmplx<double>>(label); }
    };

    // The main get() method is simply constructor for returntype

    inline returntype get(const std::string &key) {
        return returntype(key, this);
    }

    inline returntype get() { return returntype("", this); }

    /// General single-value input method, can be called with
    /// <input>.get_value<type>("key");
    /// but typically used in .get() -methods

    template <typename T>
    bool get_value(T &val, const std::string &label, bool bcast = true) {
        val = {};
        bool no_error = handle_key(label); // removes whitespace

        if (hila::myrank() == 0 && no_error) {
            std::string tok;
            if (!(get_token(tok) && is_value(tok, val))) {

                if (speaking)
                    hila::output << "Error: expecting " << type_id<T>()
                                 << " after '" << label << "'\n";

                no_error = false;
            }
        }

        if (bcast) {
            // string has to be treated separately
            if constexpr (std::is_same<T, std::string>::value) {
                broadcast(val);
                broadcast(no_error);
            } else {
                broadcast(val, no_error);
            }
        }

        return no_error;
    }

    /// Specialize the above method to Cmplx -pair:  (re,im)

    template <typename T>
    bool get_value(Cmplx<T> &val, const std::string &label, bool bcast = true) {
        val = 0;
        bool no_error = handle_key(label); // removes whitespace

        if (hila::myrank() == 0 && no_error) {
            std::string tok;
            T re, im;

            no_error =
                (match_token("(") && get_token(tok) && is_value(tok, re) &&
                 match_token(",") && get_token(tok) && is_value(tok, im) &&
                 match_token(")"));
            if (!no_error && speaking) {
                hila::output
                    << "Error: expecting complex value '(re,im)' after '"
                    << label << "'\n";
            }

            val = Cmplx<T>(re, im);
        }

        if (bcast) {
            broadcast(val, no_error);
        }
        return no_error;
    }

    /// Specialize .get_value<Vector<n,T>>() : which includes CoordinateVector

    template <int n, typename T>
    bool get_value(Vector<n, T> &val, const std::string &label,
                   bool bcast = true) {
        val = 0;
        bool no_error = true;

        if (hila::myrank() == 0) {
            no_error = get_value(val[0], label, false);
            for (int i = 1; i < n && no_error; i++) {
                no_error = get_value(val[i], ",", false);
            }

            if (!no_error && speaking) {
                hila::output << "Error: expecting " << n << " comma-separated "
                             << type_id<T>() << "s after '" << label << "'\n";
            }
        }

        if (bcast) {
            broadcast(val, no_error);
        }
        return no_error;
    }

    /// Specialization to CoordinateVector

    template <int n = NDIM>
    bool get_value(CoordinateVector &val, const std::string &label,
                   bool bcast = true) {
        Vector<n, int> iv;
        bool b = get_value(iv, label, bcast);
        val = iv;
        return b;
    }

    /// Specialize -get_value() to std::vector<>

    template <typename T>
    bool get_value(std::vector<T> &val, const std::string &label,
                   bool bcast = true) {
        val = {};
        bool no_error = true;

        if (hila::myrank() == 0) {
            T v;
            no_error = get_value(v, label, false);
            val.push_back(v);
            while (no_error && match_token(",")) {
                no_error = get_value(v, "", false);
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

    int get_item(const std::string &label,
                 const std::vector<std::string> &items, bool bcast = true);

  private:
    /// a helper method to give type name
    template <typename T> inline const char *type_id() { return nullptr; }

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

template <> inline const char *input::type_id<int>() { return "int"; }
template <> inline const char *input::type_id<long>() { return "long"; }
template <> inline const char *input::type_id<float>() { return "float"; }
template <> inline const char *input::type_id<double>() { return "double"; }
template <> inline const char *input::type_id<std::string>() {
    return "string";
}
template <> inline const char *input::type_id<Cmplx<float>>() {
    return "complex value";
}
template <> inline const char *input::type_id<Cmplx<double>>() {
    return "complex value";
}

#endif
