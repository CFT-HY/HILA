/** @file input.h */
#ifndef PARAM_INPUT_H
#define PARAM_INPUT_H

#include <string>
#include "defs.h"
#include "lattice.h"
#include "plumbing/com_mpi.h"

////////////////////////////////////////////////////////////////////////
/// hila::input - Class for parsing runtime parameter files
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
///    hila::input f("filename");   - initialization with filename opens
///                                   the file for input
///    hila::input f;
///----------------------------------------------------------------------------
///
/// get_item(): select one item from a "menu":
///
///        int hila::input::get_item(std::string key,
///                        std::vector<std::string> items,
///                        bool broadcast = true);
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
///----------------------------------------------------------------------------
///
/// get_value(): read input (alternative to get())
///
///        template <typename T>
///        bool hila::input::get_value(T & val,std::string key,
///                                    bool broadcast=true);
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
///----------------------------------------------------------------------------
///
/// close():    hila::input::close()
///
///         f.close();
///    closes the input f.  Now "f.open("file")" can be used again.
///    File is also closed when variable "f" goes out of scope.
///
///----------------------------------------------------------------------------
///
/// quiet():    hila::input::quiet(bool be_silent=true)
///
///         f.quiet();      // don't print read items to hila::out
///         f.quiet(false); // re-enable printing
///
///    By default hila::input methods print everything read to hila::out
///    for logging.  f.quiet(); disables this.
///
///------------------------------------------------------------------
///
/// NOTE: methods which broadcast to all nodes (default) must be called
///       from all nodes synchronously. These include open(), get(),
///       get_value() with bcast=true, get_item with bcast=true.
///
///       Thus;
///          if (hila::myrank() == 0) {
///              double v = f.get("a value");
///              ...
///          }
///
///       deadlocks (if there are more than 1 rank). Method
///          f.get_value(v,"a value",false);
///       can be used in this context.
///
////////////////////////////////////////////////////////////////////////

namespace hila {

/**
 * @brief hila::input - Class for parsing runtime parameter files.
 *
 * @details Input files consist normally of "key <value>" -pairs.
 *
 * It is important to note that the data structure employed by the input class is not a dictionary,
 * but technically a stack. This means that the values are read in sequentially, thus they must be
 * in the file in the order they are read in.
 *
 * \code{.cpp}
 *  hila::input f("filename");   // initialization with filename opens the file for input
 * \endcode
 *
 * Comment character '#': everything after # in input file is a comment to the end of the line.
 */
class input {

  private:
    std::ifstream inputfile;
    bool is_initialized = false;
    std::string filename;
    bool use_cin;

    std::string linebuffer;
    size_t lb_start = 0; // linebuffer start index
    bool is_line_printed;
    bool speaking = true; // false does not print information

  public:
    input() {}
    ~input() {
        close();
    }
    input(const std::string &fname) {
        open(fname);
    }

    /**
     * @brief Open file that parameters are read from
     * @details If input class is initialized with path to file that needs to be opened,
     * then the file will be opened automatically, and this method does not need to be called.
     *
     * If no argument is given then filename will be interperated as default input file name
     *
     * In the case that use_cmdline is True and if no argument is given and filename is given with
     * -i {alternative_file} cmdline argument then {alternative_file} will be the specified file to
     * be opened.
     *
     * @param fname Path to file that needs to be read
     * @param use_cmdline Default True
     * @param exit_on_error If true exit on error. Return value is passed to all MPI nodes. Default
     * True
     * @return true
     * @return false
     */
    bool open(const std::string &fname, bool use_cmdline = true, bool exit_on_error = true);
    void close();

    // make class quiet (no printouts), quiet(false) returns to normal
    void quiet(bool really = true) {
        speaking = !really;
    }

    /// returntyhpe is a special class for resolving get("label") return type
    class returntype {
      public:
        std::string label;
        input *parent;

        returntype(const std::string &str, input *a) : label(str), parent(a) {}

        /// cast operator does the conversion - disable nullptr_t cast used
        /// in hila types
        /// disable type "char" because then c++ does not know how to initialize strings
        template <typename T, std::enable_if_t<(hila::is_complex_or_arithmetic<T>::value &&
                                                !std::is_same<T, char>::value) ||
                                                   std::is_same<T, std::string>::value,
                                               int> = 0>
        operator T() {
            T val;
            if (!parent->get_value(val, label, true))
                hila::finishrun();
            return val;
        }

        template <typename T, int n>
        operator Vector<n, T>() {
            Vector<n, T> val;
            if (!parent->get_value(val, label, true))
                hila::finishrun();
            return val;
        }

        template <typename T>
        operator CoordinateVector_t<T>() {
            CoordinateVector_t<T> val;
            if (!parent->get_value(val, label, true))
                hila::finishrun();
            return val;
        }

        template <typename T, std::enable_if_t<hila::is_complex_or_arithmetic<T>::value ||
                                                   std::is_same<T, std::string>::value,
                                               int> = 0>
        operator std::vector<T>() {
            std::vector<T> val;
            if (!parent->get_value(val, label, true))
                hila::finishrun();
            return val;
        }
    };

    /**
     * @brief Get next value in stack of read in input string from parameters file.
     * @details Use as
     * \code{.cpp}
     *    var = f.get("key");
     * \endcode
     * reads in a key-value pair
     * \code{.txt}
     *         key   <value(s)>
     * \endcode
     * from the input file f, and returns the value of type of variable var.
     * The value is broadcast to all MPI nodes.  The method infers
     * the type of the returned variable from the type of the assignment.
     *
     * Key is an alphanumeric string, which may contain words separated by whitespace.
     *
     * The get method simply traverses the parsed in stack of input values in the given
     * parameters file. A functionality of the method is that the given argument key is skipped and
     * the value that this key is assigned to is then fetched and returned.
     *
     * Technically the key's can also be read in if no argument is given to get. The main
     * functionality of the key argument is to book keep by the user that all parameters are read
     * in. If the key that is wanting to be read doesn't exist, then get throws an error and the
     * program is stopped.
     *
     * With this logic one could construct a parameters file of only values:
     *
     * \code {.txt}
     * value_1
     * value_2
     * .
     * .
     * .
     * value_n
     * \endcode
     *
     * and read in the values in a loop by simply calling hila::input::get() consecutively, but this
     * is not advised.
     *
     * The type of the value to be read in is inferred by the variable the value is fetched into. If
     * the value cannot be casted into the variable, then the parser will throw an error and crash.
     *
     * ## Types which the parser supports with examples
     *
     * Multiple items are separated by commas, whitespace is not significant.
     *
     * ### Any arithmetic type (ints/floats)
     * @code{.cpp}
     *        int i = f.get("number of cats");
     * @endcode
     *    __matches__ " number of cats   5 "
     *
     * ### Complex<float/double>
     * @code{.cpp}
     *        Complex<double> phase = f.get("complex phase");
     * @endcode
     *    __matches__ "complex phase   (0.4, 0.5)"
     *
     *    complex values are given in pairs within ( , )
     *
     * ### std::string
     * @code{.cpp}
     *        std::string s = f.get("key");
     * @endcode
     *    __matches__ "key <string value>" where string value is either
     *
     *     - sequence of non-whitespace chars, delimited by whitespace, eol, ','
     *        or '#'.
     *
     *     - characters enclosed by quotes "..".  These have to pair
     *        within the same line.  Quote marks are removed.
     *
     * ### CoordinateVector
     * @code{.cpp}
     *        CoordinateVector v;
     *        v = f.get("lattice size");
     * @endcode
     *    __matches__ "lattice  size  32, 32, 32, 32"  (if NDIM == 4).
     *
     * ### Vector<T,int>
     * @code{.cpp}
     * Vector<5,int> vec = f.get("initial vector")
     * @endcode
     *    __matches__ "initial vector 1, 2, 3, 4, 5"
     *
     * ### std::vector<T>
     * @code{.cpp}
     *        std::vector<double> dvec = f.get("vec");
     * @endcode
     *    __matches__ "vec  3,4, 5.5, 7.8, 4"  and returns a vector of double values.
     *    The numbers are read until they are not followed by a comma.  If comma
     *    is the last non-whitespace character on the line, reading continues to
     *    the next line.
     *
     *  T is one of the other supported types.
     *
     * __NOTE__: The parser will crash if during reading in a multi valued line (Vector,
     * std::vector) the line ends with a ",". The parser will interpret that there is data on the
     * next line to read in and will fail during the next get call.
     *
     * @param key Parameter to be read in.
     * @return returntype Value corresponding to input key
     */
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
    bool get_value(T &val, const std::string &label, bool bcast = true) {
        val = {};
        bool no_error = handle_key(label); // removes whitespace

        if (hila::myrank() == 0 && no_error) {
            std::string tok;
            if (!(get_token(tok) && is_value(tok, val))) {

                if (speaking)
                    hila::out << "Error: expecting a value of type '" << type_id<T>() << "' after '"
                              << label << "'\n";

                no_error = false;
            }
        }

        if (bcast) {
            // string has to be treated separately
            if constexpr (std::is_same<T, std::string>::value) {
                hila::broadcast(val);
                hila::broadcast(no_error);
            } else {
                hila::broadcast2(val, no_error);
            }
        }

        return no_error;
    }

    /// Specialize the above method to Complex -pair:  (re,im)

    template <typename T>
    bool get_value(Complex<T> &val, const std::string &label, bool bcast = true) {
        val = 0;
        bool no_error = handle_key(label); // removes whitespace

        if (hila::myrank() == 0 && no_error) {
            std::string tok;
            T re, im;

            no_error =
                (match_token("(") && get_token(tok) && is_value(tok, re) && match_token(",") &&
                 get_token(tok) && is_value(tok, im) && match_token(")"));
            if (!no_error && speaking) {
                hila::out << "Error: expecting complex value '(re,im)' after '" << label << "'\n";
            }

            val = Complex<T>(re, im);
        }

        if (bcast) {
            hila::broadcast2(val, no_error);
        }
        return no_error;
    }

    /// Specialize .get_value<Vector<n,T>>() : which includes CoordinateVector

    template <int n, typename T>
    bool get_value(Vector<n, T> &val, const std::string &label, bool bcast = true) {
        val = 0;
        bool no_error = true;

        if (hila::myrank() == 0) {
            no_error = get_value(val[0], label, false);
            for (int i = 1; i < n && no_error; i++) {
                no_error = get_value(val[i], ",", false);
            }

            if (!no_error && speaking) {
                hila::out << "Error: expecting " << n << " comma-separated " << type_id<T>()
                          << "s after '" << label << "'\n";
            }
        }

        if (bcast) {
            hila::broadcast2(val, no_error);
        }
        return no_error;
    }

    /// Specialization to CoordinateVector

    template <int n = NDIM>
    bool get_value(CoordinateVector &val, const std::string &label, bool bcast = true) {
        Vector<n, int> iv;
        bool b = get_value(iv, label, bcast);
        val = iv;
        return b;
    }

    /// Specialize -get_value() to std::vector<>

    template <typename T>
    bool get_value(std::vector<T> &val, const std::string &label, bool bcast = true) {
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
                hila::out << "Error: expecting a comma-separated list of " << type_id<T>()
                          << "s after '" << label << "'\n";
            }
        }

        if (bcast) {
            hila::broadcast(val);
            hila::broadcast(no_error);
        }

        return no_error;
    }


    /**
     * @brief Get the item object
     *
     * @param label
     * @param items
     * @param bcast
     * @return int
     */
    int get_item(const std::string &label, const std::vector<std::string> &items,
                 bool bcast = true);

  private:
    /// a helper method to give type name
    template <typename T>
    inline const char *type_id() {
        return nullptr;
    }

    bool peek_token(std::string &tok);
    bool get_token(std::string &tok);
    bool match_token(const std::string &tok);

    bool scan_string(std::string &val);

    template <typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    bool is_value(const std::string &s, T &val) {
        std::istringstream ss(s);
        ss >> val;
        if (ss.fail() || ss.bad())
            return false;
        else
            return true;
    }

    bool is_value(const std::string &s, std::string &val);

    bool contains_word_list(const std::string &list, int &end_of_key);

    std::string remove_quotes(const std::string &val);

    void print_linebuf(int eok);

    bool get_line();
    bool handle_key(const std::string &c);

    bool remove_whitespace();
};

/// give here specializations of the type_id helper

template <>
inline const char *input::type_id<int>() {
    return "int";
}
template <>
inline const char *input::type_id<long>() {
    return "long";
}
template <>
inline const char *input::type_id<long long>() {
    return "long long";
}
template <>
inline const char *input::type_id<unsigned int>() {
    return "unsigned int";
}
template <>
inline const char *input::type_id<unsigned long>() {
    return "unsigned long";
}
template <>
inline const char *input::type_id<unsigned long long>() {
    return "unsigned long long";
}

template <>
inline const char *input::type_id<float>() {
    return "float";
}
template <>
inline const char *input::type_id<double>() {
    return "double";
}
template <>
inline const char *input::type_id<std::string>() {
    return "string";
}
template <>
inline const char *input::type_id<Complex<float>>() {
    return "complex value";
}
template <>
inline const char *input::type_id<Complex<double>>() {
    return "complex value";
}


} // namespace hila

#endif
