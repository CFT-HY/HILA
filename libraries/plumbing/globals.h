#ifndef HILA_GLOBAL_H_
#define HILA_GLOBAL_H_
/**
 * @file globals.h
 * @brief Definition of global variable class
 */

namespace hila {
/**
 * @brief Global variable class within hila namespace
 *
 * @details Special wrapper for "simple" types (elementary types, simple structs etc.).
 * Implemented on gpus automatically using `__constant__` memory. All global variables
 * exist for the whole program lifetime.
 *
 * ## Variable declaration:
 *
 * Declaration is possible on top (file) level only, not inside functions!
 * \code{.cpp}
 * hila::global<type> globalvar;
 * \endcode
 * For example
 * \code{.cpp}
 * hila::global<Vector<4,Complex<double>> cvec1, cvec2;
 * \endcode
 *
 * ## Variable assingment:
 * \code{.cpp}
 * globalvar.set(<value>);
 * \endcode
 *
 * Assignment is not possible inside loops
 * Only "full" assignment is possible, not e.g. by struct records.  Value must be convertible
 * to the type of the global
 *
 * ## Variable access
 *
 * The value of the variable is obtained with .get() -method
 * \code{.cpp}
 * globalvar.get();
 * \endcode
 * .get() returns a const reference to the variable, it does no copying
 * Variable can be used everywhere in the source file.
 *
 * \code{.cpp}
 * // at global file level
 * struct param_t { double a,b[2]; };
 *
 * hila::global<param_t> params;
 * params.set( {1,2,3} );
 *
 * // define some function
 * double myfunc(double dv) {
 *     return cos( dv / params.get().a );
 * }
 *
 * ...
 * // inside main() or some other function - use helper struct to assign values
 * param_t tmp;
 * tmp.a = 3*M_PI;
 * tmp.b[0] = ....
 * tmp.b[1] = ...
 * params.set(tmp);  // do the full struct assignment to global
 *
 * ...
 * Field<double> df;
 * onsites(ALL) {
 *     df[X] = myfunc(X.x()) + params.get().b[0]/params.get().b[1];
 * }
 * hila::out0 << "Parameter a is " << params().a << '\n';
 * \endcode
 *
 * ## Additional information
 *
 * "extern" and "static" can be used, if there are several source files, with the usual meaning:
 *
 * \code{.cpp}
 * extern hila::global<double> a;      // variable a (of identical type) is defined
 *                                     // in some other file
 * static hila::global<mytype> par;    // par is invisible to other files
 * \endcode
 *
 * Declarations can be enclosed in namespaces:
 * \code{.cpp}
 * namespace myglobals {
 *   hila::global<double> a, b;
 * }
 * ...
 *
 * myglobals::a.set(3);
 * hila::out0 << "a has value " << myglobals::a.get() << '\n';
 * \endcode
 *
 * Variables cannot be initialized in declaration, because (possible) GPUs are not initialized.
 */
template <typename T, typename custom = void>
class global {

    static_assert(std::is_trivial<T>::value && std::is_standard_layout<T>::value,
                  "hila::global<> expects only pod-type elements (plain old data): default "
                  "constructor, copy and delete");

  private:
    T val;

    // empty stub, will be specialized by hilapp (if needed)
    void copy_to_device() const {}

  public:
    // Get the value of the variable with operator()

    const T &get() const {
        return val;
    }

    const T &operator()() const {
        return val;
    }

    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void set(const S &rhs) {
        assert(hila::is_initialized && "Assign to global possible only after hila::initialize()");

        val = rhs;

        copy_to_device();
    }


    // assignment is possible from compatible rhs values
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void operator=(const S &rhs) {
        this->set(rhs);
    }
};
} // namespace hila

#endif
