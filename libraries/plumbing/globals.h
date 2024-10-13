#ifndef HILA_GLOBAL_H_
#define HILA_GLOBAL_H_

/**
 * hila::global<> variable class.
 * Special wrapper for "simple" types (elementary types, simple structs etc.).
 * Implemented on gpus automatically using __constant__ memory. All global variables
 * exist for the whole program lifetime.
 *
 * Variable declaration:
 * Declaration is possible on top (file) level only, not inside functions!
 *   > hila::global< type > globalvar;
 * For example  hila::global<Vector<4,Complex<double>> cvec1, cvec2;
 *
 * Variable assingment:
 *   > globalvar = <value>;
 * Assignment is not possible inside loops
 * Only "full" assignment is possible, not e.g. by field records
 *
 * The value of the variable is obtained with function-like "()"
 *   > globalvar()
 * Variable can be used everywhere in the source file.
 *
 * Example:
 *          // at global file level
 *          struct param_t { double a,b[2]; };
 *
 *          hila::global<param_t> params;
 *
 *          // define some function
 *          double myfunc(double dv) {
 *              return cos( dv / params().a );
 *          }
 *
 *          ...
 *          // inside main() or some other function - use helper struct to assign values
 *          param_t tmp;
 *          tmp.a = 3*M_PI;
 *          tmp.b[0] = ....
 *          tmp.b[1] = ...
 *          params = tmp;  // do the full struct assignment to global
 *
 *          ...
 *          Field<double> df;
 *          onsites(ALL) {
 *              df[X] = myfunc(X.x()) + params().b[0]/params().b[1];
 *          }
 *          hila::out0 << "Parameter a is " << params().a << '\n';
 *
 * "extern" and "static" can be used, if there are several source files, with the usual meaning:
 *          extern hila::global<double> a;      // variable a (of identical type) is defined
 *                                              // in some other file
 *          static hila::global<mytype> par;    // par is invisible to other files
 *
 * Declarations can be enclosed in namespaces:
 *          namespace myglobals {
 *          hila::global< double > a, b;
 *          }
 *          ...
 *
 *          myglobals::a = 3;
 *          hila::out0 << "a has value " << myglobals::a() << '\n';
 *
 * Variables cannot be initialized in declaration, because (possible) GPUs are not initialized.
 * In
 */


namespace hila {
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
    const T &operator()() const {
        return val;
    }

    // assignment is possible from compatible rhs values
    template <typename S, std::enable_if_t<hila::is_assignable<T &, S>::value, int> = 0>
    void operator=(const S &rhs) {
        assert(hila::is_initialized && "Assign to global possible only after hila::initialize()");

        val = rhs;

        copy_to_device();
    }
};
} // namespace hila

#endif
