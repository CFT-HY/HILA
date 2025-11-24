////////////////////////////////////////////////////////////////////////////////
/// @file multicanonical.h
/// @author Jaakko Hällfors
/// @brief Header for model agnostic implementation of
///        various multicanonical (muca) methods
////////////////////////////////////////////////////////////////////////////////
#ifndef MULTICANONICAL_HEADER
#define MULTICANONICAL_HEADER

// !OBS!
// Please see the doxygen documentation from github for further details. The
// corresponding .cpp file contains the actual implementations and also
// documents the internal functions (not included in doxygen) in relative
// detail.
// Furthermore, the repository will contain a short example program that
// elucidates the use of the various muca methods.

// NOTE: TODO: at the moment different iteration methods `direct`, `direct_smooth`, `canonical` will
// result into differently normalized weight functions. End results should be the same upto a some
// constant shift. This does not matter in the end since only differences matter.

#include <vector>
#include <string>

namespace hila {
// Multicanonical methods are separated to their own namespace
// namespace muca {

struct Muca;

// Quick helper function for writing values to a file
template <class... Ks>
static inline void to_file(std::ofstream &output_file, std::string fmt, Ks... input_value);

/// For generating time stamped file names
std::string append_time_stamp(const std::string &outfile_name_base);

////////////////////////////////////////////////////////////////////////////////
/// @struct weight_iteration_parameters
/// @brief An internal struct parametrising multicanonical methods.
///
/// @var string weight_iteration_parameters::weight_loc
/// @brief Path to the weight function file
///
/// @var string weight_iteration_parameters::outfile_name_base
/// @brief Prefix for the saved weight function files
/// @details
///
/// @var string weight_iteration_parameters::method
/// @brief Name of the iteration method
/// @details Current options: "direct"
/// See the documentation for the details of different methods.
///
/// @var bool weight_iteration_parameters::visuals
/// @brief Whether to print histogram during iteration
///
/// @var bool weight_iteration_parameters::hard_walls
/// @brief Whether the weight outside max_OP and min_OP is infinite
/// @details If not, the weight is assigned through constant extrapolation of
/// the nearest bin (first or last). Please start your simulations so that the
/// first configuration is within the interval, or the iteration can get stuck.
///
/// @var double weight_iteration_parameters::max_OP
/// @brief Maximum order parameter value
/// @details Only matters when creating bins from given parameters. When the
/// details of the weights are read from file this parameter is not used.
///
/// @var double weight_iteration_parameters::min_OP
/// @brief Minimum order parameter value
/// @details Only matters when creating bins from given parameters. When the
/// details of the weights are read from file this parameter is not used.
///
/// @var int weight_iteration_parameters::bin_number
/// @brief Number of OP bins
/// @details Only matters when creating bins from given parameters. When the
/// details of the weights are read from file this parameter is not used.
///
/// @var bool weight_iteration_parameters::AR_iteration
/// @brief Whether to update the weights after each call to accept_reject
///
/// @var struct direct_iteration weight_iteration_parameters::DIP
/// @brief A substruct that contains method specific parameters
///
/// @var struct canonical_iteration weight_iteration_parameters::CIP
/// @brief A substruct that contains method specific parameters
///////////////////////////////////////////////////////////////////////////////
struct weight_iteration_parameters {
    std::string weight_loc;
    std::string outfile_name_base;
    std::string method;

    bool visuals;
    bool hard_walls;
    double min_OP;
    double max_OP;
    int bin_number;
    bool AR_iteration;
};

/// type of the finish condition function
typedef bool (*finish_condition_fn)(const Muca &muca);
// Different finish condition funcs
bool all_visited(const Muca &muca);
bool first_last_visited(const Muca &muca);

/// type of the iteration function
typedef bool (*iteration_fn)(Muca &muca, const double OP);
// Different iteration funcs
// iteration fuctions
bool iterate_weight_function_direct(Muca &muca, double OP);
bool iterate_weight_function_direct_single(Muca &muca, double OP);
bool iterate_weight_function_direct_smooth(Muca &muca, double OP);
bool iterate_weight_function_canonical(Muca &muca, double OP);

/// Muca struct
struct Muca {
    /////////////////////////////////////////////////////////////
    // Intended interface:

    /// Initialises the muca computations according to the weight parameter file `wfile_name`.
    /// This function is to always be called before using any of the below functions.
    /// (see, Example parameter file in applications/multicanonical_example/muca_parameters.)
    bool initialise(const std::string wfile_name);
    /// Pointer to the iteration function, set by initialisation (or manually if custom)
    iteration_fn iterate_weights;
    /// Writes weight functions to a file
    bool write_weight_function(const std::string &W_function_filename);
    /// Gives the weight as a function of the order parameter
    double weight_function(double OP) const;
    double weight(double OP) const;
    /// Accept/reject determination for pairs of order parameter values
    bool accept_reject(const double OP_old, const double OP_new);


    /////////////////////////////////////////////////////////////
    // Internals (in normal use should be considered private)
    // private:

    // Data: ----------------------------------------------------

    /// parameter struct filled by read_weight_parameters
    weight_iteration_parameters WParam;

    /// Weight function is parametrized as piecewise linear func
    ///
    /// The order parameter (OP) range [min_OP, max_OP] is divided into bins  \n<pre>
    ///    'c' `OP_bin_centers`:    0   1   2   3    ... N                    \n
    ///                           | c | c | c | c |  ...                      \n
    ///     '|' `OP_bin_limits`:  0   1   2   3   4  ... N + 1                \n</pre>
    /// Thus limits for OP_bin_centers[i] are OP_bin_limits[i] to OP_bin_limits[i + 1]
    /// The weights are tracked on the centers of each bin ,stored in `weights`
    std::vector<double> weights;        ///< The weight values at each OP bin center
    std::vector<double> OP_bin_centers; ///< Positions of the bin centers in OP
    std::vector<double> OP_bin_limits;  ///< Limits of the bins in OP

    int weightIterationCount = 0;
    bool weightIterationFlag = true;
    finish_condition_fn finish_check; ///< Pointer to the finish condition check

    std::vector<int> OP_bin_hits;       ///< Number of hits in OP bins
    std::vector<int> OP_bin_hits_total; ///< Total accumulated hits in OP bins

    /// An internal struct for direct weight iteration methods parameters and iteration data.
    struct direct_iteration {
        /// @brief Determines the iteration condition for the iteration method.
        /// @details Current options: "all_visited", "ends_visited". The weight
        /// modification factor \f$C\f$ is decreased after set bins are visited. These
        /// can include e.g. all bins, or just the ends (the two preset options).
        /// Finishing conditions can also be determined by the user and passed to the
        /// methods through muca::direct_iteration_finish_condition(&func_pointer)
        std::string finish_condition;
        /// Number of samples before weight function update
        int sample_size;
        /// Determines how often the update condition is checked for the case
        /// direct_iteration::sample_size = 1
        int single_check_interval;
        /// Initial magnitude of the weight update
        double C_init;
        /// Minimum magnitude of the weight update
        double C_min;
        /// @brief Magnitude of the weight update
        /// @details This entry is decreased through the direct iteration method until
        /// \f$ C < C_\mathrm{min}\f$ at which point the iteration is considered
        /// complete. The magnitude of the update is such that the mean weight
        /// modification is \f$C\f$. That is, the update satisfies
        /// \f{equation}{\sum_i \delta W_i = N C,\f}
        /// where \f$N\f$. is the number of bins. For sample_size = 1 we simply add
        /// \f$ C \f$ to the single relevant bin each iteration.
        double C;
    } DI;

    /// An internal struct for canonical weight iteration methods parameters and iteration data.
    struct canonical_iteration {
        /// Number of samples before weight function update
        int sample_size;
        /// Initial number of entries in the bin totals. Larger number makes the weight less
        /// susceptible to changes in the following iterations.
        int initial_bin_hits;
        /// Minimum number of hits a bin needs to have before taken into account in the weight
        /// function update.
        int min_bin_hits;
        /// stop iteration after reaching this many updates of the weight func
        int max_iters;
        ///
        double OC_factor;
        /// Up to which iteration step the overcorrection updates can be used
        int OC_max_iter; // TODO: not used atm
        /// How often the overcorrection update is used. If the value is n, every n:th update will
        /// be an overcorrection.
        int OC_frequency; // TODO: not used atm

        int weight_update_count;         ///< How many times the weight function has been updated
        std::vector<double> can_hist;    ///< estimate of the canonical histogram `h_i^k`
        std::vector<int> hits_nsum;      ///< cumulative sum of the number of hits in bins `n_i^k`
        std::vector<int> gsum;           ///< cumulative sum of the two-bin weight factor `g_i^k`
        std::vector<double> noc_weights; ///< buffer to save the Non-OverCorrected weight function
    } CI;


    // Internal functions: --------------------------------------

    void read_weight_parameters(std::string parameter_file_name);
    // Reads weight functions from a file
    bool read_weight_function(const std::string &W_function_filename);
    void initialise_weight_vectors();
    void setup_iteration();

    // Direct iteration stuff
    /// Set the direct iteration finish condition
    void set_direct_iteration_FC(finish_condition_fn fc);
    void bin_hit_OP_value(double OP);
    void setup_equidistant_bins();
    void print_iteration_histogram() const;

    /// Set to perform the weight iteration at each call to accept_reject
    void set_continuous_iteration(bool YN);
    /// For the continuous iteration the finish condition is tracked internally
    /// and can be checked and set using the two functions below
    bool check_weight_iter_flag();
    void set_weight_iter_flag(bool YN);

    void overcorrect(std::vector<double> &Weight, std::vector<int> &n_sum);
    void recursive_weight_iteration(std::vector<double> &Weight, std::vector<int> &n,
                                    std::vector<int> &g_sum, std::vector<double> &g_log_h_sum);
};

// static helper funcs
static inline std::vector<double> get_equidistant_bin_limits(double min, double max, int N_bins);
static inline int find_OP_bin_index(double OP, const std::vector<double> &OP_bin_limits);

//} // namespace muca
} // namespace hila

#endif
