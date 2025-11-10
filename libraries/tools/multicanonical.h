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

#include <vector>
#include <string>

namespace hila {
// Multicanonical methods are separated to their own namespace
namespace muca {

struct Muca;

// Quick helper function for writing values to a file
template <class... Ks>
static inline void to_file(std::ofstream &output_file, std::string fmt, Ks... input_value);

// Generates timestamped file names
std::string generate_outfile_name(const std::string &outfile_name_base);

// Structs to encompass the various options related to
// iteration of the multicanonical weights.
////////////////////////////////////////////////////////////////////////////////
/// @struct canonical_iteration
/// @brief An internal struct parametrising the canonical weight
/// iteration method.
///
/// @var int canonical_iteration::sample_size
/// @brief Number of samples before weight function update
///
/// @var int canonical_iteration::initial_bin_hits
/// @brief Initial number of entries in the bin totals
/// @details Larger number makes the weight less susceptible to changes in the
/// following iterations.
///
/// @var int canonical_iteration::OC_max_iter
/// @brief Up to which iteration the overcorrection updates can be used
///
/// @var int canonical_iteration::OC_frequency
/// @brief How often the overcorrection update is used
/// @details If the value is n, every n:th update will be an overcorrection.
///////////////////////////////////////////////////////////////////////////////
struct canonical_iteration {
    int sample_size;
    int initial_bin_hits;
    int OC_max_iter;
    int OC_frequency;
};

////////////////////////////////////////////////////////////////////////////////
/// @struct direct_iteration
/// @brief An internal struct parametrising the direct weight iteration method.
///
/// @var string direct_iteration::finish_condition
/// @brief Determines the iteration condition for the iteration method.
/// @details Current options: "all_visited", "ends_visited". The weight
/// modification factor \f$C\f$ is decreased after set bins are visited. These
/// can include e.g. all bins, or just the ends (the two preset options).
/// Finishing conditions can also be determined by the user and passed to the
/// methods through muca::direct_iteration_finish_condition(&func_pointer)
///
/// @var int direct_iteration::sample_size
/// @brief Number of samples before weight function update
///
/// @var int direct_iteration::single_check_interval
/// @brief Determines how often the update condition is checked for the case
/// direct_iteration::sample_size = 1
///
/// @var double direct_iteration::C_init
/// @brief Initial magnitude of the weight update
///
/// @var double direct_iteration::C_min
/// @brief Minimum magnitude of the weight update
///
/// @var double direct_iteration::C
/// @brief Magnitude of the weight update
/// @details This entry is decreased through the direct iteration method until
/// \f$ C < C_\mathrm{min}\f$ at which point the iteration is considered
/// complete. The magnitude of the update is such that the mean weight
/// modification is \f$C\f$. That is, the update satisfies
/// \f{equation}{\sum_i \delta W_i = N C,\f}
/// where \f$N\f$. is the number of bins. For sample_size = 1 we simply add
/// \f$ C \f$ to the single relevant bin each iteration.
///////////////////////////////////////////////////////////////////////////////
struct direct_iteration {
    std::string finish_condition;
    int sample_size;
    int single_check_interval;
    double C_init;
    double C_min;
    double C;
};

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
    double max_OP;
    double min_OP;
    int bin_number;
    bool AR_iteration;
    struct direct_iteration DIP;
    struct canonical_iteration CIP;
};

// type of the finish condition function
typedef bool (*finish_condition_fn)(const Muca &muca);
// Different finish condition funcs
bool all_visited(const Muca &muca);
bool first_last_visited(const Muca &muca);

// type of the iteration function
typedef bool (*iteration_fn)(Muca &muca, const double OP);
// Different iteration funcs
// iteration fuctions
bool iterate_weight_function_direct(Muca &muca, double OP);
bool iterate_weight_function_direct_single(Muca &muca, double OP);
bool iterate_weight_function_direct_smooth(Muca &muca, double OP);

struct Muca {
    // parameter struct filled by read_weight_parameters
    weight_iteration_parameters WParam;

    std::vector<double> OPValues;
    std::vector<double> OPBinLimits;
    std::vector<double> WValues;

    int weightIterationCount = 0;
    bool weightIterationFlag = true;
    // Pointer to the iteration function
    iteration_fn iterate_weights;
    // Pointer to the finish condition check
    finish_condition_fn finish_check;

    // Direct iteration
    std::vector<int> N_OP_Bin;
    std::vector<int> N_OP_BinTotal;

    // Initialises the muca computations according to the weight parameter file.
    // This function is to always be called before using any of the above functions
    bool initialise(const std::string wfile_name);

    /////////////////////////////////////////////////////////////
    // Reads parameters for muca computations
    void read_weight_parameters(std::string parameter_file_name);

    // Reads weight functions from a file
    bool read_weight_function(const std::string &W_function_filename);

    // Writes weight functions to a file
    bool write_weight_function(const std::string &W_function_filename);

    // Gives the weight as a function of the order parameter
    double weight_function(double OP) const;
    double weight(double OP) const;

    // Accept/reject determination for pairs of order parameter values
    bool accept_reject(const double OP_old, const double OP_new);

    // Set the direct iteration finish condition
    void set_direct_iteration_FC(finish_condition_fn fc);

    // Set to perform the weight iteration at each call to accept_reject
    void set_continuous_iteration(bool YN);

    // For the continuous iteration the finish condition is tracked internally
    // and can be checked and set using the two functions below
    bool check_weight_iter_flag();
    void set_weight_iter_flag(bool YN);

    /////////////////////////////////////////////////////////////
    // private:
    void bin_OP_value(double OP);

    void overcorrect(std::vector<double> &Weight, std::vector<int> &n_sum);
    void recursive_weight_iteration(std::vector<double> &Weight, std::vector<int> &n,
                                    std::vector<int> &g_sum, std::vector<double> &g_log_h_sum);

    void print_iteration_histogram() const;

    void setup_equidistant_bins();
    void initialise_weight_vectors();

    void setup_iteration();
};


// static helper funcs
static inline std::vector<double> get_equidistant_bin_limits(double min, double max, int N_bins);
static inline int find_OP_bin_index(double OP, const std::vector<double> &OPBinLimits);

} // namespace muca
} // namespace hila

#endif
