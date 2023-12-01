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

namespace hila
{
// Multicanonical methods are separated to their own namespace
namespace muca
{
// Function pointer to the iteration function
typedef bool (* iteration_pointer)(const double OP);
extern iteration_pointer iterate_weights;

// Quick helper function for writing values to a file
template <class K>
void to_file(std::ofstream &output_file, std::string fmt, K input_value);

// Generates timestamped file names
std::string generate_outfile_name();

// Reads parameters for muca computations
void read_weight_parameters(std::string parameter_file_name);

// Reads weight functions from a file
void read_weight_function(std::string W_function_filename);

// Writes weight functions to a file
void write_weight_function(std::string filename);

// Gives the weight as a function of the order parameter
double weight_function(double OP);
double weight(double OP);

// Accept/reject determination for pairs of order parameter values
bool accept_reject(const double OP_old, const double OP_new);

// Set the direct iteration finish condition
void set_direct_iteration_FC(bool (* fc_pointer)(std::vector<int> &n));

// Set to perform the weight iteration at each call to accept_reject
void set_continuous_iteration(bool YN);

// For the continuous iteration the finish condition is tracked internally
// and can be checked and set using the two functions below
bool check_weight_iter_flag();
void set_weight_iter_flag(bool YN);

// Initialises the muca computations according to the weight parameter file.
// This function is to always be called before using any of the above functions
void initialise(const std::string wfile_name);

////////////////////////////////////////////////////////////////////////////////
// Static functions are internal to above methods. See .cpp file for details
////////////////////////////////////////////////////////////////////////////////

static void bin_OP_value(double OP);

static int find_OP_bin_index(double OP);

static bool all_visited(std::vector<int> &n);
static bool first_last_visited(std::vector<int> &n);

static void overcorrect(std::vector<double> &Weight, std::vector<int> &n_sum);

static void recursive_weight_iteration(std::vector<double> &Weight,
                                       std::vector<int> &n,
                                       std::vector<int> &g_sum,
                                       std::vector<double> &g_log_h_sum);

static void print_iteration_histogram();

static std::vector<double> get_equidistant_bin_limits();

static void setup_equidistant_bins();

static void initialise_weight_vectors();

static bool iterate_weight_function_direct(double OP);

static bool iterate_weight_function_direct_single(double OP);

static void setup_iteration();

}
}

#endif
