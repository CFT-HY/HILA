////////////////////////////////////////////////////////////////////////////////
/// @file multicanonical.cpp
/// @author Jaakko Hällfors
/// @brief Model agnostic implementation of various multicanonical methods
/// @details TBA
////////////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <regex>
#include <cmath>
#include "hila.h"
#include "tools/multicanonical.h"

namespace hila {
namespace muca {

// File global parameter struct filled by read_weight_parameters
static weight_iteration_parameters g_WParam;

// Initialise some static vectors for this file only.
static vector g_OPValues(1, 0);
static vector g_OPBinLimits(2, 0);
static vector g_WValues(1, 0);
static int_vector g_N_OP_Bin(1, 0);
static int_vector g_N_OP_BinTotal(1, 0);
static int g_WeightIterationCount = 0;
static bool g_WeightIterationFlag = true;

// Pointer to the iteration function
iteration_pointer iterate_weights;

// Pointer to the finish condition check
finish_condition_pointer finish_check;

////////////////////////////////////////////////////////////////////////////////
/// @brief Writes a variable to the file, given the format string.
///
/// @param output_file
/// @param fmt           format string corresponding to input_value
/// @param input_value   numerical value to write to output_file
////////////////////////////////////////////////////////////////////////////////
template <class K>
void to_file(std::ofstream &output_file, string fmt, K input_value) {
    char buffer[1024];
    sprintf(buffer, fmt.c_str(), input_value);
    if (hila::myrank() == 0)
        output_file << string(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Generates a time stamped and otherwise appropriate file name for the
///        saved weight function files.
///
/// @return      generated filename string
////////////////////////////////////////////////////////////////////////////////
string generate_outfile_name() {
    string filename = g_WParam.outfile_name_base + "_weight_function_";

    // A terrible mess to get the datetime format nice
    std::stringstream ss;
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    ss << std::put_time(&tm, "created_%Y.%m.%d_%H:%M:%S");
    string date = ss.str();

    filename = filename + date;
    return filename;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Parses the weight parameter file and fills the g_WParam struct.
/// @details Cannot fail. If parameter file is not of the correct format hila::finishrun() kills the
/// process.
///
/// @param parameter_file_name   parameter file name
////////////////////////////////////////////////////////////////////////////////
void read_weight_parameters(string parameter_file_name) {
    // Open the weight parameter file and list through the parameters.
    // See the parameter file for the roles of the parameters.
    hila::input par(parameter_file_name);

    // Generic control parameters
    string output_loc = par.get("output file location");
    string outfile_name_base = par.get("output file name base");

    string weight_loc = par.get("weight file location");
    string iter_method = par.get("iteration method");
    string hard_walls = par.get("hard walls");
    double max_OP = par.get("max OP");
    double min_OP = par.get("min OP");
    int bin_number = par.get("bin number");
    string iter_vis = par.get("iteration visuals");

    // Direct iteration parameters
    string finish_condition = par.get("finish condition");
    int DIM_sample_size = par.get("DIM sample size");
    int DIM_check_interval = par.get("DIM visit check interval");
    double add_initial = par.get("add initial");
    double add_minimum = par.get("add minimum");

    // Canonical iteration parameters
    int CIM_sample_size = par.get("CIM sample size");
    int initial_bin_hits = par.get("initial bin hits");
    int OC_max_iter = par.get("OC max iter");
    int OC_frequency = par.get("OC frequency");

    par.close();

    // clang-format off
    struct direct_iteration DIP
    {
        finish_condition,
        DIM_sample_size,
        DIM_check_interval,
        add_initial,
        add_minimum,
        add_initial
    };

    struct canonical_iteration CIP
    {
        CIM_sample_size,
        initial_bin_hits,
        OC_max_iter,
        OC_frequency
    };
    // clang-format on

    bool AR_ITER = false;

    bool visuals;
    if (iter_vis.compare("YES") == 0)
        visuals = true;
    else
        visuals = false;

    bool hwalls;
    if (hard_walls.compare("YES") == 0)
        hwalls = true;
    else
        hwalls = false;

    // clang-format off
    g_WParam =
    {
        weight_loc,
        outfile_name_base,
        iter_method,
        visuals,
        hwalls,
        max_OP,
        min_OP,
        bin_number,
        AR_ITER,
        DIP,
        CIP
    };
    // clang-format on
}

////////////////////////////////////////////////////////////////////////////////
/// @brief
/// Reads a precomputed weight function from file.
/// @details
/// The input file is to have three tab-separated columns:
/// 1. bin edges, of length n + 1
/// 2. bin centres, of length n
/// 3. weight values, of length n
///
/// Naturally, column 2 can be determined from the first one
/// but they are considered a separated input anyways.
///
/// The header can be whatever* and is always skipped. Regex finds the
/// data by finding first row with the substring "OP_value" and
/// assumes that the following contains the data as specified above.
///
/// *Avoid substring "OP_value"
///
/// Should be called only from hila::myrank()==0
///
/// @param W_function_filename
/// @return false if reading the file fails.
////////////////////////////////////////////////////////////////////////////////
bool read_weight_function(string W_function_filename) {
    printf("Loading the user supplied weight function `%s`\n", W_function_filename.c_str());

    int header_length = 1, data_length = -1;
    std::ifstream W_file;
    W_file.open(W_function_filename.c_str());

    if (!W_file.is_open()) {
        printf("ERROR: Could not open file `%s`: %s\n", W_function_filename.c_str(),
               std::strerror(errno));
        return false;
    }

    // Compute first header length by counting lines until finding
    // the column header through regex.
    string line;
    while (std::getline(W_file, line)) {
        if (std::regex_match(line, std::regex(".*OP_value.*")))
            data_length = 0;
        if (data_length < 0)
            header_length += 1;
        else
            data_length += 1;
    }
    if (!W_file.eof() && W_file.fail()) {
        printf("ERROR:(1) Failed to getline on file `%s`: %s\n", W_function_filename.c_str(),
               std::strerror(errno));
        return false;
    }

    // goto the beginning of the file
    W_file.clear(); // clear error bit set by above reaching eof
    W_file.seekg(std::ios_base::beg);

    // Skip the header and sscanf the values and add them into the vectors.
    int count = -header_length;
    data_length -= 1;
    printf("Weight function has header length of %d rows.\n", header_length);
    printf("Weight function has data length of %d rows.\n", data_length);
    printf("Reading the weight function into the program.\n");

    // Initialise the weight vectors to correct dimensions
    g_OPBinLimits = vector(data_length);
    g_OPValues = vector(data_length - 1);
    g_WValues = vector(data_length - 1);

    // Read in the values. Note that g_OPBinLimits has one more entry than
    // the weight vector.
    int line_no = 0;
    while (std::getline(W_file, line)) {
        line_no++;
        float bin_lim, weight, centre;
        // Read lower bin limits and corresponding weights.
        if (count >= 0 and count < data_length - 1) {
            if (3 != sscanf(line.c_str(), "%e\t%e\t%e", &bin_lim, &centre, &weight)) {
                printf("ERROR: Could not match `%%e\\t%%e\\t%%e` at %s:%d \n",
                       W_function_filename.c_str(), line_no);
                return false;
            }

            g_OPBinLimits[count] = bin_lim;
            g_OPValues[count] = centre;
            g_WValues[count] = weight;
        }
        // And the rightmost bin limit.
        if (count == data_length - 1) {
            if (1 != sscanf(line.c_str(), "%e", &bin_lim)) {
                printf("ERROR: Could not match `%%e` at %s:%d \n", W_function_filename.c_str(),
                       line_no);
                return false;
            }

            g_OPBinLimits[count] = bin_lim;
        }
        count += 1;
    }
    if (!W_file.eof() && W_file.fail()) {
        printf("ERROR:(2) Failed to getline on file %s:%d %s\n", W_function_filename.c_str(),
               line_no, std::strerror(errno));
        return false;
    }

    W_file.close(); // Lets not care if closing fails...

    // Fill out bin centres for the interpolator.
    for (int i = 0; i < data_length - 1; i++) {
        g_OPValues[i] = (g_OPBinLimits[i + 1] + g_OPBinLimits[i]) / 2.0;
    }
    printf("Succesfully loaded the user provided weight function.\n");
    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Reads the precomputed weight function from run_parameters struct
///        and saves it into a file.
/// @details The printing happens in an format identical to what is expected
///          by the funciton read_weight_function. See its documentation for
///          details.
///          TBA: Add string input that can contain user specified header data.
///
/// @param W_function_filename
/// @param g_WParam                    struct of weight iteration parameters
////////////////////////////////////////////////////////////////////////////////
void write_weight_function(string W_function_filename) {
    if (hila::myrank() == 0) {
        std::ofstream W_file;
        // string filename = generate_outfile_name(RP);
        W_file.open(W_function_filename.c_str());
        // write_weight_file_header(W_file, FP, RP, g_WParam);

        to_file(W_file, "%s", "OP_bin_limit\tOP_value\tWeight\n");
        int i;
        for (i = 0; i < g_OPValues.size(); ++i) {
            to_file(W_file, "%e\t", g_OPBinLimits[i]);
            to_file(W_file, "%e\t", g_OPValues[i]);
            to_file(W_file, "%e\n", g_WValues[i]);
        }
        // Remember to write the last bin upper limit
        to_file(W_file, "%e\n", g_OPBinLimits[i]);
        W_file.close();
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Returns a weight associated to the used order parameter.
/// @details The function uses supplied pairs of points to linearly interpolate
///          the function on the interval. This interpolant provides the
///          requested weights to be used as the multicanonical weight.
///
/// @param  OP   value of the order parameter
/// @return The value of the weight.
////////////////////////////////////////////////////////////////////////////////
double weight_function(double OP) {
    double val;
    // If out of range, constant extrapolation or for hard walls, num inf.
    if (OP <= g_OPValues.front()) {
        if (g_WParam.hard_walls)
            val = std::numeric_limits<double>::infinity();
        else
            val = g_WValues.front();
    } else if (OP >= g_OPValues.back()) {
        if (g_WParam.hard_walls)
            val = std::numeric_limits<double>::infinity();
        else
            val = g_WValues.back();
    }
    // Otherwise find interval, calculate slope, base index, etc.
    // Basic linear interpolation to obtain the weight value.
    else {
        auto it = std::lower_bound(g_OPValues.begin(), g_OPValues.end(), OP);
        int i = std::distance(g_OPValues.begin(), it) - 1;
        double y_value = g_WValues[i + 1] - g_WValues[i];
        double x_value = g_OPValues[i + 1] - g_OPValues[i];
        double slope = y_value / x_value;

        double xbase = g_OPValues[i];
        double ybase = g_WValues[i];

        double xdiff = OP - xbase;
        val = ybase + xdiff * slope;
    }

    return val;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief process 0 interface to "weight function" for the user accessing
///        the weights.
////////////////////////////////////////////////////////////////////////////////
double weight(double OP) {
    double val;
    if (hila::myrank() == 0) {
        val = weight_function(OP);
    }
    hila::broadcast(val);
    return val;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Sets the static g_WeightIterationFlag to given boolean.
///
/// @param YN   boolean indicating whether the iteration is to continue
////////////////////////////////////////////////////////////////////////////////
void set_weight_iter_flag(bool YN) {
    if (hila::myrank() == 0)
        g_WeightIterationFlag = YN;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Returns the value of the static g_WeightIterationFlag to user
///
/// @return State of g_WeighITerationFlag
////////////////////////////////////////////////////////////////////////////////
bool check_weight_iter_flag() {
    bool flag;
    if (hila::myrank() == 0)
        flag = g_WeightIterationFlag;
    hila::broadcast(flag);
    return flag;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Accepts/rejects a multicanonical update.
/// @details Using the values of the old and new order parameters the muca
///          update is accepted with the logarithmic probability
///          log(P) = - (W(OP_new) - W(OP_old))
///
/// @param  OP_old     current order parameter
/// @param  OP_new     order parameter of proposed configuration
/// @return Boolean indicating whether the update was accepted (true) or
///         rejected (false).
////////////////////////////////////////////////////////////////////////////////
bool accept_reject(const double OP_old, const double OP_new) {
    bool update;
    bool AR_iterate;
    // Only compute on node 0, broadcast to others
    if (hila::myrank() == 0) {
        double W_new = weight_function(OP_new);
        double W_old = weight_function(OP_old);

        // get log(exp(-delta(W))) = -delta(W)
        // (just like -delta(S) in Metropolis-Hastings)
        double log_P = -(W_new - W_old);

        // Get a random uniform from [0,1) and return a boolean indicating
        // whether the update is to be accepted.
        double rval = hila::random();
        if (::log(rval) < log_P) {
            update = true;
        } else {
            update = false;
        }

        // Get value from process 0
        AR_iterate = g_WParam.AR_iteration;
    }

    // Broadcast the update status to other processes along with the
    // weight iteration parameter
    hila::broadcast(update);
    hila::broadcast(AR_iterate);

    // Check if iteration is enabled
    if (AR_iterate) {
        bool continue_iter = iterate_weights(update ? OP_new : OP_old);
        set_weight_iter_flag(continue_iter);
    }

    return update;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Finds the index of the correc order parameter bin.
/// @details Using the bin limit vector the correct order parameter bin is
///          determined through a simple standard library search of
///          g_OPBinLimits. When the value of the given order parameter is
///          out of range, -1 is returned.
///
/// @param OP   value of the order parameter
/// @return integer index for the vector g_N_OP_Bin
////////////////////////////////////////////////////////////////////////////////
static int find_OP_bin_index(double OP) {
    // Return -1 when not in the interval
    if (OP <= g_OPBinLimits.front()) {
        return -1;
    } else if (OP >= g_OPBinLimits.back()) {
        return -1;
    }
    // Find index of minimum edge value such that edge < OP:
    auto it = std::lower_bound(g_OPBinLimits.begin(), g_OPBinLimits.end(), OP);
    int lower_limit = std::distance(g_OPBinLimits.begin(), it) - 1;
    return lower_limit;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Same as find_OP_bin_index, except uses the index to simply modify the
///        bin hit vector g_N_OP_Bin. Does not modify when outside of the range.
/// @details
///
/// @param OP   value of the order parameter
////////////////////////////////////////////////////////////////////////////////
static void bin_OP_value(double OP) {
    // Don't bin visits outside of the binned areas
    if (OP <= g_OPBinLimits.front()) {
        return;
    } else if (OP >= g_OPBinLimits.back()) {
        return;
    }
    // Find index of minimum edge value such that edge < OP:
    auto it = std::lower_bound(g_OPBinLimits.begin(), g_OPBinLimits.end(), OP);
    int lower_limit = std::distance(g_OPBinLimits.begin(), it) - 1;
    g_N_OP_Bin[lower_limit] += 1;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Checks if all the bins have been visited by.
/// @details Simply checks whether all bins have a nonzero number of entries
///
/// @param  visit   integer vector with values 1 corresponding to visits
/// @return a boolean indicating the statement
////////////////////////////////////////////////////////////////////////////////
static bool all_visited(int_vector &n) {
    int len = n.size();
    for (int i = 0; i < len; ++i) {
        if (n[i] == 0)
            return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Checks if the first and last bin have been visited
/// @details Simply checks whether all bins have a nonzero number of entries.
///
/// @param  visit   integer vector with values 1 corresponding to visits
/// @return a boolean indicating the statement
////////////////////////////////////////////////////////////////////////////////
static bool first_last_visited(int_vector &n) {
    int len = n.size();
    if ((n[0] == 0) or (n[len - 1] == 0))
        return false;
    else
        return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Sets a user provided function to the check in the "direct iteration"
/// method.
/// @details The for a given magnitude of update the "direct iteration" method
/// periodically checks whether the MCMC chain has covered enough of the
/// desired order parameter range, before reducing the update magnitude. Some
/// preset methods exist (and should suffice) but when needed, a new condition
/// can be set through this function. The input is a vector of integers
/// indicating the number of visits to each bin, and the output is a boolean
/// telling whether the desired condition has been achieved.
///
/// @param fc_pointer   A function pointer to a suitable condition function
////////////////////////////////////////////////////////////////////////////////
void set_direct_iteration_FC(bool (*fc_pointer)(int_vector &n)) {
    finish_check = fc_pointer;
}

////////////////////////////////////////////////////////////////////////////////
// Following three functions related to "canonical" iteration are currently
// not implemented properly.
////////////////////////////////////////////////////////////////////////////////
/// @brief Computes an overcorrection of W given the visits ns (NOT FUNCTIONAL).
/// @details Overcorrection checks which bins have been visited the most and
///          modifies the weights in such a manner that the frequently visited
///          bins are weighted (disproportionally) heavily. This should
///          encourage the process to visit other bins during the next run of
///          iterations.
///
/// @param Weight   vector containing the previously used weights
/// @param n_sum    vector containing the total number of visits
////////////////////////////////////////////////////////////////////////////////
static void overcorrect(vector &Weight, int_vector &n_sum) {
    int N = n_sum.size();
    vector W(N, 0);
    constexpr float C = 1;

    for (int m = 0; m < N; ++m) {
        if (n_sum[m] > 0)
            W[m] = C * ::log(n_sum[m]);
    }

    // Redefine W[0] back to zero, set new weights
    double base = W[0];
    for (int i = 0; i < N; ++i) {
        Weight[i] += W[i] - base;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Computes the next weight vector through iteration based on previous
///        weights (NOT FUNCTIONAL).
/// @details
///
/// @param Weight        vector containing the previously used weights
/// @param n             vector containing the current number of visits
/// @param g_sum         vector containing the previous local sums of n
/// @param g_log_h_sum   vector containing the sum of previous logs of the
///                      canonical weights
////////////////////////////////////////////////////////////////////////////////
static void recursive_weight_iteration(vector &Weight, int_vector &n, int_vector &g_sum,
                                       vector &g_log_h_sum) {
    const int nmin = 10;
    int N = n.size();

    // Fill out log(h)
    vector log_h(N);
    for (int m = 0; m < N; ++m) {
        if (n[m] > 0) {
            log_h[m] = ::log(n[m]) + Weight[m];
        }
        // This is only going to get multiplied by zero,
        // so the value doesn't matter.
        else
            log_h[m] = 0;
    }

    // Compute the iteration. Note that W[0] is permanently kept to zero
    vector W(N);
    W[0] = 0;
    Weight[0] = W[0];
    for (int m = 1; m < N; ++m) {
        int gm;
        // Check that both bins contain nmin hits before taking into account
        if (n[m] > nmin and n[m - 1] > nmin) {
            gm = n[m] + n[m - 1];
        } else
            gm = 0;

        g_sum[m] += gm;
        g_log_h_sum[m] += gm * (log_h[m] - log_h[m - 1]);

        W[m] = W[m - 1] + g_log_h_sum[m] / g_sum[m];

        // Modify the input weights to their new values
        Weight[m] = W[m];
    }
}

//////////////////////////////////////////////////////////////////////////////////
///// @brief Given an order parameter, iterates the weight function until the
/////        sampling becomes acceptably efficient. Save the weight function into
/////        a file for later use (NOT FUNCTIONAL).
///// @details
/////
///// @param F     struct of fields
///// @param FP    struct of field parameters
///// @param RP    struct of run parameters
///// @param g_WParam   struct of weight iteration parameters
//////////////////////////////////////////////////////////////////////////////////
// void iterate_weight_function_canonical(fields &F,
//                                        field_parameters &FP,
//                                        run_parameters &RP,
//                                        weight_iteration_parameters &g_WParam)
//{
//     int samples  = g_WParam.sample_size;
//     int n_sweeps = g_WParam.sample_steps;
//
//     double max = g_WParam.max_OP;
//     double min = g_WParam.min_OP;
//
//     int N = g_WParam.bin_number;
//     // Initialise the weight vector with the bin centres
//     // and get corresponding bin limits.
//     {
//         double dx = (max - min) / (N - 1);
//         for (int i = 0; i < N; ++i)
//         {
//             RP.OP_values[i] = min + i * dx;
//         }
//     }
//     vector limits = get_bin_limits(min, max, N);
//
//     // Initialise the storage vectors to zero:
//     int_vector n(N, 0), g_sum(N, 0), n_sum(N, 0);
//     vector g_log_h_sum(N, 0), log_h(N, 0), W_prev(N, 0);
//
//     // Get initial guesses
//     for (int i = 0; i < N; ++i)
//     {
//         n[i]     = 0;
//         g_sum[i] = g_WParam.initial_bin_hits;
//         log_h[i] = RP.W_values[i];
//     }
//
//     static int count = 1;
//     while (true)
//     {
//         float accept = 0;
//         float OP_sum = 0;
//         for (int i = 0; i < samples; i++)
//         {
//             accept += mc_update_sweeps(F, FP, RP, n_sweeps);
//             OP_sum += FP.OP_value;
//             bin_OP_values(n, limits, FP.OP_value);
//         }
//
//         for (int m = 0; m < N; m++)
//         {
//             n_sum[m] += n[m];
//         }
//
//         if (count % g_WParam.OC_frequency == 0 and count < g_WParam.OC_max_iter)
//         {
//             overcorrect(RP.W_values, n_sum);
//         }
//         else
//         {
//             recursive_weight_iteration(RP.W_values, n, g_sum, g_log_h_sum);
//         }
//
//
//         // Some mildly useful print out
//         int nmax = *std::max_element(n_sum.begin(), n_sum.end());
//         for (int m = 0; m < N; ++m)
//         {
//             std::string n_sum_log = "";
//             for (int i = 0; i < int(n_sum[m] * 50.0 / nmax); i++)
//             {
//                 n_sum_log += "|";
//             }
//             if (hila::myrank() == 0)
//             {
//                 printf("%5.3f\t%10.3f\t\t%d\t%s\n", limits[m],
//                         RP.W_values[m],
//                         n_sum[m], n_sum_log.c_str());
//             }
//             n[m] = 0;
//         }
//
//         count += 1;
//         if (count > g_WParam.max_iter) break;
//     }
// }

////////////////////////////////////////////////////////////////////////////////
/// @brief Procures a vector containing equidistant bin edges.
/// @details
///
/// @return vector containing the bin edges
////////////////////////////////////////////////////////////////////////////////
static vector get_equidistant_bin_limits() {
    double min = g_WParam.min_OP;
    double max = g_WParam.max_OP;
    int N = g_WParam.bin_number;
    vector bin_edges(N + 1);
    double diff = (max - min) / (N - 1);
    for (int i = 0; i < N + 1; ++i) {
        bin_edges[i] = min - diff / 2.0 + diff * i;
    }
    return bin_edges;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Sets up the global vectors for bin limits and centres using
///        get_equidistant_bin_limits.
/// @details
////////////////////////////////////////////////////////////////////////////////
static void setup_equidistant_bins() {
    // Get bin limits so that centre of first bin is at min_OP and
    // the last bin centre is at max_OP.
    g_OPBinLimits = get_equidistant_bin_limits();
    for (int i = 0; i < g_OPValues.size(); i++) {
        double centre = (g_OPBinLimits[i + 1] + g_OPBinLimits[i]) / 2.0;
        g_OPValues[i] = centre;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Initialises the global vectors appropriately, setting up a binning
///        if not provided by the user.
/// @details The global vectors are initialised to correct dimensions as to
///          prevent indexing errors in the iteration methods.
////////////////////////////////////////////////////////////////////////////////
static void initialise_weight_vectors() {
    // If no input weight, set up equidistant bins
    if (g_WParam.weight_loc.compare("NONE") == 0) {
        int N = g_WParam.bin_number;
        g_WValues = vector(N, 0.0);
        g_OPValues = vector(N, 0.0);
        g_OPBinLimits = vector(N + 1, 0.0);

        setup_equidistant_bins();

        g_N_OP_Bin = int_vector(N, 0);
        g_N_OP_BinTotal = int_vector(N, 0);
    }
    // Same for predetermined bins. g_OPValues, g_WValues
    // and g_OP_BinLimits have been read from the input file.
    // To prevent any accidents with the iterators, these bin vectors
    // are initialised in all cases. Should really not affect performance.
    else {
        int N = g_WValues.size();
        g_N_OP_Bin = int_vector(N, 0);
        g_N_OP_BinTotal = int_vector(N, 0);
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Given an order parameter, bins it to correct weight interval, and
///        periodically updates the weights accordingly.
/// @details This extremely simple update method
///
/// @param  OP   order parameter of the current configuration (user supplied)
/// @return boolean indicating whether the iteration is considered complete
////////////////////////////////////////////////////////////////////////////////
static bool iterate_weight_function_direct(double OP) {
    bool continue_iteration;
    if (hila::myrank() == 0) {
        int samples = g_WParam.DIP.sample_size;
        int N = g_WValues.size();

        bin_OP_value(OP);
        g_WeightIterationCount += 1;

        if (g_WeightIterationCount >= samples) {
            for (int m = 0; m < N; m++) {
                g_WValues[m] += g_WParam.DIP.C * g_N_OP_Bin[m] * N / samples;
                g_N_OP_BinTotal[m] += g_N_OP_Bin[m];
            }

            if (g_WParam.visuals)
                print_iteration_histogram();

            double base = *std::min_element(g_WValues.begin(), g_WValues.end());
            for (int m = 0; m < N; ++m) {
                // Always set minimum weight to zero. This is inconsequential
                // as only the differences matter.
                g_WValues[m] -= base;
                g_N_OP_Bin[m] = 0;
            }
            g_WeightIterationCount = 0;

            if (finish_check(g_N_OP_BinTotal)) {
                for (int m = 0; m < N; m++) {
                    g_N_OP_BinTotal[m] = 0;
                }

                g_WParam.DIP.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = " << g_WParam.DIP.C
                           << "\n";
            }
            write_weight_function("intermediate_weight.dat");
        }

        continue_iteration = true;
        if (g_WParam.DIP.C < g_WParam.DIP.C_min) {
            hila::out0 << "Muca: Reached minimum update size C = " << g_WParam.DIP.C
                       << " Weight iteration complete.\n";
            continue_iteration = false;
        }
    }
    hila::broadcast(continue_iteration);
    return continue_iteration;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Same as iterate_weight_function_direct for sample size 1.
/// @details In this special case the weights are modified after each new
///          value of the order parameter. Reduced internal complexity due to
///          this simplification.
///          Whether all bins have been visited is now checked only every
///          g_WParam.DIP.single_check_interval to prevent excessive checking
///          for the visits.
///
/// @param  OP   order parameter of the current configuration (user supplied)
/// @return boolean indicating whether the iteration is considered complete
////////////////////////////////////////////////////////////////////////////////
static bool iterate_weight_function_direct_single(double OP) {
    int continue_iteration;
    if (hila::myrank() == 0) {
        int samples = g_WParam.DIP.sample_size;
        int N = g_WValues.size();

        int bin_index = find_OP_bin_index(OP);
        // Only increment if on the min-max interval
        if (bin_index != -1)
            g_N_OP_BinTotal[bin_index] += 1;

        g_WeightIterationCount += 1;

        g_WValues[bin_index] += g_WParam.DIP.C;

        if (g_WeightIterationCount % g_WParam.DIP.single_check_interval == 0) {

            double base = *std::min_element(g_WValues.begin(), g_WValues.end());
            for (int m = 0; m < N; ++m) {
                // Always set minimum weight to zero. This is inconsequential
                // as only the differences matter.
                g_WValues[m] -= base;
                g_N_OP_Bin[m] = 0;
            }

            // Visuals
            if (g_WParam.visuals)
                print_iteration_histogram();

            // If condition satisfied, zero the totals and decrease C
            if (finish_check(g_N_OP_BinTotal)) {
                for (int m = 0; m < N; m++) {
                    g_N_OP_BinTotal[m] = 0;
                }

                g_WParam.DIP.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = " << g_WParam.DIP.C
                           << "\n";
            }

            continue_iteration = true;
            if (g_WParam.DIP.C < g_WParam.DIP.C_min) {
                hila::out0 << "Muca: Reached minimum update size C = " << g_WParam.DIP.C
                           << " Weight iteration complete.\n";
                continue_iteration = false;
            }

            write_weight_function("intermediate_weight.dat");
        }
    }
    hila::broadcast(continue_iteration);
    return continue_iteration;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Prints out a crude horisontal histogram.
/// @details Procures a crude horisontal ASCII histogram based on the g_N_OP_Bin
///          vector. The histogram bin heights are proportional to g_N_OP_Bin
///          values. This is not very expressive for large N, as it won't fit
///          the height of the screen.
////////////////////////////////////////////////////////////////////////////////
static void print_iteration_histogram() {
    int samples = g_WParam.DIP.sample_size;
    int N = g_WValues.size();
    // Find maximum bin content for normalisation
    int nmax = *std::max_element(g_N_OP_Bin.begin(), g_N_OP_Bin.end());
    // Write a column header
    printf("Order Parameter     Weight 	         Number of hits\n");

    for (int m = 0; m < N; ++m) {
        // For each bin get a number of "|":s proportional to the number of
        // hits to each bin and print it out along with relevant numerical
        // values
        std::string n_sum_hist = "";
        if (g_N_OP_BinTotal[m] > 0)
            n_sum_hist += "O";
        for (int i = 0; i < int(g_N_OP_Bin[m] * 200.0 / samples); i++) {
            n_sum_hist += "|";
        }
        printf("%-20.3f%-20.3f%d\t\t\t%s\n", g_OPValues[m], g_WValues[m], g_N_OP_Bin[m],
               n_sum_hist.c_str());
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Initialises all the variables needed for the weight iteration.
/// @details For the iteration the function pointer to iterate_weights must be
///          initialised. The condition for proceeding to the next step of the
///          iteration is also set through a function pointer.
///
///          The correct methods are chosen according to the
///          parameter file. Further, method specific variables that are
///          modified during the run are also set to the necessary values.
////////////////////////////////////////////////////////////////////////////////
static void setup_iteration() {
    // Initialise iterate_weights by pointing it at the
    // correct method
    if (g_WParam.method.compare("direct") == 0) {
        if (g_WParam.DIP.sample_size > 1) {
            iterate_weights = &iterate_weight_function_direct;
        } else {
            iterate_weights = &iterate_weight_function_direct_single;
        }
        g_WParam.DIP.C = g_WParam.DIP.C_init;
    } else {
        iterate_weights = &iterate_weight_function_direct;
        g_WParam.DIP.C = g_WParam.DIP.C_init;
    }

    // Zero the iteration counter
    g_WeightIterationCount = 0;

    // Set up the finish condition pointer for the direct method.
    if (g_WParam.DIP.finish_condition.compare("all_visited") == 0) {
        finish_check = &all_visited;
    } else if (g_WParam.DIP.finish_condition.compare("ends_visited") == 0) {
        finish_check = &first_last_visited;
    } else {
        finish_check = &all_visited;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Enable/disable continuous weight iteration
/// @details Premits the user to enable/disable continuous weight iteration at
///          each call to accept_reject. Simply modifies a flag parameter
///          that is checked in accept_reject.
///
/// @param YN   enable (true) or disable (false) the iteration
////////////////////////////////////////////////////////////////////////////////
void set_continuous_iteration(bool YN) {
    if (hila::myrank() == 0)
        g_WParam.AR_iteration = YN;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Loads parameters and weights for the multicanonical computation.
/// @details Sets up iteration variables. Can be called multiple times and must
///          be called at least once before attempting to use any of the muca
///          methods.
///
/// @param wfile_name   path to the weight parameter file
/// @return false if something fails during init.
////////////////////////////////////////////////////////////////////////////////
bool initialise(const string wfile_name) {
    // Read parameters into g_WParam struct
    read_weight_parameters(wfile_name);
    bool ret = true;
    // This is fine to do just for process zer0
    if (hila::myrank() == 0) {
        // Read pre-existing weight if given
        if (g_WParam.weight_loc.compare("NONE") != 0) {
            if (!read_weight_function(g_WParam.weight_loc)) {
                ret = false;
                goto early_return;
            }
        }
        // Initialise rest of the uninitialised vectors
        initialise_weight_vectors();
    }

    // Choose an iteration method (or the default)
    setup_iteration();

early_return:
    hila::broadcast(ret);
    return ret;
}

} // namespace muca
} // namespace hila
