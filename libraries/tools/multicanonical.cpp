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

////////////////////////////////////////////////////////////////////////////////
/// @brief Writes variables to the file, given the format std::string.
///        Quick and dirty. No error handling.
/// TODO: swith to std::format when c++20
///
/// @param output_file
/// @param fmt           format std::string corresponding to input_value
/// @param input_values   numerical values to write to output_file
////////////////////////////////////////////////////////////////////////////////
template <class... Ks>
static inline void to_file(std::ofstream &output_file, std::string fmt, Ks... input_values) {

    constexpr int buf_len = 1024;
    char buffer[buf_len];
    int res = snprintf(buffer, buf_len, fmt.c_str(), input_values...);
    if (res < 0) {
        printf("WARNING: to_file snprintf error %d\n", res);
    } else if (res >= buf_len) {
        printf("WARNING: to_file snprintf truncated %d\n", res);
    }

    output_file << std::string(buffer);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Generates a time stamped and otherwise appropriate file name for the
///        saved weight function files.
///
/// @return      generated filename std::string
////////////////////////////////////////////////////////////////////////////////
std::string generate_outfile_name(const std::string &outfile_name_base) {
    std::string filename = outfile_name_base + "_weight_function_";

    // A terrible mess to get the datetime format nice
    std::stringstream ss;
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    ss << std::put_time(&tm, "created_%Y.%m.%d_%H:%M:%S");
    std::string date = ss.str();

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
void Muca::read_weight_parameters(std::string parameter_file_name) {
    // Open the weight parameter file and list through the parameters.
    // See the parameter file for the roles of the parameters.
    hila::input par(parameter_file_name);

    // clang-format off
    // Generic control parameters
    std::string output_loc        = par.get("output file location"); // Unused

    std::string outfile_name_base = par.get("output file name base");
    std::string weight_loc        = par.get("weight file location");
    std::string iter_method       = par.get("iteration method");
    std::string hard_walls        = par.get("hard walls");
    double max_OP                 = par.get("max OP");
    double min_OP                 = par.get("min OP");
    int bin_number                = par.get("bin number");
    std::string iter_vis          = par.get("iteration visuals");

    // Direct iteration parameters
    std::string finish_condition = par.get("finish condition");
    int DIM_sample_size          = par.get("DIM sample size");
    int DIM_check_interval       = par.get("DIM visit check interval");
    double add_initial           = par.get("add initial");
    double add_minimum           = par.get("add minimum");
    struct direct_iteration DIP
    {
        finish_condition,
        DIM_sample_size,
        DIM_check_interval,
        add_initial,
        add_minimum,
        add_initial
    };

    // Canonical iteration parameters
    int CIM_sample_size  = par.get("CIM sample size");
    int initial_bin_hits = par.get("initial bin hits");
    int OC_max_iter      = par.get("OC max iter");
    int OC_frequency     = par.get("OC frequency");
    struct canonical_iteration CIP
    {
        CIM_sample_size,
        initial_bin_hits,
        OC_max_iter,
        OC_frequency
    };

    par.close();

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

    WParam =
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
bool Muca::read_weight_function(const std::string &W_function_filename) {
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
    std::string line;
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
    OPBinLimits = std::vector<double>(data_length);
    OPValues = std::vector<double>(data_length - 1);
    WValues = std::vector<double>(data_length - 1);

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

            OPBinLimits[count] = bin_lim;
            OPValues[count] = centre;
            WValues[count] = weight;
        }
        // And the rightmost bin limit.
        if (count == data_length - 1) {
            if (1 != sscanf(line.c_str(), "%e", &bin_lim)) {
                printf("ERROR: Could not match `%%e` at %s:%d \n", W_function_filename.c_str(),
                       line_no);
                return false;
            }

            OPBinLimits[count] = bin_lim;
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
        OPValues[i] = (OPBinLimits[i + 1] + OPBinLimits[i]) / 2.0;
    }
    printf("Succesfully loaded the user provided weight function.\n");
    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Reads the precomputed weight function from run_parameters struct
///        and saves it into a file.
/// @details The printing happens in an format identical to what is expected
///          by the function read_weight_function. See its documentation for
///          details.
///          Should only be called from hila::myrank()==0
///
///          TBA: Add std::string input that can contain user specified header data.
///
/// @param W_function_filename
/// @param g_WParam                    struct of weight iteration parameters
/// @return false if writing fails
////////////////////////////////////////////////////////////////////////////////
bool Muca::write_weight_function(const std::string &W_function_filename) {

    std::ofstream W_file;
    // std::string filename = generate_outfile_name(RP);
    W_file.open(W_function_filename.c_str());
    if (!W_file.is_open()) {
        printf("WARNING: Could not open file `%s` for `write_weight_function()`: %s\n",
               W_function_filename.c_str(), std::strerror(errno));
        return false;
    }
    // write_weight_file_header(W_file, FP, RP, g_WParam);

    to_file(W_file, "%s", "OP_bin_limit\tOP_value\tWeight\n");
    int i;
    for (i = 0; i < OPValues.size(); ++i) {
        to_file(W_file, "%e\t%e\t%e\n", OPBinLimits[i], OPValues[i], WValues[i]);
    }
    // Remember to write the last bin upper limit
    to_file(W_file, "%e\n", OPBinLimits[i]);
    W_file.close();
    return true;
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
double Muca::weight_function(double OP) const {
    // the Muca struct is to be initialised only on node 0
    if (hila::myrank() != 0)
        return 0;

    double val;
    // If out of range, constant extrapolation or for hard walls, num inf.
    if (OP <= OPValues.front()) {
        if (WParam.hard_walls)
            val = std::numeric_limits<double>::infinity();
        else
            val = WValues.front();
    } else if (OP >= OPValues.back()) {
        if (WParam.hard_walls)
            val = std::numeric_limits<double>::infinity();
        else
            val = WValues.back();
    }
    // Otherwise find interval, calculate slope, base index, etc.
    // Basic linear interpolation to obtain the weight value.
    else {
        auto it = std::lower_bound(OPValues.begin(), OPValues.end(), OP);
        int i = std::distance(OPValues.begin(), it) - 1;
        double y_value = WValues[i + 1] - WValues[i];
        double x_value = OPValues[i + 1] - OPValues[i];
        double slope = y_value / x_value;

        double xbase = OPValues[i];
        double ybase = WValues[i];

        double xdiff = OP - xbase;
        val = ybase + xdiff * slope;
    }

    return val;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief process 0 interface to "weight function" for the user accessing
///        the weights.
////////////////////////////////////////////////////////////////////////////////
double Muca::weight(double OP) const {
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
void Muca::set_weight_iter_flag(bool YN) {
    weightIterationFlag = YN;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Returns the value of the static g_WeightIterationFlag to user
///
/// @return State of g_WeighITerationFlag
////////////////////////////////////////////////////////////////////////////////
bool Muca::check_weight_iter_flag() {
    bool flag;
    if (hila::myrank() == 0)
        flag = weightIterationFlag;
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
bool Muca::accept_reject(const double OP_old, const double OP_new) {
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
        AR_iterate = WParam.AR_iteration;
    }

    // Broadcast the update status to other processes along with the
    // weight iteration parameter
    hila::broadcast(update);
    hila::broadcast(AR_iterate);

    // Check if iteration is enabled
    if (AR_iterate) {
        bool continue_iter = iterate_weights(*this, update ? OP_new : OP_old);
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
/// @param OPBinLimits   bin limits of the weight function hist
/// @return integer index for the vector g_N_OP_Bin
////////////////////////////////////////////////////////////////////////////////
static inline int find_OP_bin_index(double OP, const std::vector<double> &OPBinLimits) {
    // Return -1 when not in the interval
    if (OP <= OPBinLimits.front()) {
        return -1;
    } else if (OP >= OPBinLimits.back()) {
        return -1;
    }
    // Find index of minimum edge value such that edge < OP:
    auto it = std::lower_bound(OPBinLimits.begin(), OPBinLimits.end(), OP);
    int lower_limit = std::distance(OPBinLimits.begin(), it) - 1;
    return lower_limit;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Same as find_OP_bin_index, except uses the index to simply modify the
///        bin hit vector g_N_OP_Bin. Does not modify when outside of the range.
/// @details
///
/// @param OP   value of the order parameter
////////////////////////////////////////////////////////////////////////////////
inline void Muca::bin_OP_value(double OP) {
    // Don't bin visits outside of the binned areas
    if (OP <= OPBinLimits.front()) {
        return;
    } else if (OP >= OPBinLimits.back()) {
        return;
    }
    // Find index of minimum edge value such that edge < OP:
    auto it = std::lower_bound(OPBinLimits.begin(), OPBinLimits.end(), OP);
    int lower_limit = std::distance(OPBinLimits.begin(), it) - 1;
    N_OP_Bin[lower_limit] += 1;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Checks if all the bins have been visited by.
/// @details Simply checks whether all bins have a nonzero number of entries
///
/// @param  visit   integer vector with values 1 corresponding to visits
/// @return a boolean indicating the statement
////////////////////////////////////////////////////////////////////////////////
bool all_visited(const Muca &muca) {
    const auto &n = muca.N_OP_BinTotal;
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
bool first_last_visited(const Muca &muca) {
    const auto &n = muca.N_OP_BinTotal;
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
inline void Muca::set_direct_iteration_FC(finish_condition_fn fc) {
    finish_check = fc;
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
static void overcorrect(std::vector<double> &Weight, std::vector<int> &n_sum) {
    int N = n_sum.size();
    std::vector<double> W(N, 0);
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
static void recursive_weight_iteration(std::vector<double> &Weight, std::vector<int> &n,
                                       std::vector<int> &g_sum, std::vector<double> &g_log_h_sum) {
    const int nmin = 10;
    int N = n.size();

    // Fill out log(h)
    std::vector<double> log_h(N);
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
    std::vector<double> W(N);
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
//      std::vector<int> n(N, 0), g_sum(N, 0), n_sum(N, 0);
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
/// @param min weight function range minimum
/// @param max weight function range maximum
/// @param N_bins number of bins for the weight function
/// @return vector containing the bin edges
////////////////////////////////////////////////////////////////////////////////
static inline std::vector<double> get_equidistant_bin_limits(double min, double max, int N_bins) {
    std::vector<double> bin_edges(N_bins + 1);
    double diff = (max - min) / (N_bins - 1);
    for (int i = 0; i < N_bins + 1; ++i) {
        bin_edges[i] = min - diff / 2.0 + diff * i;
    }
    return bin_edges;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Sets up the global vectors for bin limits and centres using
///        get_equidistant_bin_limits.
/// @details
////////////////////////////////////////////////////////////////////////////////
inline void Muca::setup_equidistant_bins() {
    // Get bin limits so that centre of first bin is at min_OP and
    // the last bin centre is at max_OP.
    OPBinLimits = get_equidistant_bin_limits(WParam.min_OP, WParam.max_OP, WParam.bin_number);
    for (int i = 0; i < OPValues.size(); i++) {
        double centre = (OPBinLimits[i + 1] + OPBinLimits[i]) / 2.0;
        OPValues[i] = centre;
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Initialises the hist vectors appropriately, setting up a binning
///        if not provided by the user.
/// @details The vectors are initialised to correct dimensions as to
///          prevent indexing errors in the iteration methods.
////////////////////////////////////////////////////////////////////////////////
inline void Muca::initialise_weight_vectors() {
    // If no input weight, set up equidistant bins
    if (WParam.weight_loc.compare("NONE") == 0) {
        int N = WParam.bin_number;
        WValues = std::vector<double>(N, 0.0);
        OPValues = std::vector<double>(N, 0.0);
        OPBinLimits = std::vector<double>(N + 1, 0.0);

        setup_equidistant_bins();

        N_OP_Bin = std::vector<int>(N, 0);
        N_OP_BinTotal = std::vector<int>(N, 0);
    }
    // Same for predetermined bins. OPValues, WValues
    // and OP_BinLimits have been read from the input file.
    // To prevent any accidents with the iterators, these bin vectors
    // are initialised in all cases. Should really not affect performance.
    else {
        int N = WValues.size();
        N_OP_Bin = std::vector<int>(N, 0);
        N_OP_BinTotal = std::vector<int>(N, 0);
    }
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Given an order parameter, bins it to correct weight interval, and
///        periodically updates the weights accordingly.
/// @details This extremely simple update method
///
/// @param muca the multicanonical 'context'
/// @param  OP   order parameter of the current configuration (user supplied)
/// @return boolean indicating whether the iteration is considered complete
////////////////////////////////////////////////////////////////////////////////
bool iterate_weight_function_direct(Muca &muca, double OP) {
    bool continue_iteration;
    if (hila::myrank() == 0) {
        int samples = muca.WParam.DIP.sample_size;
        int N = muca.WValues.size();

        muca.bin_OP_value(OP);
        muca.weightIterationCount += 1;

        if (muca.weightIterationCount >= samples) {
            for (int m = 0; m < N; m++) {
                muca.WValues[m] += muca.WParam.DIP.C * muca.N_OP_Bin[m] * N / samples;
                muca.N_OP_BinTotal[m] += muca.N_OP_Bin[m];
            }

            if (muca.WParam.visuals)
                muca.print_iteration_histogram();

            double base = *std::min_element(muca.WValues.begin(), muca.WValues.end());
            for (int m = 0; m < N; ++m) {
                // Always set minimum weight to zero. This is inconsequential
                // as only the differences matter.
                muca.WValues[m] -= base;
                muca.N_OP_Bin[m] = 0;
            }
            muca.weightIterationCount = 0;

            if (muca.finish_check(muca)) {
                for (int m = 0; m < N; m++) {
                    muca.N_OP_BinTotal[m] = 0;
                }

                muca.WParam.DIP.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = "
                           << muca.WParam.DIP.C << "\n";
            }
            muca.write_weight_function("intermediate_weight.dat");
        }

        continue_iteration = true;
        if (muca.WParam.DIP.C < muca.WParam.DIP.C_min) {
            hila::out0 << "Muca: Reached minimum update size C = " << muca.WParam.DIP.C
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
///          WParam.DIP.single_check_interval to prevent excessive checking
///          for the visits.
///
/// @param  OP   order parameter of the current configuration (user supplied)
/// @return boolean indicating whether the iteration is considered complete
////////////////////////////////////////////////////////////////////////////////
bool iterate_weight_function_direct_single(Muca &muca, double OP) {
    int continue_iteration;
    if (hila::myrank() == 0) {
        int samples = muca.WParam.DIP.sample_size;
        int N = muca.WValues.size();

        int bin_index = find_OP_bin_index(OP, muca.OPBinLimits);
        // Only increment if on the min-max interval
        if (bin_index >= 0) {
            muca.N_OP_BinTotal[bin_index] += 1;
            muca.WValues[bin_index] += muca.WParam.DIP.C;
        }

        muca.weightIterationCount += 1;


        if (muca.weightIterationCount % muca.WParam.DIP.single_check_interval == 0) {

            double base = *std::min_element(muca.WValues.begin(), muca.WValues.end());
            for (int m = 0; m < N; ++m) {
                // Always set minimum weight to zero. This is inconsequential
                // as only the differences matter.
                muca.WValues[m] -= base;
                muca.N_OP_Bin[m] = 0;
            }

            // Visuals
            if (muca.WParam.visuals)
                muca.print_iteration_histogram();

            // If condition satisfied, zero the totals and decrease C
            if (muca.finish_check(muca)) {
                for (int m = 0; m < N; m++) {
                    muca.N_OP_BinTotal[m] = 0;
                }

                muca.WParam.DIP.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = "
                           << muca.WParam.DIP.C << "\n";
            }

            continue_iteration = true;
            if (muca.WParam.DIP.C < muca.WParam.DIP.C_min) {
                hila::out0 << "Muca: Reached minimum update size C = " << muca.WParam.DIP.C
                           << " Weight iteration complete.\n";
                continue_iteration = false;
            }

            muca.write_weight_function("intermediate_weight.dat");
        }
    }
    hila::broadcast(continue_iteration);
    return continue_iteration;
}

bool iterate_weight_function_direct_smooth(Muca &muca, double OP) {
    bool continue_iteration;
    if (hila::myrank() == 0) {
        const int samples = muca.WParam.DIP.sample_size;
        const int N = muca.WValues.size();

        // Don't bin visits outside of the binned areas
        if (OP > muca.OPBinLimits.front() && OP < muca.OPBinLimits.back()) {
            auto it = std::lower_bound(muca.OPBinLimits.begin(), muca.OPBinLimits.end(), OP);
            int i = std::distance(muca.OPBinLimits.begin(), it) - 1;

            // clang-format off
            // Helps smooth the weight function (kind of kernel density estimation)
            if(i>=2)  muca.N_OP_Bin[i-2] += 1;//#
            if(i>=1)  muca.N_OP_Bin[i-1] += 3;//###
                      muca.N_OP_Bin[i+0] += 5;//##### hit
            if(i<N-1) muca.N_OP_Bin[i+1] += 3;//###
            if(i<N-2) muca.N_OP_Bin[i+2] += 1;//#
            // clang-format on
        }

        muca.weightIterationCount += 1;

        if (muca.weightIterationCount >= samples) {
            for (int m = 0; m < N; m++) {
                muca.WValues[m] += muca.WParam.DIP.C * (muca.N_OP_Bin[m] - muca.N_OP_Bin[0]) / N;
                muca.N_OP_BinTotal[m] += muca.N_OP_Bin[m];
            }

            if (muca.WParam.visuals)
                muca.print_iteration_histogram();

            for (int m = 0; m < N; ++m) {
                muca.N_OP_Bin[m] = 0;
            }
            muca.weightIterationCount = 0;

            if (muca.finish_check(muca)) {
                for (int m = 0; m < N; m++) {
                    muca.N_OP_BinTotal[m] = 0;
                }

                muca.WParam.DIP.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = "
                           << muca.WParam.DIP.C << "\n";
            }
            muca.write_weight_function("intermediate_weight.dat");
        }

        continue_iteration = true;
        if (muca.WParam.DIP.C < muca.WParam.DIP.C_min) {
            hila::out0 << "Muca: Reached minimum update size C = " << muca.WParam.DIP.C
                       << " Weight iteration complete.\n";
            continue_iteration = false;
        }
    }
    hila::broadcast(continue_iteration);
    return continue_iteration;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Prints out a crude horisontal histogram.
/// @details Procures a crude horisontal ASCII histogram based on the N_OP_Bin
///          vector. The histogram bin heights are proportional to N_OP_Bin
///          values. This is not very expressive for large N, as it won't fit
///          the height of the screen.
////////////////////////////////////////////////////////////////////////////////
inline void Muca::print_iteration_histogram() const {
    int samples = WParam.DIP.sample_size;
    int N = WValues.size();
    // Find maximum bin content for normalisation
    int nmax = *std::max_element(N_OP_Bin.begin(), N_OP_Bin.end());
    // Write a column header
    printf("Order Parameter     Weight 	         Number of hits\n");

    for (int m = 0; m < N; ++m) {
        // For each bin get a number of "|":s proportional to the number of
        // hits to each bin and print it out along with relevant numerical
        // values
        std::string n_sum_hist = "";
        if (N_OP_BinTotal[m] > 0)
            n_sum_hist += "O";
        for (int i = 0; i < int(N_OP_Bin[m] * 200.0 / samples); i++) {
            n_sum_hist += "|";
        }
        printf("%-20.3f%-20.3f%d\t\t\t%s\n", OPValues[m], WValues[m], N_OP_Bin[m],
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
inline void Muca::setup_iteration() {
    // Initialise iterate_weights by pointing it at the
    // correct method
    if (WParam.method.compare("direct") == 0) {
        if (WParam.DIP.sample_size > 1) {
            iterate_weights = &iterate_weight_function_direct;
        } else {
            iterate_weights = &iterate_weight_function_direct_single;
        }
        WParam.DIP.C = WParam.DIP.C_init;
    } else {
        iterate_weights = &iterate_weight_function_direct;
        WParam.DIP.C = WParam.DIP.C_init;
        if (hila::myrank() == 0) {
            printf("Muca: note: input iteration method `%s` did not match any method:\n"
                   "            setting the default `direct` method\n",
                   WParam.method.c_str());
        }
    }

    // Zero the iteration counter
    weightIterationCount = 0;

    // Set up the finish condition pointer for the direct method.
    if (WParam.DIP.finish_condition.compare("all_visited") == 0) {
        finish_check = &all_visited;
    } else if (WParam.DIP.finish_condition.compare("ends_visited") == 0) {
        finish_check = &first_last_visited;
    } else {
        finish_check = &all_visited;
        if (hila::myrank() == 0) {
            printf("Muca: note: input finish condition `%s` did not match any:\n"
                   "            setting the default `all_visited` condition\n",
                   WParam.DIP.finish_condition.c_str());
        }
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
void Muca::set_continuous_iteration(bool YN) {
    if (hila::myrank() == 0)
        WParam.AR_iteration = YN;
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
bool Muca::initialise(const std::string wfile_name) {
    // Read parameters into g_WParam struct
    read_weight_parameters(wfile_name);
    bool ret = true;
    // This is fine to do just for process zer0
    if (hila::myrank() == 0) {
        // Read pre-existing weight if given
        if (WParam.weight_loc.compare("NONE") != 0) {
            if (!read_weight_function(WParam.weight_loc)) {
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
