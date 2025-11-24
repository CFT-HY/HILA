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
// namespace muca {

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
std::string append_time_stamp(const std::string &outfile_name_base) {
    std::string filename = outfile_name_base;

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
    WParam.weight_loc        = par.get("weight_file_location");
    WParam.method            = par.get("iteration_method");
    std::string hard_walls   = par.get("hard_walls");
    WParam.hard_walls        = (hard_walls.compare("YES") == 0);
    WParam.max_OP            = par.get("max_OP");
    WParam.min_OP            = par.get("min_OP");
    WParam.bin_number        = par.get("number_of_bins");
    std::string iter_vis     = par.get("iteration_visuals");
    WParam.visuals           = (iter_vis.compare("YES") == 0);
    WParam.AR_iteration = false;

    // TODO ? check that these parameters are in allowed range ?
    // Direct iteration parameters
    DI.finish_condition      = par.get("DI_finish_condition");      // all_visited
    DI.sample_size           = par.get("DI_sample_size");           // 100
    DI.single_check_interval = par.get("DI_single_check_interval"); // 100
    DI.C_init                = par.get("DI_C_init");                // 0.1
    DI.C_min                 = par.get("DI_C_min");                 // 0.001
    DI.C = DI.C_init; // set also the initial value here

    // Canonical iteration parameters
    CI.sample_size      = par.get("CI_sample_size");          // 100
    CI.initial_bin_hits = par.get("CI_initial_bin_hits");     // 5
    CI.min_bin_hits     = par.get("CI_min_bin_hits");         // 8
    CI.max_iters        = par.get("CI_max_iters");            // 400
    CI.OC_factor        = par.get("CI_overcorrect_factor");   // 2.0
    CI.OC_max_iter      = par.get("CI_overcorrect_max_iter"); //
    CI.OC_frequency     = par.get("CI_overcorrect_frequency");//

    par.close();
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
    OP_bin_limits = std::vector<double>(data_length);
    OP_bin_centers = std::vector<double>(data_length - 1);
    weights = std::vector<double>(data_length - 1);

    // Read in the values. Note that OP_bin_limits has one more entry than
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

            OP_bin_limits[count] = bin_lim;
            OP_bin_centers[count] = centre;
            weights[count] = weight;
        }
        // And the rightmost bin limit.
        if (count == data_length - 1) {
            if (1 != sscanf(line.c_str(), "%e", &bin_lim)) {
                printf("ERROR: Could not match `%%e` at %s:%d \n", W_function_filename.c_str(),
                       line_no);
                return false;
            }

            OP_bin_limits[count] = bin_lim;
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
        OP_bin_centers[i] = (OP_bin_limits[i + 1] + OP_bin_limits[i]) / 2.0;
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
///          TODO: Add std::string input that can contain user specified header data.
///
/// @param W_function_filename
/// @return false if writing fails
////////////////////////////////////////////////////////////////////////////////
bool Muca::write_weight_function(const std::string &W_function_filename) {

    std::ofstream W_file;

    W_file.open(W_function_filename.c_str());
    if (!W_file.is_open()) {
        printf("WARNING: Could not open file `%s` for `write_weight_function()`: %s\n",
               W_function_filename.c_str(), std::strerror(errno));
        return false;
    }
    // write_weight_file_header(W_file, FP, RP, g_WParam);

    to_file(W_file, "%s", "OP_bin_limit\tOP_value\tWeight\n");
    int i;
    for (i = 0; i < OP_bin_centers.size(); ++i) {
        to_file(W_file, "%e\t%e\t%e\n", OP_bin_limits[i], OP_bin_centers[i], weights[i]);
    }
    // Remember to write the last bin upper limit
    to_file(W_file, "%e\n", OP_bin_limits[i]);
    W_file.close();
    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Returns a weight associated to the used order parameter.
/// @details The function uses supplied pairs of points to linearly interpolate
///          the function on the interval. This interpolant provides the
///          requested weights to be used as the multicanonical weight.
///
///          The bin centers are used to interpolate the weight function.
///          The last bins are interpolated as constants (towards the range ends) after the bin
///          center.
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
    if (OP <= OP_bin_centers.front()) {
        if (WParam.hard_walls)
            val = std::numeric_limits<double>::infinity();
        else
            val = weights.front();
    } else if (OP >= OP_bin_centers.back()) {
        if (WParam.hard_walls)
            val = std::numeric_limits<double>::infinity();
        else
            val = weights.back();
    }
    // Otherwise find interval, calculate slope, base index, etc.
    // Basic linear interpolation to obtain the weight value.
    else {
        auto it = std::lower_bound(OP_bin_centers.begin(), OP_bin_centers.end(), OP);
        // `it` is the first bin_center for which *it >= OP
        // note: `it` is never .begin() or .end()-1 since .front() < OP <.back(), checked above
        //       thus: 0 <= i <= .size() - 2, so this is ok.
        int i = std::distance(OP_bin_centers.begin(), it) - 1;
        double dy = weights[i + 1] - weights[i];
        double dx = OP_bin_centers[i + 1] - OP_bin_centers[i];
        double slope = dy / dx;

        double xbase = OP_bin_centers[i];
        double ybase = weights[i];

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
    constexpr int update = 0;
    constexpr int AR_iterate = 1;
    // flags[0] = accept/reject, flags[1] = do iterate weights,
    std::array<bool, 2> flags;
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
            flags[update] = true;
        } else {
            flags[update] = false;
        }

        // Get value from process 0
        flags[AR_iterate] = WParam.AR_iteration;
    }

    // Broadcast the update status to other processes along with the
    // weight iteration parameter
    hila::broadcast(flags);

    // Check if iteration is enabled
    if (flags[AR_iterate]) {
        bool continue_iter = iterate_weights(*this, flags[update] ? OP_new : OP_old);
        set_weight_iter_flag(continue_iter);
    }

    return flags[update];
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Finds the index of the correc order parameter bin.
/// @details Using the bin limit vector the correct order parameter bin is
///          determined through a simple standard library search of
///          g_OP_bin_limits. When the value of the given order parameter is
///          out of range, -1 is returned.
///
/// @param OP   value of the order parameter
/// @param OP_bin_limits   bin limits of the weight function hist
/// @return integer index for the vector OP_bin_hits
////////////////////////////////////////////////////////////////////////////////
static inline int find_OP_bin_index(double OP, const std::vector<double> &OP_bin_limits) {
    // Return -1 when not in the interval
    if (OP <= OP_bin_limits.front()) {
        return -1;
    } else if (OP >= OP_bin_limits.back()) {
        return -1;
    }
    // Find index of minimum edge value such that edge < OP:
    auto it = std::lower_bound(OP_bin_limits.begin(), OP_bin_limits.end(), OP);
    int lower_limit = std::distance(OP_bin_limits.begin(), it) - 1;
    return lower_limit;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Same as find_OP_bin_index, except uses the index to simply modify the
///        bin hit vector OP_bin_hits. Does not modify when outside of the range.
/// @details
///
/// @param OP   value of the order parameter
////////////////////////////////////////////////////////////////////////////////
inline void Muca::bin_hit_OP_value(double OP) {
    // Don't bin visits outside of the binned areas
    if (OP <= OP_bin_limits.front()) {
        return;
    } else if (OP >= OP_bin_limits.back()) {
        return;
    }
    // Find index of minimum edge value such that edge < OP:
    auto it = std::lower_bound(OP_bin_limits.begin(), OP_bin_limits.end(), OP);
    int lower_limit = std::distance(OP_bin_limits.begin(), it) - 1;
    OP_bin_hits[lower_limit] += 1;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Checks if all the bins have been visited.
/// @details Simply checks whether all bins have a nonzero number of entries in
/// `muca.OP_bin_hits_total`.
///
/// @param  muca the multicanonical context
/// @return a boolean indicating the statement
////////////////////////////////////////////////////////////////////////////////
bool all_visited(const Muca &muca) {
    const auto &n = muca.OP_bin_hits_total;
    int len = n.size();
    for (int i = 0; i < len; ++i) {
        if (n[i] == 0)
            return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Checks if the first and last bin have been visited
/// @details Simply checks whether first and last entry of `muca.OP_bin_hits_total` are nonzero.
///
/// @param  muca the multicanonical context
/// @return a boolean indicating the statement
////////////////////////////////////////////////////////////////////////////////
bool first_last_visited(const Muca &muca) {
    const auto &n = muca.OP_bin_hits_total;
    return !((n.front() == 0) or (n.back() == 0));
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
    OP_bin_limits = get_equidistant_bin_limits(WParam.min_OP, WParam.max_OP, WParam.bin_number);
    for (int i = 0; i < OP_bin_centers.size(); i++) {
        double centre = (OP_bin_limits[i + 1] + OP_bin_limits[i]) / 2.0;
        OP_bin_centers[i] = centre;
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
        weights = std::vector<double>(N, 0.0);
        OP_bin_centers = std::vector<double>(N, 0.0);
        OP_bin_limits = std::vector<double>(N + 1, 0.0);

        setup_equidistant_bins();

        OP_bin_hits = std::vector<int>(N, 0);
        OP_bin_hits_total = std::vector<int>(N, 0);
    }
    // Same for predetermined bins. OP_bin_centers, weights
    // and OP_BinLimits have been read from the input file.
    // To prevent any accidents with the iterators, these bin vectors
    // are initialised in all cases. Should really not affect performance.
    else {
        int N = weights.size();
        OP_bin_hits = std::vector<int>(N, 0);
        OP_bin_hits_total = std::vector<int>(N, 0);
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
        int samples = muca.DI.sample_size;
        int N = muca.weights.size();

        muca.bin_hit_OP_value(OP);
        muca.weightIterationCount += 1;

        if (muca.weightIterationCount >= samples) {
            for (int m = 0; m < N; m++) {
                muca.weights[m] += muca.DI.C * muca.OP_bin_hits[m] * N / samples;
                muca.OP_bin_hits_total[m] += muca.OP_bin_hits[m];
            }

            if (muca.WParam.visuals)
                muca.print_iteration_histogram();

            double base = *std::min_element(muca.weights.begin(), muca.weights.end());
            for (int m = 0; m < N; ++m) {
                // Always set minimum weight to zero. This is inconsequential
                // as only the differences matter.
                muca.weights[m] -= base;
                muca.OP_bin_hits[m] = 0;
            }
            muca.weightIterationCount = 0;

            if (muca.finish_check(muca)) {
                for (int m = 0; m < N; m++) {
                    muca.OP_bin_hits_total[m] = 0;
                }

                muca.DI.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = " << muca.DI.C
                           << "\n";
            }
            muca.write_weight_function("intermediate_weight.dat");
        }

        continue_iteration = true;
        if (muca.DI.C < muca.DI.C_min) {
            hila::out0 << "Muca: Reached minimum update size C = " << muca.DI.C
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
///          DI.single_check_interval to prevent excessive checking
///          for the visits.
///
/// @param  OP   order parameter of the current configuration (user supplied)
/// @return boolean indicating whether the iteration is considered complete
////////////////////////////////////////////////////////////////////////////////
bool iterate_weight_function_direct_single(Muca &muca, double OP) {
    int continue_iteration;
    if (hila::myrank() == 0) {
        int samples = muca.DI.sample_size;
        int N = muca.weights.size();

        int bin_index = find_OP_bin_index(OP, muca.OP_bin_limits);
        // Only increment if on the min-max interval
        if (bin_index >= 0) {
            muca.OP_bin_hits_total[bin_index] += 1;
            muca.weights[bin_index] += muca.DI.C;
        }

        muca.weightIterationCount += 1;


        if (muca.weightIterationCount % muca.DI.single_check_interval == 0) {

            double base = *std::min_element(muca.weights.begin(), muca.weights.end());
            for (int m = 0; m < N; ++m) {
                // Always set minimum weight to zero. This is inconsequential
                // as only the differences matter.
                muca.weights[m] -= base;
                muca.OP_bin_hits[m] = 0;
            }

            // Visuals
            if (muca.WParam.visuals)
                muca.print_iteration_histogram();

            // If condition satisfied, zero the totals and decrease C
            if (muca.finish_check(muca)) {
                for (int m = 0; m < N; m++) {
                    muca.OP_bin_hits_total[m] = 0;
                }

                muca.DI.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = " << muca.DI.C
                           << "\n";
            }

            continue_iteration = true;
            if (muca.DI.C < muca.DI.C_min) {
                hila::out0 << "Muca: Reached minimum update size C = " << muca.DI.C
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
        const int samples = muca.DI.sample_size;
        const int N = muca.weights.size();

        // Don't bin visits outside of the binned areas
        if (OP > muca.OP_bin_limits.front() && OP < muca.OP_bin_limits.back()) {
            auto it = std::lower_bound(muca.OP_bin_limits.begin(), muca.OP_bin_limits.end(), OP);
            int i = std::distance(muca.OP_bin_limits.begin(), it) - 1;

            // clang-format off
            // Helps smooth the weight function (kind of kernel density estimation)
            if(i>=2)  muca.OP_bin_hits[i-2] += 1;//#
            if(i>=1)  muca.OP_bin_hits[i-1] += 3;//###
                      muca.OP_bin_hits[i+0] += 5;//##### hit
            if(i<N-1) muca.OP_bin_hits[i+1] += 3;//###
            if(i<N-2) muca.OP_bin_hits[i+2] += 1;//#
            // clang-format on
        }

        muca.weightIterationCount += 1;

        if (muca.weightIterationCount >= samples) {
            for (int m = 0; m < N; m++) {
                muca.weights[m] += muca.DI.C * (muca.OP_bin_hits[m] - muca.OP_bin_hits[0]) / N;
                muca.OP_bin_hits_total[m] += muca.OP_bin_hits[m];
            }

            if (muca.WParam.visuals)
                muca.print_iteration_histogram();

            for (int m = 0; m < N; ++m) {
                muca.OP_bin_hits[m] = 0;
            }
            muca.weightIterationCount = 0;

            if (muca.finish_check(muca)) {
                for (int m = 0; m < N; m++) {
                    muca.OP_bin_hits_total[m] = 0;
                }

                muca.DI.C /= 1.5;
                hila::out0 << "Muca: Decreasing update size. New update size C = " << muca.DI.C
                           << "\n";
            }
            muca.write_weight_function("intermediate_weight.dat");
        }

        continue_iteration = true;
        if (muca.DI.C < muca.DI.C_min) {
            hila::out0 << "Muca: Reached minimum update size C = " << muca.DI.C
                       << " Weight iteration complete.\n";
            continue_iteration = false;
        }
    }
    hila::broadcast(continue_iteration);
    return continue_iteration;
}

/// 'Canonical' iteration method based on https://arxiv.org/abs/hep-lat/9804019 p.19
///     TODO: doublecheck how edge bins are handeled
///     TODO: better stop iter criterion
bool iterate_weight_function_canonical(Muca &muca, double OP) {
    bool continue_iteration = true;
    if (hila::myrank() == 0) {

        muca.weightIterationCount += 1;

        // Accumulate bin hits
        if (OP > muca.OP_bin_limits.front() && OP < muca.OP_bin_limits.back()) {
            auto it = std::lower_bound(muca.OP_bin_limits.begin(), muca.OP_bin_limits.end(), OP);
            int i = std::distance(muca.OP_bin_limits.begin(), it) - 1;

            muca.OP_bin_hits[i] += 1;
            double w = muca.weight_function(OP);
            muca.CI.can_hist[i] += ::exp(w);
        }

        // update weights
        if (muca.weightIterationCount >= muca.CI.sample_size) {
            muca.weightIterationCount = 0;

            // before updating weights swap weights to the non-corrected one
            std::swap(muca.CI.noc_weights, muca.weights);

            // normalisation necesary if bins with different lengths
            for (size_t i = 1; i < muca.CI.can_hist.size(); ++i) {
                double bin_width = muca.OP_bin_limits[i] - muca.OP_bin_limits[i - 1];
                muca.CI.can_hist[i - 1] /= bin_width;
            }

            int min_hits = muca.CI.min_bin_hits;
            std::vector<double> w_prev = muca.weights; // copy

            muca.CI.hits_nsum[0] += muca.OP_bin_hits[0];
            for (size_t i = 1; i < muca.CI.can_hist.size(); ++i) {
                muca.CI.hits_nsum[i] += muca.OP_bin_hits[i];

                size_t gi = 0;
                if (muca.OP_bin_hits[i] >= min_hits && muca.OP_bin_hits[i - 1] >= min_hits) {
                    gi = muca.OP_bin_hits[i] + muca.OP_bin_hits[i - 1];
                } else {
                    // no update, but need to still do this if [i-1] was updated..
                    muca.weights[i] = w_prev[i] - w_prev[i - 1] + muca.weights[i - 1];
                    continue;
                }

                // [hep-lat/9804019] eq. (6.12):
                double ln = -(double)gi * ::log(muca.CI.can_hist[i - 1] / muca.CI.can_hist[i]);
                size_t gi_sum = muca.CI.gsum[i] + gi;
                muca.weights[i] = muca.weights[i - 1] +
                                  ((w_prev[i] - w_prev[i - 1]) * muca.CI.gsum[i] + ln) / (gi_sum);
                muca.CI.gsum[i] = gi_sum;
            }

            // reset counters
            for (size_t i = 0; i < muca.OP_bin_hits.size(); ++i) {
                muca.OP_bin_hits[i] = 0;
                muca.CI.can_hist[i] = 0;
            }

            // Overcorrect [hep-lat/9804019] eq. (6.13)
            // Overcorrection checks which bins have been visited the most and
            // modifies the weights in such a manner that the frequently visited
            // bins are weighted (disproportionally) heavily. This should
            // encourage the process to visit other bins during the next run of
            // iterations.
            // overcorrected weights are used only for the next updates.
            // The weight function to be modified after the next updates will be the non
            // overcorrected one. So we save the wights before overcorrection.
            std::memcpy(muca.CI.noc_weights.data(), muca.weights.data(),
                        muca.CI.noc_weights.size() * sizeof(muca.CI.noc_weights[0]));
            double C = muca.CI.OC_factor;
            double d0 = muca.OP_bin_limits[1] - muca.OP_bin_limits[0];
            for (size_t i = 1; i < muca.weights.size(); ++i) {
                double bin_width = muca.OP_bin_limits[i] - muca.OP_bin_limits[i - 1];
                muca.weights[i - 1] +=
                    C * ::log((d0 / bin_width) * muca.CI.hits_nsum[i - 1] / muca.CI.hits_nsum[0]);
            }

            muca.write_weight_function("intermediate_weight.dat");

            // check if we should stop
            muca.CI.weight_update_count++;
            if (muca.CI.weight_update_count > muca.CI.max_iters)
                continue_iteration = false;
            // print progress
            if (muca.CI.weight_update_count % (muca.CI.max_iters / 100) == 0)
                printf("Muca: iteration %d / %d\n", muca.CI.weight_update_count, muca.CI.max_iters);
        }
    }
    hila::broadcast(continue_iteration);
    return continue_iteration;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Prints out a crude horisontal histogram.
/// @details Procures a crude horisontal ASCII histogram based on the OP_bin_hits
///          vector. The histogram bin heights are proportional to OP_bin_hits
///          values. This is not very expressive for large N, as it won't fit
///          the height of the screen.
////////////////////////////////////////////////////////////////////////////////
inline void Muca::print_iteration_histogram() const {
    int samples = DI.sample_size;
    int N = weights.size();
    // Find maximum bin content for normalisation
    int nmax = *std::max_element(OP_bin_hits.begin(), OP_bin_hits.end());
    // Write a column header
    printf("Order Parameter     Weight 	         Number of hits\n");

    for (int m = 0; m < N; ++m) {
        // For each bin get a number of "|":s proportional to the number of
        // hits to each bin and print it out along with relevant numerical
        // values
        std::string n_sum_hist = "";
        if (OP_bin_hits_total[m] > 0)
            n_sum_hist += "O";
        for (int i = 0; i < int(OP_bin_hits[m] * 200.0 / samples); i++) {
            n_sum_hist += "|";
        }
        printf("%-20.3f%-20.3f%d\t\t\t%s\n", OP_bin_centers[m], weights[m], OP_bin_hits[m],
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
        if (DI.sample_size > 1) {
            iterate_weights = &iterate_weight_function_direct;
        } else {
            iterate_weights = &iterate_weight_function_direct_single;
        }
        DI.C = DI.C_init;
    } else if (WParam.method.compare("direct_smooth") == 0) {
        iterate_weights = &iterate_weight_function_direct_smooth;
        DI.C = DI.C_init;
    } else if (WParam.method.compare("canonical") == 0) {
        iterate_weights = &iterate_weight_function_canonical;

        // initialise stuff specific for canonical iteration
        constexpr int n_initial = 1;
        CI.weight_update_count = 0;
        CI.hits_nsum = std::vector<int>(OP_bin_hits.size(), n_initial);
        CI.gsum = std::vector<int>(OP_bin_hits.size(), CI.initial_bin_hits);
        CI.can_hist = std::vector<double>(OP_bin_hits.size(), 0);
        CI.noc_weights = weights; // copy

    } else {
        iterate_weights = &iterate_weight_function_direct;
        DI.C = DI.C_init;
        if (hila::myrank() == 0) {
            printf("Muca: note: input iteration method `%s` did not match any method:\n"
                   "            setting the default `direct` method\n",
                   WParam.method.c_str());
        }
    }

    // Zero the iteration counter
    weightIterationCount = 0;

    // Set up the finish condition pointer for the direct method.
    if (DI.finish_condition.compare("all_visited") == 0) {
        finish_check = &all_visited;
    } else if (DI.finish_condition.compare("ends_visited") == 0) {
        finish_check = &first_last_visited;
    } else {
        finish_check = &all_visited;
        if (hila::myrank() == 0) {
            printf("Muca: note: input finish condition `%s` did not match any:\n"
                   "            setting the default `all_visited` condition\n",
                   DI.finish_condition.c_str());
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
/// @details Sets up iteration variables. NOTE: these are set up only on `hila::myrank()==0`.
///          Can be called multiple times and must be called at least once before attempting to use
///          any of the muca methods.
///
/// @param wfile_name   path to the weight parameter file
/// @return false if something fails during init.
////////////////////////////////////////////////////////////////////////////////
bool Muca::initialise(const std::string wfile_name) {
    // Read parameters into WParam struct
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

//} // namespace muca
} // namespace hila
