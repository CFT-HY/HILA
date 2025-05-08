
#include <time.h>
#include <chrono>

#include "defs.h"

//////////////////////////////////////////////////////////////////
// Time related routines (runtime - timing - timelimit)
// Check timing.h for details
//////////////////////////////////////////////////////////////////

#include "com_mpi.h"

// these includes need to be outside namespace hila
#include <csignal>
#include <cstring>

namespace hila {

// This stores the start time of the program
static double start_time = -1.0;

/////////////////////////////////////////////////////////////////
/// Timer routines - for high-resolution event timing
/////////////////////////////////////////////////////////////////

// store all timers in use
std::vector<timer *> timer_list = {};

// initialize timer to this timepoint
void timer::init(const char *tag) {
    if (tag != nullptr)
        label = tag;
    reset();
    // Store it on 1st use!
}

// remove the timer also from the list
void timer::remove() {
    for (auto it = timer_list.begin(); it != timer_list.end(); ++it) {
        if (*it == this) {
            timer_list.erase(it);
            return;
        }
    }
}

void timer::reset() {
    t_start = t_total = 0.0;
    count = 0;
    is_on = is_error = false;
}

void timer::error() {
    if (!is_error) {
        hila::out0 << " **** Timer '" << label
                   << "' error, unbalanced start/stop.  Removing from statistics\n";
    }
    is_error = true;
}

double timer::start() {
    if (is_on)
        error();
    is_on = true;

    // Move storing the timer ptr here, because if timer is initialized
    // in the global scope the timer_list is possibly initialized later!
    if (count == 0) {
        timer_list.push_back(this);
    }

#ifdef GPU_SYNCHRONIZE_TIMERS
    gpuStreamSynchronize(0);
#endif

    t_start = hila::gettime();
    return t_start;
}

double timer::stop() {
    if (!is_on)
        error();
    is_on = false;

#ifdef GPU_SYNCHRONIZE_TIMERS
    gpuStreamSynchronize(0);
#endif

    double e = hila::gettime();
    t_total += (e - t_start);
    count++;
    return e;
}

timer_value timer::value() {
    timer_value r;
    r.time = t_total;
    r.count = count;
    return r;
}

void timer::report(bool print_not_timed) {
    if (hila::myrank() == 0) {
        char line[202];

        // time used during the counter activity
        double ttime = gettime();
        if (count > 0 && !is_error) {
            if (t_total / count > 0.1) {
                std::snprintf(line, 200, "%-20s: %14.5f %14ld %12.5f s  %9.6f\n", label.c_str(),
                              t_total, (long)count, t_total / count, t_total / ttime);
            } else if (t_total / count > 1e-4) {
                std::snprintf(line, 200, "%-20s: %14.5f %14ld %12.5f ms %9.6f\n", label.c_str(),
                              t_total, (long)count, 1e3 * t_total / count, t_total / ttime);
            } else {
                std::snprintf(line, 200, "%-20s: %14.5f %14ld %12.5f us %9.6f\n", label.c_str(),
                              t_total, (long)count, 1e6 * t_total / count, t_total / ttime);
            }
            hila::out << line;
        } else if (!is_error && print_not_timed) {
            std::snprintf(line, 200, "%-20s: no timed calls made\n", label.c_str());
            hila::out << line;
        } else if (is_error) {
            std::snprintf(line, 200, "%-20s: error:unbalanced start/stop\n", label.c_str());
            hila::out << line;
        }
    }
}

void report_timers() {
    if (hila::myrank() == 0) {
        if (timer_list.size() > 0) {

#if defined(CUDA) || defined(HIP)
#if defined(GPU_SYNCHRONIZE_TIMERS)
            hila::out << "TIMERS: synchronized to GPU kernel execution (GPU_SYNCHRONIZE_TIMERS "
                         "defined)\n";
#else
            hila::out << "TIMERS: GPU_SYNCHRONIZE_TIMERS not defined, fine-grained timing "
                         "likely to be incorrect\n";
#endif
#endif

            hila::out << "TIMER REPORT:             total(sec)          calls       "
                         "time/call  fraction\n";
            hila::out << "------------------------------------------------------------"
                         "-----------------\n";

            for (auto tp : timer_list) {
                tp->report();
            }

            hila::out << "------------------------------------------------------------"
                         "-----------------\n";
        } else {
            hila::out << "No timers defined\n";
        }
    }
}

/////////////////////////////////////////////////////////////////
/// Use clock_gettime() to get the accurate time
/// (alternative: use gettimeofday()  or MPI_Wtime())
/// gettime returns the time in secs since program start

double gettime() {
    struct timespec tp;

    if (start_time == -1.0)
        inittime();

    clock_gettime(CLOCK_MONOTONIC, &tp);
    return (((double)tp.tv_sec - start_time) + 1.0e-9 * (double)tp.tv_nsec);
}

void inittime() {
    if (start_time == -1.0) {
        start_time = 0.0;
        start_time = gettime();
    }
}

//////////////////////////////////////////////////////////////////
/// Routines for checking remaining cpu-time
/// void setup_timelimit():   set the time limit to watch
/// bool time_to_finish(): is called periodically on a point where exit can be done.
/// It uses the max of the time intervals for the estimate for one further round.
/// If not enough time returns true, else false
///

static double timelimit = 0;

///
/// Setup time limit with seconds

void setup_timelimit(const double secs) {
    timelimit = secs;
    hila::broadcast(timelimit);
    hila::out0 << "Time limit is " << timelimit << " seconds\n";
}

///
/// setup time limit from the time given in timestr
/// Format is d-h:m:s, not required to be "normalized" to std ranges
/// Fields can be unused from the largest fields onwards, i.e. simplest case only seconds
/// string "slurm" indicates that we call slurm 'squeue' to obtain the time limit

void setup_timelimit(const std::string &timestr) {
    //
    constexpr int timelimit_buf_size = 100;

    int status = 0;
    if (hila::myrank() == 0) {
        const char *str = timestr.c_str();
        char buf[timelimit_buf_size];

        if (timestr == "slurm") {

            const char cmd[] = "squeue -h --job ${SLURM_JOB_ID} -O TimeLeft";
            std::FILE *fp = popen(cmd, "r");

            if (fp && fgets(buf, timelimit_buf_size - 1, fp)) {
                buf[timelimit_buf_size - 1] = 0;
                // zero extra spaces and lf at the end of the buf
                for (int i = std::strlen(buf) - 1; i >= 0 && std::isspace(buf[i]); i--)
                    buf[i] = 0;

                str = buf;
                hila::out0 << "Got time limit with command '" << cmd << '\n';
            } else {
                hila::out0 << "COULD NOT GET TIME FROM squeue COMMAND\n";
                status = -1; // exit the program
            }
            pclose(fp);
        }

        if (status == 0) {

            unsigned d{0}, h{0}, m{0}, s{0};
            // use short circuiting of || here to stop parsing on 1st match
            // zeroing the incorrectly read time variables in the same chain
            int nargs = 5;
            if (std::sscanf(str, "%u-%u:%u:%u", &d, &h, &m, &s) == --nargs ||
                std::sscanf(str, "%u:%u:%u", &h, &m, &s) == --nargs ||
                std::sscanf(str, "%u:%u", &m, &s) == --nargs ||
                std::sscanf(str, "%u", &s) == --nargs) {

                if (nargs < 4) d = 0;
                if (nargs < 3) h = 0;
                if (nargs < 2) m = 0;

                timelimit = s + 60.0 * (m + 60.0 * (h + 24.0 * d));
                hila::out0 << "Time limit is " << str << " = " << timelimit << " seconds\n";

            } else {
                hila::out0 << "INVALID TIMELIMIT -t ARGUMENT " << str << '\n';
                status = -1; // exit the program
            }
        }
    }
    hila::broadcast(status);
    if (status == -1)
        hila::terminate(0);

    hila::broadcast(timelimit);
}


//////////////////////////////////////////////////////////////////////////////////////////////////
/// Set up signal handling - signal SIGUSR1 causes the time_to_finish function to return
/// true.  Function hila::setup_signal_handler(); is called at the beginning of the program.
/// Function hila::signal_status() returns signal value (!=0) if signal has been received.
/// Signal is _not_ broadcast across nodes
///
/// Could also trap SIGTERM, but that would mean programs which do not handle the signal
/// would not be killable by ctrl-C
///

static volatile std::sig_atomic_t received_signal;

void signal_handler(int signal) {
    received_signal = signal;
}

void setup_signal_handler() {
    // std::signal(SIGTERM, signal_handler);
    std::signal(SIGUSR1, signal_handler);
}

int signal_status() {
    return received_signal;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Check the cpu time limit or if signal SIGUSR1 is up - this function is meant to be called
/// periodically, and returns true if it is time to exit or signal has been raised.
///
/// Use case: on main loop check 'hila::time_to_finish()' periodically, and if it returns true
/// checkpoint/exit.
///
/// This makes sense only if the program can checkpoint or otherwise clean up.
///
/// Time limit or signaling are alternative ways to ensure clean exit. Both can be used
/// at the same time. Time limit has the advantage that if the program does periodic
/// checkpointing, it automatically adjusts the grace time to allow for checkpointing
/// at the end. With signal the grace time must be estimated in advance.
///
/// Time limit can be given with program command line argument
///   -t <time>  or  -t slurm
/// or calling function  hila::setup_timelimit(time), see above
///
/// Signal can be set in slurm submit script with
///   #SBATCH --signal=SIGUSR1@180
/// where the last number is the time in seconds when the signal SIGUSR1 is sent before
/// the run time expires.  This time must allow for the periodic check interval and the time
/// the cleanup takes.
///
/// Signal can also be sent to slurm jobs from terminal session with
///   $ scancel --signal=SIGUSR1  <jobid>
/// and, of course, for normal "terminal" runs with
///   $ kill -s SIGUSR1  <pid0>
/// Note: signal must be sent to MPI rank 0. It can be sent to all ranks too.
///


bool time_to_finish() {
    static double max_interval = 0.0;
    static double previous_time = 0.0;
    bool finish;

    // is signal up?
    int signal = signal_status();

    if (hila::myrank() == 0) {
        if (signal != 0) {
            finish = true;
            hila::out0 << "FINISH UP ON SIGNAL SIGUSR1\n";

        } else if (timelimit == 0.0) {
            // no signal nor time limit set
            finish = false;

        } else {

            double this_time = gettime();
            if (this_time - previous_time > max_interval)
                max_interval = this_time - previous_time;
            previous_time = this_time;

            // Give 2 min margin for the exit - perhaps needed for writing etc.
            if (timelimit - this_time < max_interval + 2 * 60.0)
                finish = true;
            else
                finish = false;

            // hila::out << "TIMECHECK: " << this_time << "s used, " << timelimit - this_time
            //           << "s remaining\n";

            if (finish)
                hila::out << "CPU TIME LIMIT, EXITING THE PROGRAM\n";
        }
    }
    hila::broadcast(finish);
    return finish;
}

/*****************************************************
 * Time stamp
 */

void timestamp(const char *msg) {
    if (hila::myrank() == 0) {
        int p = hila::out0.precision();
        std::time_t ct = std::time(NULL);
        if (msg != NULL)
            hila::out << msg;
        std::string d = ctime(&ct);
        d.resize(d.size() - 1); // take away \n at the end
        hila::out0 << " -- date " << d << "  run time " << std::setprecision(4) << hila::gettime()
                   << "s" << std::endl;
        hila::out0.precision(p);
        hila::out0.flush();
    }
}

void timestamp(const std::string &msg) {
    hila::timestamp(msg.c_str());
}

} // namespace hila
