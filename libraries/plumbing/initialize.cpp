
#include <cstring>
#include "hila.h"

// define these global var here - somehow NULL needed for ostream
std::ostream hila::out(NULL);
std::ostream hila::out0(NULL);
std::ofstream hila::output_file;
bool hila::about_to_finish = false;
bool hila::check_input = false;
int hila::check_with_nodes;
const char *hila::input_file;
logger_class hila::log;


void vector_type_info();

#include <limits.h>
#include <errno.h>

///////////////////////////////////////////////////////////////////////////////
/// very simple cmdline arg interpreter
///   bool          cmdline.get_option("-opt");   - true if '-opt' is present
///   const char *  cmdline.get_cstring("-flag"); - returns char * to string
///                                                 following -flag
///   int           cmdline.get_int("-I");        - returns (long) int after -I
///   int           cmdline.get_onoff("-t");      - returns 1 for '-t on',
///                                                 -1 if off and 0 if does not appear
/// Args are removed after reading them
///
///////////////////////////////////////////////////////////////////////////////
class cmdlineargs {
  private:
    int argc;
    const char **argv;

  public:
    cmdlineargs(int argc0, char **argv0) {
        argc = argc0;
        argv = (const char **)malloc(argc * sizeof(const char *));
        for (int i = 0; i < argc; i++)
            argv[i] = argv0[i];
    }

    ~cmdlineargs() {
        free(argv);
    }

    bool get_option(const char *flag) {
        int flaglen = strlen(flag);
        assert(flaglen > 0);

        for (int i = 1; i < argc; i++) {
            const char *p = argv[i];

            if (std::strcmp(p, flag) == 0) {
                argc--;
                for (; i < argc; i++)
                    argv[i] = argv[i + 1];
                return (true);
            }
        }

        return false;
    }

    const char *get_cstring(const char *flag) {
        int flaglen = strlen(flag);
        assert(flaglen > 0);

        for (int i = 1; i < argc; i++) {
            const char *p = argv[i];

            // OK if p starts with flag
            if (std::strcmp(p, flag) == 0) {
                if (i > argc - 2) {
                    hila::out0 << "Expecting an argument after command line parameter '" << flag
                               << "'\n";
                    hila::terminate(0);
                }
                p = argv[i + 1];
                argc -= 2;
                for (; i < argc; i++)
                    argv[i] = argv[i + 2];
                return (p);
            }
        }
        return (nullptr);
    }

    long get_int(const char *flag) {
        const char *p = get_cstring(flag);
        char *end;

        if (p == nullptr)
            return LONG_MAX; // not found

        long val = strtol(p, &end, 10);
        if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN)) || (errno != 0 && val == 0) ||
            end == p || *end != 0) {
            hila::out0 << "Expect a number (integer) after command line parameter '" << flag
                       << "'\n";
            hila::terminate(0);
        }
        return val;
    }

    /// returns 1=on, -1=off, 0 not found
    int get_onoff(const char *flag) {
        const char *p = get_cstring(flag);
        if (p == nullptr)
            return 0;
        if (std::strcmp(p, "on") == 0)
            return 1;
        if (std::strcmp(p, "off") == 0)
            return -1;
        hila::out0 << "Command line argument " << flag << " requires value on/off\n";
        hila::terminate(0);
        return 0; // gets rid of a warning of no return value
    }

    int items() {
        return argc - 1; // don't count argv[0]
    }

    void error_if_args_remain() {
        if (argc < 2)
            return;
        if (hila::myrank() == 0) {
            hila::out << "Unknown command line argument:\n";
            for (int i = 1; i < argc; i++) {
                hila::out << "    " << argv[i] << '\n';
            }
            // clang-format off
            hila::out
                << "Recognized:\n"
                << "  -t <seconds>    : cpu time limit\n"
                << "  -o <name>       : output filename (default: stdout)\n"
                << "  -i <name>       : input filename (overrides the 1st hila::input() name)\n"
                << "                    use '-i -' for standard input\n"
                << "  -check          : check input & layout with <nodes>-nodes & exit\n"
                << "                    only with 1 real MPI node (without mpirun)\n"
                << "  -n nodes        : number of nodes used in layout check, only relevant with -check\n"
                << "  -partitions n   : number of partitioned lattice streams\n"
                << "  -sync on/off    : synchronize partition runs (default=no)\n";
            // clang-format on
        }
        hila::terminate(0);
    }
};

/////////////////////////////////////////////////////////////////////////////////
/// Initial setup routines
/////////////////////////////////////////////////////////////////////////////////

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
#include <malloc.h>
#endif

void setup_partitions(cmdlineargs &cl);

void hila::initialize(int argc, char **argv) {

#if (defined(__GNUC__) && !defined(DARWIN) && !defined(_MAC_OSX_)) // || defined(__bg__)
    /* First, adjust malloc so that glibc free() does not
     * release space to the system, increasing the performance
     * of the glib malloc substantially.  The memory use is cyclic,
     * so we can just sit on the max memory.
     */
    // mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
    /* HACK: don't release memory by calling sbrk */
    mallopt(M_TRIM_THRESHOLD, -1);
#endif

    // initialize MPI so that hila::myrank() etc. works
    initialize_communications(argc, &argv);

    // Default output file - we're happy with this unless partitions
    // or otherwise indicated
    // This channels outf to std::cout
    hila::out.rdbuf(std::cout.rdbuf());

    // set the timing so that gettime() returns time from this point
    hila::inittime();


    // open hila::out0 only for node 0
    if (hila::myrank() == 0)
        hila::out0.rdbuf(std::cout.rdbuf());


    // Init command line - after MPI has been started, so
    // that all nodes do this
    cmdlineargs commandline(argc, argv);

    // check the "-check" -input early
    // do it only with 1 node
    if (lattice.nodes.number == 1) {

        if (commandline.get_option("-check")) {
            long nodes = commandline.get_int("-n");
            if (nodes == LONG_MAX)
                nodes = 1;

            hila::check_input = true;
            if (nodes <= 0)
                nodes = 1;
            hila::check_with_nodes = nodes;
            hila::out << "****** INPUT AND LAYOUT CHECK ******" << std::endl;

            // reset node variables
            lattice.mynode.rank = 0;
            lattice.nodes.number = hila::check_with_nodes;
        }
    }

#if defined(CUDA) || defined(HIP)
    if (!hila::check_input) {
        initialize_gpu(lattice.mynode.rank);
    }
#endif

    setup_partitions(commandline);

    // check the output file if partitions not used
    if (hila::partitions.number() == 1) {
        int do_exit = 0;
        if (hila::myrank() == 0) {
            if (const char *name = commandline.get_cstring("-o")) {
                // Open file for append
                if (std::strlen(name) == 0) {
                    hila::out << "Filename must be given with option '-o'\n";
                    do_exit = 1;
                } else if (!hila::check_input) {
                    hila::output_file.open(name, std::ios::out | std::ios::app);
                    if (hila::output_file.fail()) {
                        hila::out << "Cannot open output file " << name << '\n';
                        do_exit = 1;
                    } else {
                        hila::out.flush();
                        hila::out.rdbuf(
                            hila::output_file.rdbuf()); // output now points to output_redirect

                        if (hila::myrank() == 0)
                            hila::out0.rdbuf(hila::out.rdbuf());
                    }
                }
            }
        }
        hila::broadcast(do_exit);
        if (do_exit)
            hila::finishrun();
    }

    if (hila::myrank() == 0) {
        print_dashed_line(u8"HILA ⩩ lattice framework");
        hila::out << "Running program " << argv[0] << "\n";
        hila::out << "with command line arguments '";
        for (int i = 1; i < argc; i++)
            hila::out << argv[i] << ' ';
        hila::out << "'\n";
        hila::out << "Code version: ";
#if defined(GIT_SHA_VALUE)
#define xstr(s) makestr(s)
#define makestr(s) #s
        hila::out << "git SHA " << xstr(GIT_SHA_VALUE) << '\n';
#else
        hila::out << "no git information available\n";
#endif
        hila::out << "Compiled " << __DATE__ << " at " << __TIME__ << '\n';

        hila::out << "with options: EVEN_SITES_FIRST";
#ifndef EVEN_SITES_FIRST
        hila::out << "=0";
#endif
#ifdef SPECIAL_BOUNDARY_CONDITIONS
        hila::out << " SPECIAL_BOUNDARY_CONDITIONS";
#endif
        hila::out << '\n';

        hila::timestamp("Starting");
    }

    long cputime = commandline.get_int("-t");
    if (cputime != LONG_MAX) {
        hila::out0 << "CPU time limit " << cputime << " seconds\n";
        hila::setup_timelimit(cputime);
    } else {
        hila::out0 << "No runtime limit given\n";
    }


    if ((hila::input_file = commandline.get_cstring("-i"))) {
        if (std::strlen(hila::input_file) == 0) {
            hila::out0 << "Filename must be given with '-i <name>'\n"
                       << "Or use '-i stdin' to use standard input\n";
            hila::finishrun();
        }

        hila::out0 << "Input file from command line: " << hila::input_file << '\n';
    }

    // error out if there are more cmdline options
    commandline.error_if_args_remain();

#if defined(OPENMP)
    hila::out0 << "Using option OPENMP - with " << omp_get_max_threads() << " threads\n";
#endif


#if defined(CUDA) || defined(HIP)
#if defined(GPU_AWARE_MPI)
    hila::out0 << "Using GPU_AWARE_MPI\n";
#else
    hila::out0 << "Not using GPU_AWARE_MPI\n";
#endif
    if (!hila::check_input)
        gpu_device_info();
#endif


#ifdef AVX
    vector_type_info();
#endif

    /* basic static node variables */
#if defined(CUDA) && !defined(PIZDAINT)
    // localhost_info(&g_local_nodeid, &g_num_local_nodes);
#endif

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
    hila::out0 << "GNU c-library performance: not returning allocated memory\n";
#endif
}


///////////////////////////////////////////////////////////////
/// Force quit for multinode processes -- kill all nodes
/// No synchronisation done
///////////////////////////////////////////////////////////////
void hila::terminate(int status) {
    hila::timestamp("Terminate");
    print_dashed_line();
    hila::about_to_finish = true; // avoid destructors
    if (is_comm_initialized()) {
        finish_communications();
    }
    exit(1);
}

////////////////////////////////////////////////////////////////
/// Print message and force quit
////////////////////////////////////////////////////////////////

void hila::error(const char *msg) {
    hila::out0 << "Error: " << msg << '\n';
    hila::terminate(0);
}

void hila::error(const std::string &msg) {
    hila::error(msg.c_str());
}

////////////////////////////////////////////////////////////////
/// Normal, controlled exit - all nodes must call this.
/// Prints timing information and information about
/// communications
////////////////////////////////////////////////////////////////

void hila::finishrun() {
    report_timers();

    for (const lattice_struct *latp : lattices) {


        int64_t gathers = latp->n_gather_done;
        int64_t avoided = latp->n_gather_avoided;

        if (gathers + avoided > 0) {
            hila::out0 << " COMMS from node 0: " << gathers << " done, " << avoided << "("
                       << 100.0 * avoided / (avoided + gathers) << "%) optimized away\n";
        } else {
            hila::out0 << " No communications done from node 0\n";
        }

    }


#if defined(CUDA) || defined(HIP)
    gpuMemPoolReport();
#endif

    if (hila::partitions.number() > 1) {
        hila::timestamp("Waiting to sync partitions");
    }

    // hip seems to want this?
    FFT_delete_plans();

    hila::synchronize();
    hila::timestamp("Finishing");

    hila::about_to_finish = true;

    finish_communications();

    print_dashed_line();
    exit(0);
}

/******************************************************
 * Open parameter file - moved here in order to
 * enable partition division if requested
 */
#if 0

FILE *open_parameter_file()
{
  static char parameter[] = "parameter";
  FILE *fil = NULL;

  if (mynode == 0) {
#ifdef SUBLATTICES
    if (n_partitions > 1) {
      char parameter_name[50];
      /* First, try opening parameter99 etc. */
      sprintf(parameter_name,"%s%d",parameter,this_partition);
      fil = fopen(parameter_name,"r");

      if (fil != NULL) {
        fprintf(outf," READING PARAMETERS FROM %s\n",parameter_name);
      }
    }
#endif
    if (fil == NULL) {
      fil = fopen(parameter,"r");
      if (fil == NULL) {
        halt(" ** No parameter file?");
      }
    }
  } // mynode == 0
  return( fil );
}

#endif

/******************************************************
 * Sublattice division
 * Handle command line arguments
 *   partitions=nn
 *   sync=yes / sync=no
 *   out=name
 * here
 */

void setup_partitions(cmdlineargs &commandline) {

    // get partitions cmdlinearg first
    long lnum = commandline.get_int("-partitions");
    if (lnum <= 0) {
        hila::out0 << "partitions=<number> command line argument value must be positive "
                      "integer (or argument omitted)\n";
        hila::finishrun();
    }
    if (lnum == LONG_MAX) {
        hila::partitions._number = 1;
    } else {
        hila::partitions._number = lnum;
    }

    if (hila::partitions.number() == 1)
        return;

    hila::out0 << " Dividing nodes into " << hila::partitions.number() << " partitions\n";

    if (hila::number_of_nodes() % hila::partitions.number()) {
        hila::out0 << "** " << hila::number_of_nodes() << " nodes not evenly divisible into "
                   << hila::partitions.number() << " partitions\n";
        hila::finishrun();
    }

#if defined(BLUEGENE_LAYOUT)
    hila::partitions._mylattice = bg_layout_partitions(hila::partitions.number());
#else // generic
    hila::partitions._mylattice =
        (hila::myrank() * hila::partitions.number()) / hila::number_of_nodes();
    /* and divide system into partitions */
    if (!hila::check_input)
        split_into_partitions(hila::partitions.mylattice());
#endif

    const char *p = commandline.get_cstring("-o");
    std::string fname;
    if (p != nullptr)
        fname = p + std::to_string(hila::partitions.mylattice());
    else
        fname = DEFAULT_OUTPUT_NAME + std::to_string(hila::partitions.mylattice());

    // now need to open output file

    hila::out.flush(); // this should be cout at this stage

    // all nodes open the file -- perhaps not?  Leave only node 0
    int do_exit = 0;
    if (hila::myrank() == 0 && !hila::check_input) {
        hila::output_file.open(fname, std::ios::out | std::ios::app);
        if (hila::output_file.fail()) {
            std::cout << "Cannot open output file " << fname << '\n';
            do_exit = 1;
        }
    }

    hila::broadcast(do_exit);

    if (do_exit)
        hila::finishrun();

    hila::out.flush();
    if (!hila::check_input) {
        hila::out.rdbuf(hila::output_file.rdbuf());
        // output now points to output_redirect
        if (hila::myrank() == 0) {
            hila::out0.rdbuf(hila::out.rdbuf());
        }
    }
    hila::out0 << " ---- SPLIT " << hila::number_of_nodes() << " nodes into "
               << hila::partitions.number() << " partitions, this " << hila::partitions.mylattice()
               << " ----\n";


    /* Default sync is no */
    if (commandline.get_onoff("-sync") == 1) {
        hila::partitions._sync = true;
        hila::out0 << "Synchronising partition trajectories\n";
    } else {
        hila::partitions._sync = false;
        hila::out0 << "Not synchronising the partition trajectories\n"
                   << "Use '-sync on' command line argument to override\n";
    }
}

/////////////////////////////////////////////////////////////////////////////

#ifdef AVX
void vector_type_info() {

    hila::out0 << "Using VCL vector class with instruction set level INSTRSET=" << INSTRSET
               << " <=> ";

    switch (INSTRSET) {
    case 2:
        hila::out0 << "SSE2";
        break;
    case 3:
        hila::out0 << "SSE3";
        break;
    case 4:
        hila::out0 << "SSSE3";
        break;
    case 5:
        hila::out0 << "SSE4.1";
        break;
    case 6:
        hila::out0 << "SSE4.2";
        break;
    case 7:
        hila::out0 << "AVX";
        break;
    case 8:
        hila::out0 << "AVX2";
        break;
    case 9:
        hila::out0 << "AVX512F";
        break;
    case 10:
        hila::out0 << "AVX512BW/DQ/VL";
        break;
    default:
        hila::out0 << "Unknown";
        break;
    }
    hila::out0 << '\n';
    if (INSTRSET < 8)
        hila::out0 << " (You probably should use options '-mavx2 -fmad' in compilation)\n";
}


#endif

void print_dashed_line(const std::string &text) {
    static constexpr int linelength = 60;

    if (hila::myrank() == 0) {

        if (text.size() == 0) {
            for (int i = 0; i < linelength; i++)
                hila::out << '-';

        } else {

            hila::out << "----- " << text << ' ';
            for (int i = 7 + text.size(); i < linelength; i++)
                hila::out << '-';
        }
        hila::out << '\n';
    }
}
