
#include <cstring>
#include "hila.h"

// define these global var here - somehow NULL needed for ostream
std::ostream hila::out(NULL);
std::ostream hila::out0(NULL);
std::ofstream hila::output_file;
bool hila::about_to_finish = false;
bool hila::is_initialized = false;
bool hila::check_input = false;
int hila::check_with_nodes;
const char *hila::input_file;
logger_class hila::log;


void vector_type_info();

#include <limits.h>
#include <errno.h>

int get_onoff(std::string flag)
{
    // Check if flag has been set
    if (hila::cmdline.flag_set(flag.c_str()))
    {
        std::string opt = hila::cmdline.get_string(flag.c_str());
        if (opt.compare("on") == 0)
            return 1;
        else if (opt.compare("off") == 0)
            return -1;
        else
        {
            hila::out0 << "Command line argument " << flag << " requires value on/off\n";
            hila::terminate(0);
            return 0;
        }
    }
    else
        return 0;
}

/////////////////////////////////////////////////////////////////////////////////
/// Initial setup routines
/////////////////////////////////////////////////////////////////////////////////

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
#include <malloc.h>
#endif

// #define DEBUG_NAN

#ifdef DEBUG_NAN
#include <fenv.h>
#endif

void setup_partitions();

/**
 * @brief Read in command line arguments. Initialise default stream and MPI communication
 *
 * @param argc Number of command line arguments
 * @param argv List of command line arguments
 */

void hila::initialize(int argc, char **argv) {

#if (defined(__GNUC__) && !defined(DARWIN) && !defined(_MAC_OSX_)) // || defined(__bg__)
    /* First, adjust malloc so that glibc free() does not
     * release space to the system, increasing the performance
     * of the glib malloc substantially.  The memory use is cyclic,
     * so we can just sit on the max memory.
     */
    mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
    /* HACK: don't release memory by calling sbrk */
    mallopt(M_TRIM_THRESHOLD, -1);

#ifdef DEBUG_NAN
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
#endif

    // initialize MPI so that hila::myrank() etc. works
    initialize_communications(argc, &argv);

    // set the initialized flag
    hila::is_initialized = true;

    // Default output file - we're happy with this unless partitions
    // or otherwise indicated
    // This channels outf to std::cout
    hila::out.rdbuf(std::cout.rdbuf());

    // set the timing so that gettime() returns time from this point
    hila::inittime();

    // open hila::out0 only for node 0
    if (hila::myrank() == 0)
        hila::out0.rdbuf(std::cout.rdbuf());

    // Set the inbuilt command-line flags and their corresponding help texts
    hila::cmdline.add_flag("-t","cpu time limit");
    hila::cmdline.add_flag("-o","output filename (default: stdout)");
    hila::cmdline.add_flag("-i","input filename (overrides the 1st hila::input() name)\nuse '-i -' for standard input");
    hila::cmdline.add_flag("-device","in GPU runs using only 1 GPU, choose this GPU number (default 0)");
    hila::cmdline.add_flag("-check","check input & layout with <nodes>-nodes & exit\nonly with 1 real MPI node (don't use mpirun)");
    hila::cmdline.add_flag("-n","number of nodes used in layout check, only relevant with -check");
    hila::cmdline.add_flag("-partitions","number of partitioned lattice streams");
    hila::cmdline.add_flag("-sync","synchronize partition runs (on/off) (default = off)");

    // Init command line - after MPI has been started, so
    // that all nodes do this. First feed argc and argv to the
    // global cmdline class instance and parse for the preset flags.
    hila::cmdline.initialise_args(argc, argv);
    // The values can now be requested from hila::cmdline.

    // check the "-check" -input early
    // do it only with 1 node
    if (lattice.nodes.number == 1) {
        // Check whether '-check' was found and only then search for '-n'
        if (hila::cmdline.flag_present("-check")) {
            long nodes;
            if (hila::cmdline.flag_present("-n")) {
                nodes = hila::cmdline.get_int("-n");
            }
            else
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
        long device;
        if (hila::cmdline.flag_set("-device"))
            device = hila::cmdline.get_int("-device");
        else
            device = 0;
        hila::out0 << "Chose device " << device << "\n";

        initialize_gpu(lattice.mynode.rank, device);
    }
#endif

    setup_partitions();

    // check the output file if partitions not used
    if (hila::partitions.number() == 1) {
        int do_exit = 0;
        if (hila::myrank() == 0) {
            if (hila::cmdline.flag_present("-o")) {
                // Quits if '-o' was left without an argument
                std::string name;
                if (hila::cmdline.flag_set("-o"))
                    name = hila::cmdline.get_string("-o");
                else
                {
                    hila::out0 << "The name of the output file must be provided after flag '-o'!\n";
                    do_exit = 1;
                }
                // If found, open the file for the output
                if (!hila::check_input) {
                    hila::output_file.open(name, std::ios::out | std::ios::app);
                    if (hila::output_file.fail()) {
                        hila::out << "Cannot open output file " << name << '\n';
                        do_exit = 1;
                    } else {
                        hila::out0 << "Output is now directed to the file '"
                                   << name << "'.\n";
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
        print_dashed_line("HILA lattice framework");
        hila::out0 << "Running program " << argv[0] << "\n";
        hila::out0 << "with command line arguments '";
        for (int i = 1; i < argc; i++)
            hila::out0 << argv[i] << ' ';
        hila::out0 << "'\n";
        hila::out0 << "Code version: ";
#if defined(GIT_SHA_VALUE)
#define xstr(s) makestr(s)
#define makestr(s) #s
        hila::out0 << "git SHA " << xstr(GIT_SHA_VALUE) << '\n';
#else
        hila::out0 << "no git information available\n";
#endif
        hila::out0 << "Compiled " << __DATE__ << " at " << __TIME__ << '\n';

        hila::out0 << "with options: EVEN_SITES_FIRST";
#ifndef EVEN_SITES_FIRST
        hila::out0 << "=0";
#endif
#ifdef SPECIAL_BOUNDARY_CONDITIONS
        hila::out0 << " SPECIAL_BOUNDARY_CONDITIONS";
#endif
        hila::out0 << '\n';

        hila::timestamp("Starting");
    }

    // Check if flag set and parse
    if (hila::cmdline.flag_present("-t")) {
        // Following quits if '-t' is given without an integer argument
        long cputime = hila::cmdline.get_int("-t");
        if (cputime > 0) {
            hila::out0 << "CPU time limit " << cputime << " seconds\n";
            hila::setup_timelimit(cputime);
        }
    }
    else {
        hila::out0 << "No runtime limit given\n";
    }


    hila::input_file = nullptr;
    if (hila::cmdline.flag_present("-i"))
    {
        // Quits if '-i' given without a string argument
        // Copy to a static variable to preserve the memory address
        static const std::string input_string = hila::cmdline.get_string("-i");
        hila::input_file = input_string.c_str();
        hila::out0 << "Input file from command line: " << hila::input_file << "\n";
    }


#if defined(OPENMP)
    hila::out0 << "Using option OPENMP - with " << omp_get_max_threads() << " threads\n";
#endif


#if defined(CUDA) || defined(HIP)
    hila::out0 << "Using thread blocks of size " << N_threads << " threads\n";

#if defined(GPU_AWARE_MPI)
    hila::out0 << "Using GPU_AWARE_MPI\n";
#else
    hila::out0 << "Not using GPU_AWARE_MPI\n";
#endif

#if !defined(GPU_VECTOR_REDUCTION_THREAD_BLOCKS) || GPU_VECTOR_REDUCTION_THREAD_BLOCKS <= 0
    hila::out0 << "ReductionVector with atomic operations (GPU_VECTOR_REDUCTION_THREAD_BLOCKS=0)\n";
#else
    hila::out0 << "ReductionVector with " << GPU_VECTOR_REDUCTION_THREAD_BLOCKS
               << " thread blocks\n";
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

/**
 * @brief Normal, controlled exit - all nodes must call this.
 * Prints timing information and information about communications
 */
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

void setup_partitions() {

    // get partitions cmdlinearg first
    if (hila::cmdline.flag_present("-partitions")) {
        // Following quits if '-partitions' is given without an integer argument
        long lnum = hila::cmdline.get_int("-partitions");
        if (lnum <= 0) {
            hila::out0 << "partitions=<number> command line argument value must be positive "
                          "integer (or argument omitted)\n";
            hila::finishrun();
        }
        else
            hila::partitions._number = lnum;
    }
    else
        hila::partitions._number = 1;

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
    std::string fname;
    if (hila::cmdline.flag_present("-o"))
    {
        std::string opt = hila::cmdline.get_string("-o");
        fname = opt + std::to_string(hila::partitions.mylattice());
    }
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
    if (hila::cmdline.flag_present("-sync"))
    {
        std::string onoffopt = hila::cmdline.get_string("-sync");
        if (get_onoff(onoffopt) == 1)
        {
            hila::partitions._sync = true;
            hila::out0 << "Synchronising partition trajectories\n";
        }
    }
    else
    {
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
