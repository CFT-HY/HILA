
#include <cstring>
#include "defs.h"
#include "field.h"
#ifdef USE_MPI
#include "com_mpi.h"
#endif

// define these global var here - somehow NULL needed for ostream
std::ostream hila::output(NULL);
std::ofstream hila::output_file;
int hila::my_rank_n;
bool hila::about_to_finish = false;

// let us house the sublattices-struct here

sublattices_struct sublattices;


void vector_type_info();

///////////////////////////////////////////////////////////////////////////////
// very simple cmdline arg interpreter
// p = cmdline.get_cstring("par="); 
// returns a string of the stuff following "par="
// n = cmdline.get_int("par=");
// returns (long) int containing the value following "par="
//
// NO SPACES are allowed between par and = and the value
///////////////////////////////////////////////////////////////////////////////

#include <limits.h>
#include <errno.h>

class cmdlineargs {
 private:
  int argc;
  const char **argv;

 public:
  cmdlineargs(int argc0, char **argv0) {
    argc = argc0;
    argv = (const char **)malloc(argc*sizeof(const char *));
    for (int i=0; i<argc; i++) argv[i] = argv0[i];
  }

  ~cmdlineargs() {
    free(argv);
  }

  const char * get_cstring(const char *flag) {
    int flaglen = strlen(flag);
    assert(flaglen > 0);

    for (int i=1; i<argc; i++) {
      const char *p = argv[i];

      // OK if p starts with flag, last char is = or next is 0
      if (strncmp(p,flag,flaglen) == 0 && 
          (*(p+flaglen-1) == '=' || *(p+flaglen) == 0)) {
        p += flaglen;
        argc-- ;
        for ( ; i<argc; i++) argv[i] = argv[i+1];
        return(p);  
      }
    }
    return(nullptr);
  }

  long get_int(const char *flag) {
    const char *p = get_cstring(flag);
    char *end;
    if (!p) return LONG_MAX;  // use LONG_MAX to denote no value

    long val = strtol(p, &end, 10);
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0) || end == p || *end != 0) {
      output0 << "Expect a number (integer) after command line parameter '" << flag << "'\n";
      hila::terminate(0);      
    }
    return val;
  }

  /// returns 1=yes, -1=no, 0=not found
  int get_yesno(const char *flag) {
    const char *p = get_cstring(flag);
    if (!p) return 0;
    if (std::strcmp(p,"yes") == 0) return 1;
    if (std::strcmp(p,"no") == 0) return -1;
    output0 << "Command line argument " << flag << " requires value yes/no\n";
    hila::terminate(0);
    return 0;   // gets rid of a warning of no return value
  }

  int items() {
    return argc-1;   // don't count argv[0]
  }

  const char * get_item(int i) {
    if (i>=0 && i<argc-1) 
      return argv[i+1];
    else 
      return nullptr;
  }

  void error_if_args_remain() {
    if (argc < 2) return;
    output0 << "Unknown command line arguments:\n";
    for (int i=1; i<argc; i++) {
      output0 << "    " << argv[i] << '\n';
    }
    output0 << "Recognized:\n";
    output0 << "  timelimit=<seconds>     cpu time limit\n";
    output0 << "  output=<name>           name of output file (default: stdout/output)\n";
    output0 << "  sublattices=<n>         number of sublattices\n";
    output0 << "  sync=yes/no             synchronize sublattice runs (default=no)\n";

    hila::terminate(0);
  }

};


/////////////////////////////////////////////////////////////////////////////////
/// Initial setup routines
/////////////////////////////////////////////////////////////////////////////////
  

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
#include <malloc.h>
#endif


void setup_sublattices(cmdlineargs & cl);


void hila::initialize(int argc, char **argv)
{

#if (defined(__GNUC__) && !defined(DARWIN) && !defined(_MAC_OSX_)) // || defined(__bg__)
  /* First, adjust malloc so that glibc free() does not 
   * release space to the system, increasing the performance 
   * of the glib malloc substantially.  The memory use is cyclic,
   * so we can just sit on the max memory.
   */
  //mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
  /* HACK: don't release memory by calling sbrk */
  mallopt( M_TRIM_THRESHOLD, -1 );
#endif
  
  // Default output file - we're happy with this unless sublattices
  // or otherwise indicated
  // This channels outf to std::cout
  hila::output.rdbuf( std::cout.rdbuf() );

  // set the timing so that gettime() returns time from this point
  inittime();

  // initialize MPI (if needed) so that hila::myrank() etc. works
  initialize_communications( argc, &argv );

#ifdef CUDA
  initialize_cuda( lattice->this_node.rank );
#endif

  /// Handle commandline args here
  cmdlineargs commandline(argc, argv);

  setup_sublattices(commandline);

  // check the output file if sublattices not used
  if (sublattices.number == 1) {
    int do_exit = 0;
    if (hila::myrank() == 0) {
      if (const char *name = commandline.get_cstring("output=")) {
        // Open file for append
        if (std::strlen(name) == 0) {
          hila::output << "Filename must be given with output=<name>\n"; 
          do_exit = 1;
        } else {
          hila::output_file.open(name, std::ios::out | std::ios::app);
          if (hila::output_file.fail()) {
            hila::output << "Cannot open output file " << name << '\n';        
            do_exit = 1;
          } else {
            hila::output.flush();
            hila::output.rdbuf( hila::output_file.rdbuf() );  // output now points to output_redirect
          }
        }
      }
    }
    broadcast(do_exit);
    if (do_exit) hila::terminate(0);
  }

  if (hila::myrank() == 0) {
    hila::output << "------------- Hila lattice program --------------\n";
    hila::output << "Running target " << argv[0] << "\n";
    hila::output << "with command line arguments '";
    for (int i=1; i<argc; i++) hila::output << argv[i] << ' ';
    hila::output << "'\n";
    hila::output << "Code version: ";
    #if defined(GIT_SHA_VALUE)
    #define xstr(s) makestr(s)
    #define makestr(s) #s
    hila::output << "git SHA " << xstr(GIT_SHA_VALUE) << '\n';
    #else
    hila::output << "no git information available\n";
    #endif
    hila::output << "  Compiled " << __DATE__ << " at " << __TIME__ << '\n';

    timestamp("Starting");
  }

  long cputime = commandline.get_int("timelimit=");
  if (cputime != LONG_MAX) {
    output0 << "CPU time limit " << cputime << " seconds\n";
    setup_timelimit(cputime);
  } else {
    output0 << "No runtime limit given\n";
  }

  // error out if there are more cmdline options
  commandline.error_if_args_remain();

#ifdef CUDA
  cuda_device_info();
#endif

#ifdef AVX
  vector_type_info();
#endif

  /* basic static node variables */
#if defined(CUDA) && !defined(PIZDAINT)
  // localhost_info(&g_local_nodeid, &g_num_local_nodes);
#endif

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
  output0 << "GNU c-library performance:\n using sbrk instead of mmap; not returning memory\n";
#endif


}


/**************************************************
 * random number generators
 */

void initialize_prn(long seed)
{

#ifndef SITERAND

  int n = hila::myrank();
  if (sublattices.number > 1)
  n += sublattices.mylattice * numnodes();

  if (seed == 0) {
    if (hila::myrank() == 0) {
      seed = time(NULL);
      seed = seed^(seed<<26)^(seed<<9);
      output0 << "Random seed from time " << seed << '\n';
    }
    broadcast(seed);
  }
  seed += 1121*n;
  seed = seed ^ ((511*n)<<18);

  output0 << " Using node random numbers, seed for node 0: " << seed << '\n';

  seed_mersenne(seed);
  // warm it up
  for (int i=0; i<90000; i++) mersenne();

  // taus_initialize();

#else 
  // Now SITERAND is defined
  // This is usually used only for occasional benchmarking, where identical output
  // independent of the node number is desired

  output0 << "*** SITERAND is in use!\n";

  random_seed_arr =
    (unsigned short (*)[3])memalloc(3*node.sites*sizeof(unsigned short));
  forallsites(i) {
    random_seed_arr[i][0] = (unsigned short)(seed +   site[i].index);
    random_seed_arr[i][1] = (unsigned short)(seed + 2*site[i].index);
    random_seed_arr[i][2] = (unsigned short)(seed + 3*site[i].index);
  }

  random_seed_ptr = random_seed_arr[0];

#endif
}


// void deinit_rndgen (void)
// {
//     taus_deinit();
// }


/* version of exit for multinode processes -- kill all nodes */
void hila::terminate(int status)
{
  timestamp("Terminate");
  print_dashed_line();
  hila::about_to_finish = true;   // avoid destructors 
  if( is_comm_initialized() ){
    abort_communications(status);
  }
  exit(status);
}

void hila::error(const char * msg) {
  output0 << "Error: " << msg << '\n';
  hila::terminate(0);
}

void hila::error(const std::string &msg) {
  hila::error(msg.c_str());
}

// Normal, controlled exit of the program
void hila::finishrun()
{
  report_timers();

  for( lattice_struct * lattice : lattices ){

    unsigned long long gathers = lattice->n_gather_done;
    unsigned long long avoided = lattice->n_gather_avoided;
    if (lattice->node_rank() == 0) {
      output0 << " COMMS from node 0: " << gathers << " done, "
              << avoided << "(" 
              << 100.0*avoided/(avoided+gathers)
              << "%) optimized away\n";
    }
  }
  if (sublattices.number > 1) {
    timestamp("Waiting to sync sublattices...");
  }  
  synchronize();
  timestamp("Finishing");

  hila::about_to_finish = true;

  finish_communications();

  print_dashed_line();
  exit(0);

}




/******************************************************
 * Open parameter file - moved here in order to
 * enable sublattice division if requested
 */
#if 0

FILE *open_parameter_file()
{
  static char parameter[] = "parameter";
  FILE *fil = NULL;

  if (this_node == 0) {
#ifdef SUBLATTICES
    if (n_sublattices > 1) {
      char parameter_name[50];
      /* First, try opening parameter99 etc. */
      sprintf(parameter_name,"%s%d",parameter,this_sublattice);
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
  } // this_node == 0
  return( fil );
}

#endif

/******************************************************
 * Sublattice division
 * Handle command line arguments
 *   sublattices=nn
 *   sync=yes / sync=no
 *   out=name
 * here
 */

void setup_sublattices(cmdlineargs & commandline)
{

 // get sublattices cmdlinearg first
  long lnum = commandline.get_int("sublattices=");
  if (lnum <= 0) {
    output0 << "sublattices=<number> command line argument value must be positive integer (or argument omitted)\n";
    hila::finishrun();
  }
  if (lnum == LONG_MAX) {
    sublattices.number = 1;
  } else {
    sublattices.number = lnum;
  }

  if (sublattices.number == 1) return;

  output0 << " Dividing nodes into " << sublattices.number << " sublattices\n";

  if (numnodes() % sublattices.number) {
    output0 << "** " << numnodes() << " nodes not evenly divisible into " 
            << sublattices.number <<  " sublattices\n";
    hila::finishrun();
  }

#if defined(BLUEGENE_LAYOUT)
  sublattices.mylattice = bg_layout_sublattices( sublattices.number );
#else // generic
  sublattices.mylattice = (hila::myrank()*sublattices.number) / numnodes();
  /* and divide system into sublattices */
  split_into_sublattices( sublattices.mylattice );
#endif

  const char * p = commandline.get_cstring("output=");
  std::string fname;
  if (p != nullptr) 
    fname = p + std::to_string( sublattices.mylattice ); 
  else 
    fname = DEFAULT_OUTPUT_NAME + std::to_string( sublattices.mylattice );

  // now need to open output file

  hila::output.flush();   // this should be cout at this stage
    
  // all nodes open the file -- perhaps not?  Leave only node 0
  int do_exit = 0;
  if (hila::myrank() == 0) {
    hila::output_file.open(fname, std::ios::out | std::ios::app);
    if (hila::output_file.fail()) {
      std::cout << "Cannot open output file " << fname << '\n';        
      do_exit = 1;
    } else {
      hila::output.flush();
      hila::output.rdbuf( hila::output_file.rdbuf() );  // output now points to output_redirect
      hila::output << " ---- SPLIT " << numnodes() << " nodes into " 
                   << sublattices.number << " sublattices, this " << sublattices.mylattice << " ----\n";
    }
  }
  broadcast(do_exit);
  if (do_exit) hila::terminate(0);

  /* Default sync is no */
  if (commandline.get_yesno("sync=") == 1) {    
    sublattices.sync = true;
    output0 << "Synchronising sublattice trajectories\n";
  } else {
    sublattices.sync = false;
    output0 << "Not synchronising the sublattice trajectories\n"
            << "Use sync=yes command line argument to override\n";
  }
  
}


#ifdef AVX
void vector_type_info() {

  output0 << "Using 'vectorclass' with instruction set level INSTRSET=" 
          << INSTRSET << " <=> ";

  switch (INSTRSET) {
    case 2: output0 << "SSE2"; break;
    case 3: output0 << "SSE3"; break;
    case 4: output0 << "SSSE3"; break;
    case 5: output0 << "SSE4.1"; break;
    case 6: output0 << "SSE4.2"; break;
    case 7: output0 << "AVX"; break;
    case 8: output0 << "AVX2"; break;
    case 9: output0 << "AVX512F"; break;
    case 10: output0 << "AVX512BW/DQ/VL"; break;
    default: output0 << "Unknown"; break;
  }
  output0 << '\n';
  if (INSTRSET < 8) output0 << " (You probably should use options '-mavx2 -fmad' in compilation)\n";
}

#endif