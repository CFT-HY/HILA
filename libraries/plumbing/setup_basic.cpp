
#include "field.h"
#include <cstring>


// let us house the sublattices-struct here

struct sublattices_struct {
  unsigned number,mylattice;
  bool sync;
} sublattices;

// This file is largely copied from our old c program: TODO:convert?

///////////////////////////////////////////////////////////////////////////////
// very simple cmdline arg interpreter
// p = cmdline.get_arg("par="); 
// returns pointer to string following "par="
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
  cmdlineargs(int argc0, const char **argv0) {
    argc = argc0;
    argv = (const char **)malloc(argc*sizeof(const char *));
    for (int i=0; i<argc; i++) argv[i] = argv0[i];
  }

  ~cmdlineargs() {
    free(argv);
  }

  const char *get_arg(const char *flag) {
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
    const char *p = get_arg(flag);
    char *end;
    if (!p) return LONG_MAX;  // use LONG_MAX to denote no value

    long val = strtol(p, &end, 10);
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0) || end == p || *end != 0) {
      if (mynode() == 0) {
        hila::output << "Expect a number (integer) after command line parameter '" << flag << "'\n";
      }
      exit(0);      
    }
  }

  int items() {
    return argc-1;   // don't count argv[0]
  }

}


/* SETUP ROUTINES */

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
#include "malloc.h"
#endif


void setup_sublattices(cmdlineargs & cl);


void initial_setup(int & argc, char **argv)
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
  hila::output( std::cout.rdbuf() );

  // initialize MPI (if needed) so that mynode() etc. works
  initialize_machine( argc, &argv );

  cmdlineargs commandline(argc, argv);

  setup_sublattices(commandline);

  // check the output file if sublattices not used
  if (sublattices.number == 1 && mynode() == 0) {
    if (char *name = commandline.get_arg("output=")) {
      // Open file for append
      hila::output_redirect.open(name, ios::out | ios::app);
      hila::output( hila::output_file.rdbuf() );  // output now points to output_redirect
      if (hila::output.fail()) {
        std::cout << "hila: cannot open output file " << name << '\n';        
        exit(1);
      }
    }
  }

  /* set the timing up */
  timing.init();

  /* basic static node variables */
#if defined(CUDA) && !defined(PIZDAINT)
  localhost_info(&g_local_nodeid, &g_num_local_nodes);
#endif

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
  output0 << "GNU c-library performance:\n using sbrk instead of mmap; not returning memory\n";
#endif

  output0 << " ### Running target " << (*argvp)[0] << " ###\n";

}


/**************************************************
 * random number generators
 */

void initialize_prn(long seed)
{
  int n;
  int i;

#ifndef SITERAND

  int n = mynode();
  if (sublattices.number > 1)
  n += sublattices.mylattice * numnodes();

  if (seed == 0) {
    if (mynode() == 0) {
      seed = time(NULL);
      seed = seed^(seed<<26)^(seed<<9);
      output0 << "Random seed from time " << seed << '\n';
    }
    broadcast_field(&seed,sizeof(long));
  }
  seed += 1121*n;
  seed = seed ^ ((511*n)<<18);

  output0 << " Using node random numbers, seed for node 0: " << seed << '\n';

  seed_mersenne(seed);
  // warm it up
  for (i=0; i<90000; i++) dran();

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


/******************************************************
 * Open parameter file - moved here in order to
 * enable sublattice division if requested
 */

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
    finishrun();
  }
  if (lnum == LONG_MAX) {
    sublattices.number = 1;
  } else {
    sublattices.number = lnum;
  }

  if (sublattices.number == 1) return;

  output0 << " Dividing nodes into " << sublattice.number << " sublattices\n";

  if (numnodes() % sublattices.number) {
    output0 << "** " << numnodes() << " nodes not evenly divisible into " 
            << sublattices.number <<  " sublattices\n";
    finishrun();
  }

#if defined(BLUEGENE_LAYOUT)
  sublattices.mylattice = bg_layout_sublattices( sublattices.number );
#else // generic
  sublattices.mylattice = (mynode()*sublattices.number) / numnodes();
  /* and divide system into sublattices */
  split_into_sublattices( sublattices.mylattice );
#endif

  const char * p = commandline.get_arg("output=");
  std::string fname;
  if (p != nullptr) fname = p + sublattices.mylattice; 
  else fname = DEFAULT_OUTPUT_NAME + sublattices.mylattice;

  // now need to open output file

  hila::output.flush();   // this should be cout at this stage
    
  // all nodes open the file -- perhaps not?  Leave only node 0
  if (mynode() == 0) {
    hila::output_redirect.open(fname, ios::out | ios::app);
    hila::output( hila::output_redirect.rdbuf() );  // output now points to output_redirect
    if (hila::output.fail()) {
      std::cout << "Cannot open output file " << fname << '\n';        
      exit(1);
    }
      
    hila::output << " ---- SPLIT " << numnodes() << " nodes into " 
                 << sublattices.number << " sublattices, this " << sublattices.mylattice << " ----\n";
  }
  
  /* Default sync is no */
  sublattices.sync = false;
  if (commandline.get_arg("sync=yes") != nullptr) {
    sublattices.sync = true;
    output0 << " Synchronising sublattice trajectories\n";
  } else {
    process_cmdline("sync=no",argc,argv);
    output0 << "Not synchronising the sublattice trajectories\n"
            << "Use sync=yes command line argument to override\n";
  }
  
}

#endif // SUBLATTICES





