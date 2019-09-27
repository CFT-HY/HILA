
#include "field.h"
#include <cstring>

// This file is copied from our old c program: TODO:convert?

/***********************************
 * Cmdline arg operations
 * This routine returns char ptr to the argument, marked by
 * string flag.  e.g.
 *  p = process_cmdline("par=", argc, argv);
 * returns pointer to the string following string "par="
 * NO SPACES ARE ALLOWED!
 */

char * process_cmdline( char *flag, int & argc, char *argv[] )
{
  int i;
  char *p;

  for (i=1,p=nullptr; i<argc && p==nullptr; i++) {
    p = strstr(argv[i],flag);
    if (p != nullptr) {
      p += strlen(flag);
      argc-- ;
      for ( ; i<argc; i++) argv[i] = argv[i+1];
      return(p);
    }
  }
  return(nullptr);
}

/* SETUP ROUTINES */

#if (defined(__GNUC__) && !defined(DARWIN)) // || defined(__bg__)
#include "malloc.h"
#endif


void initial_setup(int & argc, char ***argvp)
{

#if (defined(__GNUC__) && !defined(DARWIN) && !defined(_MAC_OSX_)) // || defined(__bg__)
  /* First, adjust malloc so that glibc free() does not 
   * release space to the system, increasing the performance 
   * of the glib malloc substantially.  The memory use is cyclic,
   * so we can just sit on the max memory.
   */
  mallopt( M_MMAP_MAX, 0 );  /* don't use mmap */
  /* HACK: don't release memory by calling sbrk */
  mallopt( M_TRIM_THRESHOLD, -1 );
#endif

  // Default output file - we're happy with this unless sublattices
  // This channels outf to std::cout
  hila::output( std::cout.rdbuf() );

  /* Machine initialization first -- open up MPI nodes, copying
   * argcp, argv to all
   *
   * Also sets up sublattices (setup_sublattices) inside
   */

  initialize_machine( argc, argvp );

#ifndef SUBLATTICES
  /* Check output file name -- NULL: stdout
   * NOTE: if using sublattices, setup_sublattices took care of this
   * already
   */

  if (mynode() == 0) {
    if (char *name = process_cmdline("out=",argc,*argvp)) {
      // Open file for append
      hila::output_redirect.open(name, ios::out | ios::app);
      hila::output( hila::output_file.rdbuf() );  // output now points to output_redirect
      if (hila::output.fail()) {
        std::cout << "HiLa: cannot open output file " << name << '\n';        
        exit(1);
      }
    }
  }
#endif  // not SUBLATTICES

  /* set the timing up */
  timing.init();

  /* basic static node variables */
  number_of_nodes = numnodes();
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

  n = mynode();
#ifdef SUBLATTICES
  n += this_sublattice * numnodes();
#endif

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

  output0 << " Using node random numbers, seed for node 0: " << seed << '\n'';

  seed_mersenne(seed);
  // warm it up
  for (i=0; i<90000; i++) dran();
  taus_initialize();

#else // Now SITERAND is defined
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


void deinit_rndgen (void)
{
    taus_deinit();
}


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

#ifdef SUBLATTICES

void setup_sublattices(int & argc, char ** argv)
{
  char *p;

  sublattices.number = 0;
  /* Open new subvolumes here */

  // fprintf(outf,"Node %d, nargs %d, arg0 %s\n", this_node, *argcp, argv[0]);

  char *p = process_cmdline("sublattices=",argc,argv);
  if (p) sublattices.number = (int)strtol(p,NULL,10);

  if (sublattices.number > 1) {
    output0 << " Dividing nodes into " << sublattice.number << " sublattices\n";
  } else {
    output0 << " NO SUBLATTICES.  Use sublattices=nn command line argument to change\n";
  }


  if (numnodes() % sublattices.number) {
    output0 << "** " << numnodes() << " nodes not evenly divisible into " 
            << sublattices.number <<  " sublattices\n";
    finishrun();
  }

  /* Default sync is no */
  sublattices.sync = false;
  if (process_cmdline("sync=yes",argc,argv) != nullptr) {
    sublattices.sync = true;
    output0 << " Synchronising sublattice trajectories\n";
  } else {
    process_cmdline("sync=no",argc,argv);
    output0 << "Not synchronising the sublattice trajectories\n"
            << "Use sync=yes command line argument to override\n";
  }

  #if defined(BLUEGENE_LAYOUT)
  sublattices.mylattice = bg_layout_sublattices( sublattices.number );
  #else // generic
  sublattices.mylattice = (mynode()*sublattices.number) / n;
  /* and divide system into sublattices */
  split_into_sublattices( sublattices.mylattice );
  #endif

  /* Set output file here (crazy, but needs to be done rapidly! 
   */

  p = process_cmdline("out=",argc,argv);
    std::string fname;
    if (p! = nullptr) fname = p + sublattices.mylattice; 
    else fname = "output" + sublattices.mylattice;

    hila::output.flush();   // this should be cout at this stage
    
    // all nodes open the file -- perhaps not?  Leave only node 0
    if (mynode() == 0) {
      hila::output_redirect.open(fname, ios::out | ios::app);
      hila::output( hila::output_redirect.rdbuf() );  // output now points to output_redirect
      if (hila::output.fail()) {
        std::cout << "HiLa: cannot open output file " << fname << '\n';        
        exit(1);
      }
      
      hila::output << " ---- SPLIT " << numnodes() << " nodes into " 
                   << sublattices.number << " sublattices, this " << sublattices.mylattice << " ----\n";
      if (!sublattices.sync) {
        hila::output << "      Not synchronizing sublattice trajectories\n";
        hila::output << "      Use sync=yes command line argument to override\n";
      } else {
        hila::output << "      Synchronising sublattice trajectories\n";
      }
    }

  
}

#endif // SUBLATTICES





