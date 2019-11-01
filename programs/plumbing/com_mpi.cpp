#include "../plumbing/comm_mpi.h"

MPI_Comm mpi_comm_lat;





/* Machine initialization */
#include <sys/types.h>
void initialize_machine(int & argc, char ***argvp)
{
  /* Init MPI */
  static bool mpi_initialized = false;
  if( !mpi_initialized ){
    MPI_Init(&argc, argvp);
    mpi_initialized = true;

    /* default comm is the world */
    mpi_comm_lat = MPI_COMM_WORLD;

#ifdef SUBLATTICES
    setup_sublattices( argcp, *argvp );
#endif


  }
}




// Guarantee 64 bits for these - 32 can overflow!
unsigned long long n_gather_done = 0, n_gather_avoided = 0;




/* version of exit for multinode processes -- kill all nodes */
void terminate(int status)
{
  output0 << "Node " << mynode() << ", status = " << status << "\n";
#ifdef TIMERS
  time_stamp("Terminate");
#endif
  MPI_Abort( mpi_comm_lat, 0);
  exit(status);
}


/* clean exit from all nodes */
void finishrun()
{

#ifdef TIMERS
  report_comm_timers();
#endif

  if (mynode() == 0) {
    output0 << " COMMS from node 0: " << n_gather_done << "done,"
            <<  n_gather_avoided << "(" 
            << 100.0*n_gather_avoided/(n_gather_avoided+n_gather_done)
            << "%) optimized away\n";
  }

#ifdef TIMERS
#if defined(SUBLATTICES)
  if (n_sublattices > 1) {
    time_stamp("Waiting to sync sublattices...");
  }
#else
  time_stamp("Finishing");
#endif
#endif

  MPI_Finalize();
  exit(0);
}



/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="MPI (portable)";
char * machine_type(){
  return(name);
}

/* Return my node number */
int mynode()
{
  int node;
  MPI_Comm_rank( mpi_comm_lat, &node );
  return(node);
}

/* Return number of nodes */
int numnodes()
{
  int nodes;
  MPI_Comm_size( mpi_comm_lat, &nodes );
  return(nodes);
}

