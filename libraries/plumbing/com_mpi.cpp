
#include "plumbing/globals.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/comm_mpi.h"

#ifdef USE_MPI


/* Keep track of whether MPI has been initialized */
static bool mpi_initialized = false;

/* Machine initialization */
#include <sys/types.h>
void initialize_machine(int &argc, char ***argv)
{
  /* Init MPI */
  if( !mpi_initialized ){
    MPI_Init(&argc, argv);
    mpi_initialized = true;

#ifdef SUBLATTICES
    setup_sublattices( argcp, *argv );
#endif

#ifdef CUDA
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    initialize_cuda( rank );
#endif

  }
}

// check if MPI is on
bool is_comm_initialized(void) {
  return mpi_initialized;
}

/* version of exit for multinode processes -- kill all nodes */
void terminate(int status)
{
  if( !mpi_initialized ){
    int node;
    MPI_Comm_rank( lattices[0]->mpi_comm_lat, &node );
    output0 << "Node " << node << ", status = " << status << "\n";
#ifdef TIMERS
    time_stamp("Terminate");
#endif

    mpi_initialized = false;
    MPI_Abort( lattices[0]->mpi_comm_lat, 0);
  }
  exit(status);
}


/* clean exit from all nodes */
void finishrun()
{
  for( lattice_struct * lattice : lattices ){

#ifdef TIMERS
    report_comm_timers();
#endif

    unsigned long long gathers = lattice->n_gather_done;
    unsigned long long avoided = lattice->n_gather_avoided;
    if (lattice->node_rank() == 0) {
      output0 << " COMMS from node 0: " << gathers << " done, "
              << avoided << "(" 
              << 100.0*avoided/(avoided+gathers)
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
  }

  synchronize();
  // turn off mpi -- this is needed to avoid mpi calls in destructors
  mpi_initialized = false;

  MPI_Finalize();
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
  MPI_Comm_rank( lattice->mpi_comm_lat, &node );
  return(node);
}

/* Return number of nodes */
int numnodes()
{
  int nodes;
  MPI_Comm_size( lattice->mpi_comm_lat, &nodes );
  return(nodes);
}


void synchronize(){
  synchronize_threads();
  MPI_Barrier(lattice->mpi_comm_lat);
}



#endif // USE_MPI
