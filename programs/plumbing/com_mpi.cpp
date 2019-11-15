
#include "../plumbing/globals.h"
#include "../plumbing/lattice.h"
#include "../plumbing/field.h"
#include "../plumbing/comm_mpi.h"


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

  }
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
    if (lattice->node_number() == 0) {
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

  exit(MPI_Finalize());
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
  MPI_Comm_rank( lattices[0]->mpi_comm_lat, &node );
  return(node);
}

/* Return number of nodes */
int numnodes()
{
  int nodes;
  MPI_Comm_size( lattices[0]->mpi_comm_lat, &nodes );
  return(nodes);
}

