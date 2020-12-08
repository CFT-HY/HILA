
#ifdef USE_MPI  
// protect with USE_MPI so can be included even in no-MPI codes


#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/com_mpi.h"
#include "plumbing/timing.h"

// declare MPI timers here too - these were externs

timer start_send_timer, wait_send_timer, post_receive_timer, wait_receive_timer,
      synchronize_timer, reduction_timer, broadcast_timer, send_timer, 
      cancel_send_timer, cancel_receive_timer, sublattice_sync_timer;

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

    // global var lattice exists, assign the mpi comms there
    lattice->mpi_comm_lat = MPI_COMM_WORLD;

    MPI_Comm_rank( lattice->mpi_comm_lat, &lattice->this_node.rank );
    MPI_Comm_size( lattice->mpi_comm_lat, &lattice->nodes.number );

#ifdef CUDA
    initialize_cuda( lattice->this_node.rank );
#endif

  }
  hila::my_rank_n = lattice->this_node.rank;

  start_send_timer.init("MPI start send");
  wait_send_timer.init("MPI wait send");
  post_receive_timer.init("MPI post receive");
  wait_receive_timer.init("MPI wait receive");
  synchronize_timer.init("MPI synchronize");
  reduction_timer.init("MPI reduction");
  broadcast_timer.init("MPI broadcast");
  send_timer.init("MPI send field");
  cancel_send_timer.init("MPI cancel send");
  cancel_receive_timer.init("MPI cancel receive");
  sublattice_sync_timer.init("Sublattice sync");

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

    timestamp("Terminate");

    mpi_initialized = false;
    MPI_Abort( lattices[0]->mpi_comm_lat, 0);
  }
  exit(status);
}

void error(const char * msg) {
  output0 << "Error: " << msg << '\n';
  terminate(0);
}

void error(const std::string &msg) {
  error(msg.c_str());
}


/* clean exit from all nodes */
void finishrun()
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
  } else {
    timestamp("Finishing");
  }

  synchronize();
  // turn off mpi -- this is needed to avoid mpi calls in destructors
  mpi_initialized = false;

  MPI_Finalize();
}

// broadcast specialization
void broadcast(std::string & var) {
  broadcast_timer.start();
  int size = var.size();
  broadcast(size);

  if (hila::myrank()!=0) {
    var.resize(size,' ');
  }
  // copy directy to data() buffer
  MPI_Bcast((void *)var.data(),size, MPI_BYTE, 0, lattice->mpi_comm_lat);
}

void broadcast(std::vector<std::string> & list) {
  int size = list.size();
  broadcast(size);
  list.resize(size);

  for (auto & s : list) {
    broadcast(s);
  }
}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[]="MPI (portable)";
char * machine_type()
{
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


void synchronize()
{
  synchronize_timer.start();
  synchronize_threads();
  MPI_Barrier(lattice->mpi_comm_lat);
  synchronize_timer.stop();
}


// Split the communicator to subvolumes, using MPI_Comm_split
// New MPI_Comm is the global mpi_comm_lat
// NOTE: no attempt made here to reorder the nodes

void split_into_sublattices( int this_lattice ) 
{
  if ( MPI_Comm_split( MPI_COMM_WORLD, this_lattice, 0, &(lattice->mpi_comm_lat) )
       != MPI_SUCCESS ) {
    output0 << "MPI_Comm_split() call failed!\n";
    finishrun();
  }
  // reset also the rank and numbers -fields
  MPI_Comm_rank( lattice->mpi_comm_lat, &lattice->this_node.rank );
  MPI_Comm_size( lattice->mpi_comm_lat, &lattice->nodes.number );
}

#if 0

// Switch comm frame global-sublat
// for use in io_status_file

void reset_comm(bool global)
{
  static MPI_Comm mpi_comm_saved;
  static int set = 0;

  g_sync_sublattices();
  if (global) {
    mpi_comm_saved = lattice->mpi_comm_lat;
    lattice->mpi_comm_lat = MPI_COMM_WORLD;
    set = 1;
  } else {
    if (!set) halt("Comm set error!");
    lattice->mpi_comm_lat = mpi_comm_saved;
    set = 0;
  }
  this_node = mynode();
}

#endif



#endif // USE_MPI
