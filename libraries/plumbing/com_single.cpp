
#ifndef USE_MPI
// protect with USE_MPI so can be included even in no-MPI codes

#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"

static bool initialized = false;

/* Machine initialization */
#include <sys/types.h>
void initialize_machine(int &argc, char ***argv) {
    /* Init MPI */
    if (!initialized) {
        lattice->mynode.rank = 0;
        lattice->nodes.number = 1;

#if defined(CUDA) || defined(HIP)
        initialize_cuda(0);
#endif
        initialized = true;
    }
    hila::my_rank_n = 0;
}

// check if is intialized
bool is_comm_initialized(void) { return initialized; }

/* version of exit for multinode processes -- kill all nodes */
void abort_communications(int status) {}

/* clean exit from all nodes */
void finish_communications() {}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
char *machine_type() {
    static char name[] = "SINGLE (no MPI)";
    return name;
}

/* Return my node number */
int hila::myrank() { return 0; }

/* Return number of nodes */
int hila::number_of_nodes() { return 1; }

void hila::synchronize() { hila::synchronize_threads(); }

// Split the communicator to subvolumes, using MPI_Comm_split
// New MPI_Comm is the global mpi_comm_lat
// NOTE: no attempt made here to reorder the nodes

void split_into_sublattices(int this_lattice) {}

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
  mynode = hila::myrank();
}

#endif

#endif // USE_MPI
