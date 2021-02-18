
#ifdef USE_MPI
// protect with USE_MPI so can be included even in no-MPI codes

#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/com_mpi.h"
#include "plumbing/timing.h"

// declare MPI timers here too - these were externs

hila::timer start_send_timer("MPI start send"), wait_send_timer("MPI wait send"),
    post_receive_timer("MPI post receive"), wait_receive_timer("MPI wait receive"),
    synchronize_timer("MPI synchronize"), reduction_timer("MPI reduction"),
    reduction_wait_timer("MPI reduction wait"), broadcast_timer("MPI broadcast"),
    send_timer("MPI send field"), cancel_send_timer("MPI cancel send"),
    cancel_receive_timer("MPI cancel receive"),
    sublattice_sync_timer("sublattice sync");

/* Keep track of whether MPI has been initialized */
static bool mpi_initialized = false;

/* Machine initialization */
#include <sys/types.h>
void initialize_communications(int &argc, char ***argv) {
    /* Init MPI */
    if (!mpi_initialized && !hila::check_input) {
        MPI_Init(&argc, argv);
        mpi_initialized = true;

        // global var lattice exists, assign the mpi comms there
        lattice->mpi_comm_lat = MPI_COMM_WORLD;

        MPI_Comm_rank(lattice->mpi_comm_lat, &lattice->mynode.rank);
        MPI_Comm_size(lattice->mpi_comm_lat, &lattice->nodes.number);
    }
    if (hila::check_input) {
        lattice->mynode.rank = 0;
        lattice->nodes.number = hila::check_with_nodes;
    }
}

// check if MPI is on
bool is_comm_initialized(void) { return mpi_initialized; }

/* version of exit for multinode processes -- kill all nodes */
void abort_communications(int status) {
    if (mpi_initialized) {
        mpi_initialized = false;
        MPI_Abort(lattices[0]->mpi_comm_lat, 0);
    }
}

/* clean exit from all nodes */
void finish_communications() {
    // turn off mpi -- this is needed to avoid mpi calls in destructors
    mpi_initialized = false;
    hila::about_to_finish = true;

    MPI_Finalize();
}

// broadcast specialization
void broadcast(std::string &var) {

    if (hila::check_input)
        return;

    int size = var.size();
    broadcast(size);

    if (hila::myrank() != 0) {
        var.resize(size, ' ');
    }
    // copy directy to data() buffer
    broadcast_timer.start();
    MPI_Bcast((void *)var.data(), size, MPI_BYTE, 0, lattice->mpi_comm_lat);
    broadcast_timer.stop();
}

void broadcast(std::vector<std::string> &list) {

    if (hila::check_input)
        return;

    int size = list.size();
    broadcast(size);
    list.resize(size);

    for (auto &s : list) {
        broadcast(s);
    }
}

/* BASIC COMMUNICATIONS FUNCTIONS */

/* Tell what kind of machine we are on */
static char name[] = "MPI (portable)";
char *machine_type() { return (name); }

/// Return my node number - take care to return
/// the previous node number if mpi is being
/// torn down (used in destructors)

int hila::myrank() {
    static int node = 0;

    if (!mpi_initialized || hila::check_input)
        return node;

    MPI_Comm_rank(lattice->mpi_comm_lat, &node);
    return node;
}

/// Return number of nodes or "pseudo-nodes"
int numnodes() {
    if (hila::check_input)
        return hila::check_with_nodes;

    int nodes;
    MPI_Comm_size(lattice->mpi_comm_lat, &nodes);
    return (nodes);
}

void synchronize() {
    synchronize_timer.start();
    synchronize_threads();
    MPI_Barrier(lattice->mpi_comm_lat);
    synchronize_timer.stop();
}

// Split the communicator to subvolumes, using MPI_Comm_split
// New MPI_Comm is the global mpi_comm_lat
// NOTE: no attempt made here to reorder the nodes

void split_into_sublattices(int this_lattice) {

    if (hila::check_input)
        return;

    if (MPI_Comm_split(MPI_COMM_WORLD, this_lattice, 0, &(lattice->mpi_comm_lat)) !=
        MPI_SUCCESS) {
        output0 << "MPI_Comm_split() call failed!\n";
        hila::finishrun();
    }
    // reset also the rank and numbers -fields
    MPI_Comm_rank(lattice->mpi_comm_lat, &lattice->mynode.rank);
    MPI_Comm_size(lattice->mpi_comm_lat, &lattice->nodes.number);
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
  mynode = hila::myrank();
}

#endif

#endif // USE_MPI
