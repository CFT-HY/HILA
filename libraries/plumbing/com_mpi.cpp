
#include "plumbing/defs.h"
#include "plumbing/lattice.h"
#include "plumbing/field.h"
#include "plumbing/com_mpi.h"
#include "plumbing/timing.h"

#ifdef GPU_CCL
#include "plumbing/backend_gpu/defs.h" 
#endif

#include "datatypes/extended.h"


// declare MPI timers here too - these were externs

hila::timer start_send_timer("MPI start send");
hila::timer wait_send_timer("MPI wait send");
hila::timer post_receive_timer("MPI post receive");
hila::timer wait_receive_timer("MPI wait receive");
hila::timer synchronize_timer("MPI synchronize");
hila::timer reduction_timer("MPI reduction");
hila::timer reduction_wait_timer("MPI reduction wait");
hila::timer broadcast_timer("MPI broadcast");
hila::timer send_timer("MPI send field");
hila::timer drop_comms_timer("MPI wait drop_comms");
hila::timer partition_sync_timer("partition sync");

// let us house the partitions-struct here

hila::partitions_struct hila::partitions;

/* Keep track of whether MPI has been initialized */
static bool mpi_initialized = false;

////////////////////////////////////////////////////////////
/// Reductions: do automatic coalescing of reductions
/// if the type is float or double
/// These functions should not be called "by hand"

// buffers - first vector holds the reduction buffer,
// second the pointers to where distribute results
static std::vector<double> double_reduction_buffer;
static std::vector<double *> double_reduction_ptrs;
static int n_double = 0;

static std::vector<float> float_reduction_buffer;
static std::vector<float *> float_reduction_ptrs;
static int n_float = 0;

// static var holding the allreduce state
static bool allreduce_on = true;

void hila_reduce_double_setup(double *d, int n) {

    // ensure there's enough space
    if (n + n_double > double_reduction_buffer.size()) {
        double_reduction_buffer.resize(n + n_double + 2);
        double_reduction_ptrs.resize(n + n_double + 2);
    }

    for (int i = 0; i < n; i++) {
        double_reduction_buffer[n_double + i] = d[i];
        double_reduction_ptrs[n_double + i] = d + i;
    }

    n_double += n;
}

void hila_reduce_float_setup(float *d, int n) {

    // ensure there's enough space
    if (n + n_float > float_reduction_buffer.size()) {
        float_reduction_buffer.resize(n + n_float + 2);
        float_reduction_ptrs.resize(n + n_float + 2);
    }

    for (int i = 0; i < n; i++) {
        float_reduction_buffer[n_float + i] = d[i];
        float_reduction_ptrs[n_float + i] = d + i;
    }

    n_float += n;
}

void hila_reduce_sums() {

    if (n_double > 0) {
        std::vector<double> work(n_double);

        reduction_timer.start();

        if (allreduce_on) {
            MPI_Allreduce((void *)double_reduction_buffer.data(), (void *)work.data(), n_double,
                          MPI_DOUBLE, MPI_SUM, lattice->mpi_comm_lat);
            for (int i = 0; i < n_double; i++)
                *(double_reduction_ptrs[i]) = work[i];

        } else {
            MPI_Reduce((void *)double_reduction_buffer.data(), work.data(), n_double, MPI_DOUBLE,
                       MPI_SUM, 0, lattice->mpi_comm_lat);
            if (hila::myrank() == 0)
                for (int i = 0; i < n_double; i++)
                    *(double_reduction_ptrs[i]) = work[i];
        }

        n_double = 0;

        reduction_timer.stop();
    }

    if (n_float > 0) {
        std::vector<float> work(n_float);

        reduction_timer.start();

        if (allreduce_on) {
            MPI_Allreduce((void *)float_reduction_buffer.data(), work.data(), n_float, MPI_FLOAT,
                          MPI_SUM, lattice->mpi_comm_lat);
            for (int i = 0; i < n_float; i++)
                *(float_reduction_ptrs[i]) = work[i];

        } else {
            MPI_Reduce((void *)float_reduction_buffer.data(), work.data(), n_float, MPI_FLOAT,
                       MPI_SUM, 0, lattice->mpi_comm_lat);
            if (hila::myrank() == 0)
                for (int i = 0; i < n_float; i++)
                    *(float_reduction_ptrs[i]) = work[i];
        }

        n_float = 0;

        reduction_timer.stop();
    }
}

/// set allreduce on (default) or off on the next reduction
void hila::set_allreduce(bool on) {
    allreduce_on = on;
}

bool hila::get_allreduce() {
    return allreduce_on;
}

////////////////////////////////////////////////////////////////////////


/* Machine initialization */
#include <sys/types.h>
void hila::initialize_communications(int &argc, char ***argv) {
    /* Init MPI */
    if (!mpi_initialized) {

#ifndef OPENMP
        MPI_Init(&argc, argv);

#else

        int provided;
        MPI_Init_thread(&argc, argv, MPI_THREAD_FUNNELED, &provided);
        if (provided < MPI_THREAD_FUNNELED) {
            if (hila::myrank() == 0)
                hila::out << "MPI could not provide MPI_THREAD_FUNNELED, exiting\n";
            MPI_Finalize();
            exit(1);
        }

#endif



        mpi_initialized = true;

        // global var lattice exists, assign the mpi comms there
        lattice.ptr()->mpi_comm_lat = MPI_COMM_WORLD;

        MPI_Comm_rank(lattice->mpi_comm_lat, &lattice.ptr()->mynode.rank);
        MPI_Comm_size(lattice->mpi_comm_lat, &lattice.ptr()->nodes.number);

#ifdef GPU_CCL
        hila::initialize_gccl_communications();
#endif
    }
}

// check if MPI is on
bool hila::is_comm_initialized(void) {
    return mpi_initialized;
}

/* version of exit for multinode processes -- kill all nodes */
void hila::abort_communications(int status) {
    if (mpi_initialized) {
        mpi_initialized = false;
        MPI_Abort(lattice->mpi_comm_lat, 0);
    }
}

/* clean exit from all nodes */
void hila::finish_communications() {
    // turn off mpi -- this is needed to avoid mpi calls in destructors
    mpi_initialized = false;
    hila::about_to_finish = true;

    MPI_Finalize();
}

// broadcast specialization
void hila::broadcast(std::string &var, int rank) {

    if (hila::check_input)
        return;

    int size = var.size();
    hila::broadcast(size, rank);

    if (hila::myrank() != rank) {
        var.resize(size, ' ');
    }
    // copy directy to data() buffer
    broadcast_timer.start();
    MPI_Bcast((void *)var.data(), size, MPI_BYTE, rank, lattice->mpi_comm_lat);
    broadcast_timer.stop();
}

void hila::broadcast(std::vector<std::string> &list, int rank) {

    if (hila::check_input)
        return;

    int size = list.size();
    hila::broadcast(size, rank);
    list.resize(size);

    for (auto &s : list) {
        hila::broadcast(s, rank);
    }
}

/* BASIC COMMUNICATIONS FUNCTIONS */

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
int hila::number_of_nodes() {
    if (hila::check_input)
        return hila::check_with_nodes;

    int nodes;
    MPI_Comm_size(lattice->mpi_comm_lat, &nodes);
    return (nodes);
}

void hila::synchronize() {
    synchronize_timer.start();
    hila::synchronize_threads();
    MPI_Barrier(lattice->mpi_comm_lat);
    synchronize_timer.stop();
}

void hila::barrier() {
    synchronize_timer.start();
    MPI_Barrier(lattice->mpi_comm_lat);
    synchronize_timer.stop();
}


///  Get message tags cyclically -- defined outside classes, so that it is global and
///  unique

#define MSG_TAG_MIN 100
#define MSG_TAG_MAX (500) // standard says that at least 32767 tags available

int get_next_msg_tag() {
    static int tag = MSG_TAG_MIN;
    ++tag;
    if (tag > MSG_TAG_MAX)
        tag = MSG_TAG_MIN;
    return tag;
}


/// Split the communicator to subvolumes, using MPI_Comm_split
/// New MPI_Comm is the global mpi_comm_lat
/// NOTE: no attempt made here to reorder the nodes

void hila::split_into_partitions(int this_lattice) {

    if (hila::check_input)
        return;

    if (MPI_Comm_split(MPI_COMM_WORLD, this_lattice, 0, &(lattice.ptr()->mpi_comm_lat)) != MPI_SUCCESS) {
        hila::out0 << "MPI_Comm_split() call failed!\n";
        hila::finishrun();
    }
    // reset also the rank and numbers -fields
    MPI_Comm_rank(lattice->mpi_comm_lat, &lattice.ptr()->mynode.rank);
    MPI_Comm_size(lattice->mpi_comm_lat, &lattice.ptr()->nodes.number);
}

void hila::synchronize_partitions() {
    if (partitions.number() > 1)
        MPI_Barrier(MPI_COMM_WORLD);
}

MPI_Datatype MPI_ExtendedPrecision_type;
MPI_Op MPI_ExtendedPrecision_sum_op;

void create_extended_MPI_type() {
    ExtendedPrecision dummy;
    int block_lengths[3] = {1, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    MPI_Aint base;
    MPI_Get_address(&dummy, &base);
    MPI_Get_address(&dummy.value, &displacements[0]);
    MPI_Get_address(&dummy.compensation, &displacements[1]);
    MPI_Get_address(&dummy.compensation2, &displacements[2]);
    displacements[0] -= base;
    displacements[1] -= base;
    displacements[2] -= base;

    MPI_Type_create_struct(3, block_lengths, displacements, types, &MPI_ExtendedPrecision_type);
    MPI_Type_commit(&MPI_ExtendedPrecision_type);
}

void extended_sum_op(void *in, void *inout, int *len, MPI_Datatype *datatype) {
    ExtendedPrecision *in_data = (ExtendedPrecision *)in;
    ExtendedPrecision *inout_data = (ExtendedPrecision *)inout;

    for (int i = 0; i < *len; i++) {
        inout_data[i] += in_data[i];
    }
}

void create_extended_MPI_operation() {
    MPI_Op_create(&extended_sum_op, true, &MPI_ExtendedPrecision_sum_op);
}


/**
 * @brief Custom MPI reduction for extended type that performs Kahan summation
 *
 * @param value Input extended dataa
 * @param send_count Number of extended variables to reduce (default 1)
 * @param allreduce If true, performs MPI_Allreduce, otherwise MPI_Reduce
 */
void reduce_node_sum_extended(ExtendedPrecision *value, int send_count, bool allreduce) {

    if (hila::check_input)
        return;

    static bool init_extended_type_and_operation = true;
    if (init_extended_type_and_operation) {
        create_extended_MPI_type();
        create_extended_MPI_operation();
        init_extended_type_and_operation = false;
    }

    std::vector<ExtendedPrecision> recv_data(send_count);
    reduction_timer.start();
    if (allreduce) {
        MPI_Allreduce((void *)value, (void *)recv_data.data(), send_count,
                      MPI_ExtendedPrecision_type, MPI_ExtendedPrecision_sum_op,
                      lattice->mpi_comm_lat);
        for (int i = 0; i < send_count; i++)
            value[i] = recv_data[i];
    } else {
        MPI_Reduce((void *)value, (void *)recv_data.data(), send_count, MPI_ExtendedPrecision_type,
                   MPI_ExtendedPrecision_sum_op, 0, lattice->mpi_comm_lat);
        if (hila::myrank() == 0)
            for (int i = 0; i < send_count; i++)
                value[i] = recv_data[i];
    }
    reduction_timer.stop();
}
