#ifndef HILA_HILAPP_MPI_H_
#define HILA_HILAPP_MPI_H_


////////////////////////////////////////////////////////////////////////////
/// This header introduces MPI functions for hilapp, not for final compile
/// Useful because (statically compiled) hilapp does not automatically know
/// the location of mpi.h, and saves the trouble of finding it
////////////////////////////////////////////////////////////////////////////


// Selection of MPI datatypes - if needed add here
enum MPI_Datatype : int {
    MPI_BYTE,
    MPI_CHAR,
    MPI_SHORT,
    MPI_INT,
    MPI_LONG,
    MPI_FLOAT,
    MPI_DOUBLE,
    MPI_LONG_DOUBLE,
    MPI_C_DOUBLE_COMPLEX,
    MPI_C_FLOAT_COMPLEX,
    MPI_UNSIGNED,
    MPI_UNSIGNED_LONG,
    MPI_INT64_T,
    MPI_UINT64_T,
    MPI_2INT,
    MPI_LONG_INT,
    MPI_FLOAT_INT,
    MPI_DOUBLE_INT,
    MPI_LONG_DOUBLE_INT
};

enum MPI_Op : int { MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN, MPI_MAXLOC, MPI_MINLOC };

typedef void *MPI_Comm;
typedef void *MPI_Request;
typedef struct ompi_status_public_t MPI_Status;
typedef void *MPI_Comm;
typedef int MPI_Fint;
typedef int MPI_Aint;
typedef void *MPI_Errhandler;
#define MPI_IN_PLACE nullptr
#define MPI_COMM_WORLD nullptr
#define MPI_STATUS_IGNORE nullptr
#define MPI_ERRORS_RETURN nullptr
#define MPI_REQUEST_NULL nullptr
#define MPI_SUCCESS 1

enum MPI_thread_level : int {
    MPI_THREAD_SINGLE,
    MPI_THREAD_FUNNELED,
    MPI_THREAD_SERIALIZED,
    MPI_THREAD_MULTIPLE
};

struct ompi_status_public_t {
    /* These fields are publicly defined in the MPI specification.
       User applications may freely read from these fields. */
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
    /* The following two fields are internal to the Open MPI
       implementation and should not be accessed by MPI applications.
       They are subject to change at any time.  These are not the
       droids you're looking for. */
    int _cancelled;
    size_t _ucount;
};

int MPI_Init(int *argc, char ***argv);

int MPI_Init_thread(int *argc, char ***argv, int threadlevel, int *provided);

int MPI_Comm_rank(MPI_Comm comm, int *rank);

int MPI_Comm_size(MPI_Comm comm, int *size);

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);

int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);

int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm);

int MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, int root, MPI_Comm comm, MPI_Request *request);

int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm);

int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm, MPI_Request *request);

int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm);

int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request);

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status);

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request);

int MPI_Wait(MPI_Request *request, MPI_Status *status);

int MPI_Waitall(int count, MPI_Request array_of_requests[],
                MPI_Status *array_of_statuses);

int MPI_Barrier(MPI_Comm comm);

int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request);

int MPI_Cancel(MPI_Request *request);

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);

int MPI_Test_cancelled(const MPI_Status *status, int *flag);

int MPI_Request_free(MPI_Request *request);

int MPI_Abort(MPI_Comm comm, int errorcode);

MPI_Fint MPI_Comm_c2f(MPI_Comm comm);

int MPI_Finalize();

int MPI_Get_address(const void *location, MPI_Aint *address);

int MPI_Type_create_struct(int count,
                           const int array_of_blocklengths[],
                           const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[],
                           MPI_Datatype *newtype);

int MPI_Type_commit(MPI_Datatype *datatype);

typedef void MPI_User_function(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype);

int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op);


#endif