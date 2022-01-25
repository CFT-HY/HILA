#ifndef HILAPP_MPI_H_
#define HILAPP_MPI_H_


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
    MPI_UINT64_T
};

enum MPI_Op : int { MPI_SUM, MPI_PROD, MPI_MAX, MPI_MIN };

typedef void *MPI_Comm;
typedef void *MPI_Request;
typedef int MPI_Status;
typedef void *MPI_Comm;
#define MPI_IN_PLACE nullptr
#define MPI_COMM_WORLD nullptr
#define MPI_SUCCESS 1

enum MPI_thread_level : int {
    MPI_THREAD_SINGLE,
    MPI_THREAD_FUNNELED,
    MPI_THREAD_SERIALIZED,
    MPI_THREAD_MULTIPLE
};


int MPI_Init(int *argc, char ***argv);

int MPI_Init_thread(int *argc, char ***argv, int threadlevel, int *provided);

int MPI_Comm_rank(MPI_Comm comm, int *rank);

int MPI_Comm_size(MPI_Comm comm, int *size);

int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);

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

int MPI_Abort(MPI_Comm comm, int errorcode);

int MPI_Finalize();

#endif