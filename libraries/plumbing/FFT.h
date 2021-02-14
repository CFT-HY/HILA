#ifndef FFT_H
#define FFT_H

#include "plumbing/field.h"
#include "datatypes/cmplx.h"
#include "plumbing/timing.h"
#include "fftw3.h"

#ifdef USE_MPI
/// Run Fast Fourier Transform on the field to each direction
/// This is done by collecting a column of elements to each node,
/// running the Fourier transform on the column and redistributing
/// the result
/// Input and result are passed by reference. They may be the same.
/// The field must be complex and the underlying complex type is supplied
/// by the complex_type template argument

static hila::timer FFT_timer("FFT"), FFT_MPI_timer(" MPI in FFT"); // initialized 1st time used
static hila::timer fftw_execute_timer("FFTW execute"), fftw_plan_timer("FFTW plan");
static hila::timer sitelist_timer("fft sitelist");
static hila::timer fft_copy_payload_timer("copy payload");
static hila::timer fft_buf_timer("copy fftw buffers");
static hila::timer fft_place_timer("place payload");

template <typename T, typename complex_type>
inline void FFT_field_complex(Field<T> &input, Field<T> &result,
                              fft_direction fftdir = fft_direction::forward) {

    lattice_struct *lattice = input.fs->lattice;
    Field<T> *read_pointer = &input; // Read from input on first time, then work in result
    size_t local_volume = lattice->mynode.volume();
    int elements = sizeof(T) / sizeof(complex_type);

    FFT_timer.start();

    // Make store the result is allocated and mark it changed
    result.check_alloc();
    result.mark_changed(ALL);

    // Allocate buffers for the MPI
    char *mpi_send_buffer = (char *)malloc(local_volume * sizeof(T));
    char *mpi_recv_buffer = (char *)malloc(local_volume * sizeof(T));

    // Run transform in all directions
    foralldir(dir) {
        // Get the number of sites per column on this node and on all nodes
        size_t node_column_size = lattice->mynode.size[dir];
        size_t column_size = lattice->size(dir);
        // Count columns on this rank
        int cols = 1;
        foralldir(d2) if (d2 != dir) cols *= lattice->mynode.size[d2];

        // Get the MPI column in this direction
        lattice_struct::mpi_column_struct mpi_column = lattice->get_mpi_column(dir);
        MPI_Comm column_communicator = mpi_column.column_communicator;
        const std::vector<int> &nodelist = mpi_column.nodelist;
        int my_column_rank = mpi_column.my_column_rank;
        int nnodes = nodelist.size(); // Nodes in this column
        int cpn = cols / nnodes;      // Columns per node

        assert(cols % nnodes == 0 &&
               "This FFT requires columns per node is evenly divisible");
        // Amount of data that is communicated from a single node to another
        int block_size = cpn * node_column_size * sizeof(T);
        // Size of a full column
        int col_size = node_column_size * sizeof(T);

        // Variables needed for constructing the columns of sites
        const std::vector<node_info> &allnodes = lattice->nodelist();
        CoordinateVector min = allnodes[lattice->node_rank()].min;
        CoordinateVector size = allnodes[lattice->node_rank()].size;

        // FFTW buffers
        fftw_complex *in, *out;
        in = (fftw_complex * RESTRICT) fftw_malloc(sizeof(fftw_complex) * column_size);
        out = (fftw_complex * RESTRICT) fftw_malloc(sizeof(fftw_complex) * column_size);

        int fdir = (fftdir == fft_direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD;
        fftw_plan_timer.start();
        fftw_plan plan = fftw_plan_dft_1d(column_size, in, out, fdir, FFTW_ESTIMATE);
        fftw_plan_timer.stop();

        sitelist_timer.start();
        // Construct lists of sites for each column
        std::vector<std::vector<unsigned>> sitelist(cols);
        for (int c = 0; c < cols; c++) {
            // Get the index in the other directions
            int cc = c;
            CoordinateVector thiscol = min;
            foralldir(d2) if (d2 != dir) {
                thiscol[d2] += cc % size[d2];
                cc /= size[d2];
            }

            // And build the list of sites
            sitelist[c].resize(node_column_size);
            for (int i = 0; i < node_column_size; i++) {
                thiscol[dir] = min[dir] + i;
                sitelist[c][i] = lattice->site_index(thiscol);
            }
        }
        sitelist_timer.stop();

        // Gather a number of columns to each node
        MPI_Request my_request;
        MPI_Request other_reqs[nnodes];
        int ireq = 0;
        for (int r = 0; r < nnodes; r++) {
            fft_copy_payload_timer.start();
            char *sendbuf = mpi_send_buffer + block_size * r;
            for (int l = 0; l < cpn; l++) {
                int c = r * cpn + l;
                read_pointer->fs->payload.gather_elements((T *)(sendbuf + col_size * l),
                                                          sitelist[c].data(),
                                                          sitelist[c].size(), lattice);
            }
            fft_copy_payload_timer.stop();
            FFT_MPI_timer.start();
            MPI_Request request;
            MPI_Igather(sendbuf, cpn * node_column_size * sizeof(T), MPI_BYTE,
                        mpi_recv_buffer, cpn * node_column_size * sizeof(T), MPI_BYTE, r,
                        column_communicator, &request);
            if (r == my_column_rank) {
                my_request = request;
            } else {
                other_reqs[ireq++] = request;
            }
            FFT_MPI_timer.stop();
        }

        // Wait for my data
        FFT_MPI_timer.start();
        MPI_Status status;
        MPI_Wait(&my_request, &status);
        FFT_MPI_timer.stop();

        // now that we have columns, run FFT on each
        for (int l = 0; l < cpn; l++) {          // Columns
            for (int e = 0; e < elements; e++) { // Complex elements / field element

                fft_buf_timer.start();
                for (int s = 0; s < nnodes;
                     s++) { // Cycle over sender nodes to collect the data
                    complex_type *RESTRICT field_elem =
                        (complex_type *)(mpi_recv_buffer + block_size * s + col_size * l);
                    for (int t = 0; t < node_column_size; t++) {
                        in[t + node_column_size * s][0] = field_elem[e + elements * t].re;
                        in[t + node_column_size * s][1] = field_elem[e + elements * t].im;
                    }
                }
                fft_buf_timer.stop();

                fftw_execute_timer.start();
                // Run the fft
                fftw_execute(plan);
                fftw_execute_timer.stop();

                fft_buf_timer.start();
                // Put the transformed data back in place
                for (int s = 0; s < nnodes; s++) {
                    complex_type *RESTRICT field_elem =
                        (complex_type *)(mpi_recv_buffer + block_size * s + col_size * l);
                    for (int t = 0; t < node_column_size; t++) {
                        field_elem[e + elements * t].re =
                            out[t + node_column_size * s][0];
                        field_elem[e + elements * t].im =
                            out[t + node_column_size * s][1];
                    }
                }
                fft_buf_timer.stop();
            }
        }

        // clear the other requests -- wonder if these are needed?
        if (nnodes > 1) {
            FFT_MPI_timer.start();
            MPI_Status statarr[nnodes - 1];
            MPI_Waitall(nnodes - 1, other_reqs, statarr);
            FFT_MPI_timer.stop();
        }

        // Now reverse the gather operation. After this each node will have its original
        // local sites The scatter operation cannot be asynchronous, but this coul dbe
        // implemented with a gather operations instead.
        for (int s = 0; s < nnodes; s++) {
            FFT_MPI_timer.start();
            char *sendbuf = mpi_send_buffer + block_size * s;
            MPI_Scatter(mpi_recv_buffer, cpn * node_column_size * sizeof(T), MPI_BYTE,
                        sendbuf, cpn * node_column_size * sizeof(T), MPI_BYTE, s,
                        column_communicator);

            FFT_MPI_timer.stop();

            fft_place_timer.start();
            // Place the new data into field memory
            for (int l = 0; l < cpn; l++) {
                int c = s * cpn + l;
                result.fs->payload.place_elements((T *)(sendbuf + col_size * l),
                                                  sitelist[c].data(), sitelist[c].size(),
                                                  lattice);
            }
            fft_place_timer.stop();
        }

        read_pointer = &result; // From now on we work in result
        fftw_destroy_plan(plan);
        fftw_free(in);
        fftw_free(out);
    }

    free(mpi_send_buffer);
    free(mpi_recv_buffer);

    FFT_timer.stop();
}

#else // No mpi

template <typename T, typename C>
inline void FFT_field_complex(Field<T> &input, Field<T> &result,
                              fft_direction fdir = fft_direction::forward) {}

#endif

template <> inline void Field<Complex<double>>::FFT(fft_direction fdir) {
    FFT_field_complex<Complex<double>, Complex<double>>(*this, *this, fdir);
}

/// Match a given type T to it's underlying complex type
template <typename T, class Enable = void> struct complex_base {};

/// Match to a complex type
template <> struct complex_base<Complex<float>> { using type = Complex<float>; };

/// Match to a complex type
template <> struct complex_base<Complex<double>> { using type = Complex<double>; };

/// Match templated class B to it's underlying complex type
template <template <typename B> class C, typename B> struct complex_base<C<B>> {
    using type = typename complex_base<B>::type;
};

/// Match templated class B to it's underlying complex type
template <template <int a, typename B> class C, int a, typename B>
struct complex_base<C<a, B>> {
    using type = typename complex_base<B>::type;
};

/// Match templated class B to it's underlying complex type
template <template <int a, int b, typename B> class C, int a, int b, typename B>
struct complex_base<C<a, b, B>> {
    using type = typename complex_base<B>::type;
};

/// Run fourier transform on a complex field
// Called with any type T with a Complex type nested in the lowest level
template <typename T>
void FFT_field(Field<T> &input, Field<T> &result,
               fft_direction fdir = fft_direction::forward) {
    FFT_field_complex<T, typename complex_base<T>::type>(input, result, fdir);
}

#endif