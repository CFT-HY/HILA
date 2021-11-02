/// Static variables and functions for hila fft routines


#include "plumbing/defs.h"
#include "plumbing/timing.h"
#include "plumbing/fft.h"


#ifdef USE_MPI
#include <mpi.h>
#endif


hila::timer fft_timer("FFT total time");
hila::timer fft_plan_timer("  FFT plan");
hila::timer fft_MPI_timer("  MPI in FFT");
hila::timer fft_execute_timer("  FFT execute");
hila::timer fft_reshuffle_timer("  data reshuffle");
hila::timer fft_buffer_timer("  copy fft buffers");
hila::timer fft_collect_timer("  copy payload");
hila::timer fft_save_timer("  save result");

// a couple of global fft variables - easiest this way
std::vector<fftnode_struct> hila_fft_comms[NDIM];
unsigned hila_fft_my_columns[NDIM]; // how many columns does this node take care of
size_t fft_recv_buf_size[NDIM];

#if defined(CUDA) || defined(HIP)
// save plan 

#endif


size_t fft_get_buffer_offsets(const Direction dir, const size_t elements,
                              CoordinateVector &offset, CoordinateVector &nmin) {

    offset[dir] = 1;
    nmin = lattice->mynode.min;

    size_t element_offset = lattice->mynode.size[dir];
    size_t s = element_offset * elements;

    foralldir(d) if (d != dir) {
        offset[d] = s;
        s *= lattice->mynode.size[d];
    }

    return element_offset;
}

/// THis is to be called before fft to Direction dir


/////////////////////////////////////////////////////////////////////////////////////

void init_fft_direction(Direction dir) {

    if (hila_fft_comms[dir].size() == 0) {
        // basic structs not yet set, do it here

        hila_fft_comms[dir].resize(lattice->nodes.n_divisions[dir]);

        int nodenumber = 0;
        for (const node_info &n : lattice->nodes.nodelist) {
            bool is_in_column = true;
            foralldir(d) if (d != dir && n.min[d] != lattice->mynode.min[d]) {
                is_in_column = false;
                break; // breaks out of foralldir
            }

            // store the nodes in the hila_fft_comms -list in the right order -
            // nodes may be reordered by some weird layout
            if (is_in_column) {
                fftnode_struct fn;
                fn.node = nodenumber;
                fn.size_to_dir = n.size[dir];
                for (int i = 0; i < lattice->nodes.n_divisions[dir]; i++) {
                    if (n.min[dir] == lattice->nodes.divisors[dir][i]) {
                        hila_fft_comms[dir].at(i) = fn;
                    }
                }
            }
            ++nodenumber;
        }

        size_t total_columns = lattice->mynode.sites / lattice->mynode.size[dir];

        size_t nodes = hila_fft_comms[dir].size();

        // column offset and number are used for sending
        size_t i = 0;
        for (fftnode_struct &fn : hila_fft_comms[dir]) {
            fn.column_offset =
                ((i * total_columns) / nodes) * lattice->mynode.size[dir];
            fn.column_number = (((i + 1) * total_columns) / nodes) -
                               fn.column_offset / lattice->mynode.size[dir];

            if (fn.node == hila::myrank()) {
                hila_fft_my_columns[dir] = fn.column_number;
            }
            i++;
        }

        fft_recv_buf_size[dir] = 0;

        for (fftnode_struct &fn : hila_fft_comms[dir]) {
            fn.recv_buf_size = hila_fft_my_columns[dir] * fn.size_to_dir;

            // how big array?
            if (fn.node != hila::myrank())
                fft_recv_buf_size[dir] += fn.recv_buf_size;
        }

    } // setup
}
