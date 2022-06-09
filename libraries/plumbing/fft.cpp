/// Static variables and functions for hila fft routines


#include "plumbing/defs.h"
#include "plumbing/timing.h"
#include "plumbing/fft.h"


hila::timer fft_timer("FFT total time");
hila::timer fft_plan_timer("  FFT plan");
hila::timer pencil_MPI_timer("  MPI in FFT");
hila::timer fft_execute_timer("  FFT execute");
hila::timer pencil_reshuffle_timer("  data reshuffle");
hila::timer fft_buffer_timer("  copy fft buffers");
hila::timer pencil_collect_timer("  copy payload");
hila::timer pencil_save_timer("  save result");

// a couple of global fft variables - easiest this way
std::vector<pencil_struct> hila_pencil_comms[NDIM];
unsigned hila_fft_my_columns[NDIM]; // how many columns does this node take care of
size_t pencil_recv_buf_size[NDIM];

// static variable to hold fft plans
#if (defined(HIP) || defined(CUDA)) && !defined(HILAPP)
hila_saved_fftplan_t hila_saved_fftplan;
#endif

// Delete saved plans if required - no-op for non-gpu
void FFT_delete_plans() {
#if (defined(HIP) || defined(CUDA)) && !defined(HILAPP)
    hila_saved_fftplan.delete_plans();
#endif
}



size_t pencil_get_buffer_offsets(const Direction dir, const size_t elements,
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

void init_pencil_direction(Direction dir) {

    if (hila_pencil_comms[dir].size() == 0) {
        // basic structs not yet set, do it here

        hila_pencil_comms[dir].resize(lattice->nodes.n_divisions[dir]);

        int nodenumber = 0;
        for (const node_info &n : lattice->nodes.nodelist) {
            bool is_in_column = true;
            foralldir(d) if (d != dir && n.min[d] != lattice->mynode.min[d]) {
                is_in_column = false;
                break; // breaks out of foralldir
            }

            // store the nodes in the hila_pencil_comms -list in the right order -
            // nodes may be reordered by some weird layout
            if (is_in_column) {
                pencil_struct fn;
                fn.node = nodenumber;
                fn.size_to_dir = n.size[dir];
                for (int i = 0; i < lattice->nodes.n_divisions[dir]; i++) {
                    if (n.min[dir] == lattice->nodes.divisors[dir][i]) {
                        hila_pencil_comms[dir].at(i) = fn;
                    }
                }
            }
            ++nodenumber;
        }

        size_t total_columns = lattice->mynode.sites / lattice->mynode.size[dir];

        size_t nodes = hila_pencil_comms[dir].size();

        // column offset and number are used for sending
        size_t i = 0;
        for (pencil_struct &fn : hila_pencil_comms[dir]) {
            fn.column_offset =
                ((i * total_columns) / nodes) * lattice->mynode.size[dir];
            fn.column_number = (((i + 1) * total_columns) / nodes) -
                               fn.column_offset / lattice->mynode.size[dir];

            if (fn.node == hila::myrank()) {
                hila_fft_my_columns[dir] = fn.column_number;
            }
            i++;
        }

        pencil_recv_buf_size[dir] = 0;

        for (pencil_struct &fn : hila_pencil_comms[dir]) {
            fn.recv_buf_size = hila_fft_my_columns[dir] * fn.size_to_dir;

            // how big array?
            if (fn.node != hila::myrank())
                pencil_recv_buf_size[dir] += fn.recv_buf_size;
        }

    } // setup
}
