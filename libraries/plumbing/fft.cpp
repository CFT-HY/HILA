/// Static variables and functions for hila fft routines


#include "plumbing/defs.h"
#include "plumbing/timing.h"
#include "plumbing/fft.h"


hila::timer fft_timer("FFT total time");
hila::timer fft_plan_timer(" FFT plan");
hila::timer pencil_MPI_timer(" MPI for pencils");
hila::timer fft_execute_timer(" FFT execute");
hila::timer pencil_reshuffle_timer(" pencil reshuffle");
hila::timer fft_buffer_timer(" copy fft buffers");
hila::timer pencil_collect_timer(" copy pencils");
hila::timer pencil_save_timer(" save pencils");

hila::timer binning_timer("bin field time");

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

    foralldir (d)
        if (d != dir) {
            offset[d] = s;
            s *= lattice->mynode.size[d];
        }

    return element_offset;
}

/// THis is to be called before fft to Direction dir


/////////////////////////////////////////////////////////////////////////////////////

void init_pencil_direction(Direction dir) {

    if (lattice->fftdata == nullptr) {
        lattice.ptr()->fftdata = new hila::fftdata_struct;
    }
    
    hila::fftdata_struct &fft = *(lattice.ptr()->fftdata);

    if (fft.hila_pencil_comms[dir].size() == 0) {
        // basic structs not yet set, do it here

        fft.hila_pencil_comms[dir].resize(lattice->nodes.n_divisions[dir]);

        int nodenumber = 0;
        for (int nodenumber = 0; nodenumber < lattice->nodes.number; ++nodenumber) {

            const node_info n = lattice->nodes.nodeinfo(nodenumber);

            bool is_in_column = true;
            foralldir (d)
                if (d != dir && n.min[d] != lattice->mynode.min[d]) {
                    is_in_column = false;
                    break; // breaks out of foralldir
                }

            // store the nodes in the hila_pencil_comms -list in the right order -
            // nodes may be reordered by some weird layout
            if (is_in_column) {
                hila::pencil_struct fn;
                fn.node = nodenumber;
                fn.size_to_dir = n.size[dir];
                for (int i = 0; i < lattice->nodes.n_divisions[dir]; i++) {
                    if (n.min[dir] == lattice->nodes.divisors[dir][i]) {
                        fft.hila_pencil_comms[dir].at(i) = fn;
                    }
                }
            }
        }

        size_t total_columns = lattice->mynode.volume / lattice->mynode.size[dir];

        size_t nodes = fft.hila_pencil_comms[dir].size();

        // column offset and number are used for sending
        size_t i = 0;
        for (auto &fn : fft.hila_pencil_comms[dir]) {
            fn.column_offset = ((i * total_columns) / nodes) * lattice->mynode.size[dir];
            fn.column_number =
                (((i + 1) * total_columns) / nodes) - fn.column_offset / lattice->mynode.size[dir];

            if (fn.node == hila::myrank()) {
                fft.hila_fft_my_columns[dir] = fn.column_number;
            }
            i++;
        }

        fft.pencil_recv_buf_size[dir] = 0;

        for (auto &fn : fft.hila_pencil_comms[dir]) {
            fn.recv_buf_size = fft.hila_fft_my_columns[dir] * fn.size_to_dir;

            // how big array?
            if (fn.node != hila::myrank())
                fft.pencil_recv_buf_size[dir] += fn.recv_buf_size;
        }

    } // setup
}
