#ifndef FFTW_TRANSFORM_H
#define FFTW_TRANSFORM_H

/// Define hila_fft<>::transform() -functions and
/// scatter() / gather() -functions for fftw
///
/// This is not a standalone header, it is meant to be #include'd from
/// fft.h .

/// transform does the actual fft.
template <typename cmplx_t> void hila_fft<cmplx_t>::transform() {
    assert(0 && "Don't call this!");
}

template <> inline void hila_fft<Complex<double>>::transform() {
    extern unsigned hila_fft_my_columns[NDIM];
    extern hila::timer fft_plan_timer, fft_buffer_timer, fft_execute_timer;

    size_t n_fft = hila_fft_my_columns[dir] * elements;

    int transform_dir =
        (fftdir == fft_direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD;

    fft_plan_timer.start();

    // allocate here fftw plans.  TODO: perhaps store, if plans take appreciable time?
    // Timer will tell the proportional timing

    fftw_complex *fftwbuf =
        (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * lattice->size(dir));
    fftw_plan fftwplan = fftw_plan_dft_1d(lattice->size(dir), fftwbuf, fftwbuf,
                                          transform_dir, FFTW_ESTIMATE);

    fft_plan_timer.stop();

    for (size_t i = 0; i < n_fft; i++) {
        // collect stuff from buffers

        fft_buffer_timer.start();

        fftw_complex *cp = fftwbuf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(cp, rec_p[j] + i * rec_size[j], sizeof(fftw_complex) * rec_size[j]);
            cp += rec_size[j];
        }

        fft_buffer_timer.stop();

        // do the fft
        fft_execute_timer.start();

        fftw_execute(fftwplan);

        fft_execute_timer.stop();

        fft_buffer_timer.start();

        cp = fftwbuf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(rec_p[j] + i * rec_size[j], cp, sizeof(fftw_complex) * rec_size[j]);
            cp += rec_size[j];
        }

        fft_buffer_timer.stop();
    }

    fftw_destroy_plan(fftwplan);
    fftw_free(fftwbuf);
}

template <> inline void hila_fft<Complex<float>>::transform() {

    extern hila::timer fft_plan_timer, fft_buffer_timer, fft_execute_timer;
    extern unsigned hila_fft_my_columns[NDIM];

    size_t n_fft = hila_fft_my_columns[dir] * elements;

    int transform_dir =
        (fftdir == fft_direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD;

    fft_plan_timer.start();

    // allocate here fftw plans.  TODO: perhaps store, if plans take appreciable time?
    // Timer will tell the proportional timing

    fftwf_complex *fftwbuf =
        (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * lattice->size(dir));
    fftwf_plan fftwplan = fftwf_plan_dft_1d(lattice->size(dir), fftwbuf, fftwbuf,
                                            transform_dir, FFTW_ESTIMATE);

    fft_plan_timer.stop();

    for (size_t i = 0; i < n_fft; i++) {
        // collect stuff from buffers

        fft_buffer_timer.start();

        fftwf_complex *cp = fftwbuf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(cp, rec_p[j] + i * rec_size[j], sizeof(fftwf_complex) * rec_size[j]);
            cp += rec_size[j];
        }

        fft_buffer_timer.stop();

        // do the fft
        fft_execute_timer.start();

        fftwf_execute(fftwplan);

        fft_execute_timer.stop();

        fft_buffer_timer.start();

        cp = fftwbuf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(rec_p[j] + i * rec_size[j], cp, sizeof(fftwf_complex) * rec_size[j]);
            cp += rec_size[j];
        }

        fft_buffer_timer.stop();
    }

    fftwf_destroy_plan(fftwplan);
    fftwf_free(fftwbuf);
}

////////////////////////////////////////////////////////////////////
/// send column data to nodes

template <typename cmplx_t> void hila_fft<cmplx_t>::gather_data() {

#ifdef USE_MPI

    extern hila::timer fft_MPI_timer;
    fft_MPI_timer.start();

    // post receive and send
    int n_comms = hila_fft_comms[dir].size() - 1;

    MPI_Request sendreq[n_comms], recreq[n_comms];
    MPI_Status stat[n_comms];

    int i = 0;
    int j = 0;
    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            size_t siz = fn.recv_buf_size * elements;
            if (siz >= (1ULL << 30)) {
                hila::output << "Too large MPI message in FFT! Size " << siz
                             << " complex numbers\n";
                hila::terminate(1);
            }

            MPI_Irecv(rec_p[j], (int)siz, mpi_cmplx_t, fn.node, WRK_GATHER_TAG,
                      lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
        j++;
    }

    i = 0;
    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            cmplx_t *p = send_buf + fn.column_offset * elements;
            int n = fn.column_number * elements * lattice->mynode.size[dir];

            MPI_Isend(p, n, mpi_cmplx_t, fn.node, WRK_GATHER_TAG, lattice->mpi_comm_lat,
                      &sendreq[i]);
            i++;
        }
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq, stat);
        MPI_Waitall(n_comms, sendreq, stat);
    }

    fft_MPI_timer.stop();

#endif
}

//////////////////////////////////////////////////////////////////////////////////////
/// inverse of gather_data

template <typename cmplx_t> void hila_fft<cmplx_t>::scatter_data() {

#ifdef USE_MPI

    extern hila::timer fft_MPI_timer;
    fft_MPI_timer.start();

    int n_comms = hila_fft_comms[dir].size() - 1;

    MPI_Request sendreq[n_comms], recreq[n_comms];
    MPI_Status stat[n_comms];

    int i = 0;

    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {
            cmplx_t *p = send_buf + fn.column_offset * elements;
            int n = fn.column_number * elements * lattice->mynode.size[dir];

            MPI_Irecv(p, n, mpi_cmplx_t, fn.node, WRK_SCATTER_TAG,
                      lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
    }

    i = 0;
    int j = 0;
    for (auto &fn : hila_fft_comms[dir]) {
        if (fn.node != hila::myrank()) {

            MPI_Isend(rec_p[j], fn.recv_buf_size * elements, mpi_cmplx_t, fn.node,
                      WRK_SCATTER_TAG, lattice->mpi_comm_lat, &sendreq[i]);

            i++;
        }
        j++;
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq, stat);
        MPI_Waitall(n_comms, sendreq, stat);
    }

    fft_MPI_timer.stop();
#endif
}

#endif