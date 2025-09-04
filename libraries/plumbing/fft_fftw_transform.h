#ifndef FFTW_TRANSFORM_H
#define FFTW_TRANSFORM_H

/// Define hila_fft<>::transform() -functions and
/// scatter() / gather() -functions for fftw
///
/// This is not a standalone header, it is meant to be #include'd from
/// fft.h .

/// transform does the actual fft.
// template <typename cmplx_t>
// void hila_fft<cmplx_t>::transform() {
//     assert(0 && "Don't call this!");
// }

template <typename cmplx_t>
inline void hila_fft<cmplx_t>::transform() {

    static_assert(std::is_same<cmplx_t, Complex<double>>::value ||
                      std::is_same<cmplx_t, Complex<float>>::value,
                  "Only double or float fields in FFT");

    constexpr bool is_double = std::is_same<cmplx_t, Complex<double>>::value;

    // define fftw types for double and float
    using my_fftw_complex = typename std::conditional<is_double, fftw_complex, fftwf_complex>::type;
    using my_fftw_plan = typename std::conditional<is_double, fftw_plan, fftwf_plan>::type;

    extern unsigned hila_fft_my_columns[NDIM];
    extern hila::timer fft_plan_timer, fft_buffer_timer, fft_execute_timer;

    size_t n_fft = hila_fft_my_columns[dir] * elements;

    int transform_dir = (fftdir == fft_direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD;

    fft_plan_timer.start();

    // allocate here fftw plans.  TODO: perhaps store, if plans take appreciable time?
    // Timer will tell the proportional timing

    my_fftw_complex *fftwbuf;
    my_fftw_plan fftwplan;
    if constexpr (is_double) {
        fftwbuf = (my_fftw_complex *)fftw_malloc(sizeof(my_fftw_complex) * lattice.size(dir));
        fftwplan =
            fftw_plan_dft_1d(lattice.size(dir), fftwbuf, fftwbuf, transform_dir, FFTW_ESTIMATE);
    } else {
        fftwbuf = (my_fftw_complex *)fftwf_malloc(sizeof(my_fftw_complex) * lattice.size(dir));
        fftwplan =
            fftwf_plan_dft_1d(lattice.size(dir), fftwbuf, fftwbuf, transform_dir, FFTW_ESTIMATE);
    }

    fft_plan_timer.stop();

    for (size_t i = 0; i < n_fft; i++) {
        // collect stuff from buffers

        fft_buffer_timer.start();

        my_fftw_complex *cp = fftwbuf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(cp, rec_p[j] + i * rec_size[j], sizeof(my_fftw_complex) * rec_size[j]);
            cp += rec_size[j];
        }

        fft_buffer_timer.stop();

        // do the fft
        fft_execute_timer.start();

        if constexpr (is_double) {
            fftw_execute(fftwplan);
        } else {
            fftwf_execute(fftwplan);
        }

        fft_execute_timer.stop();

        fft_buffer_timer.start();

        cp = fftwbuf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(rec_p[j] + i * rec_size[j], cp, sizeof(my_fftw_complex) * rec_size[j]);
            cp += rec_size[j];
        }

        fft_buffer_timer.stop();
    }

    if constexpr (is_double) {
        fftw_destroy_plan(fftwplan);
        fftw_free(fftwbuf);
    } else {
        fftwf_destroy_plan(fftwplan);
        fftwf_free(fftwbuf);
    }
}

// template <>
// inline void hila_fft<Complex<float>>::transform() {

//     extern hila::timer fft_plan_timer, fft_buffer_timer, fft_execute_timer;
//     extern unsigned hila_fft_my_columns[NDIM];

//     size_t n_fft = hila_fft_my_columns[dir] * elements;

//     int transform_dir = (fftdir == fft_direction::forward) ? FFTW_FORWARD : FFTW_BACKWARD;

//     fft_plan_timer.start();

//     // allocate here fftw plans.  TODO: perhaps store, if plans take appreciable time?
//     // Timer will tell the proportional timing

//     fftwf_complex *fftwbuf =
//         (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * lattice.size(dir));
//     fftwf_plan fftwplan =
//         fftwf_plan_dft_1d(lattice.size(dir), fftwbuf, fftwbuf, transform_dir, FFTW_ESTIMATE);

//     fft_plan_timer.stop();

//     for (size_t i = 0; i < n_fft; i++) {
//         // collect stuff from buffers

//         fft_buffer_timer.start();

//         fftwf_complex *cp = fftwbuf;
//         for (int j = 0; j < rec_p.size(); j++) {
//             memcpy(cp, rec_p[j] + i * rec_size[j], sizeof(fftwf_complex) * rec_size[j]);
//             cp += rec_size[j];
//         }

//         fft_buffer_timer.stop();

//         // do the fft
//         fft_execute_timer.start();

//         fftwf_execute(fftwplan);

//         fft_execute_timer.stop();

//         fft_buffer_timer.start();

//         cp = fftwbuf;
//         for (int j = 0; j < rec_p.size(); j++) {
//             memcpy(rec_p[j] + i * rec_size[j], cp, sizeof(fftwf_complex) * rec_size[j]);
//             cp += rec_size[j];
//         }

//         fft_buffer_timer.stop();
//     }

//     fftwf_destroy_plan(fftwplan);
//     fftwf_free(fftwbuf);
// }

////////////////////////////////////////////////////////////////////
/// send column data to nodes

template <typename cmplx_t>
void hila_fft<cmplx_t>::gather_data() {


    extern hila::timer pencil_MPI_timer;
    pencil_MPI_timer.start();

    // post receive and send
    int n_comms = hila_pencil_comms[dir].size() - 1;

    std::vector<MPI_Request> sendreq(n_comms), recreq(n_comms);
    std::vector<MPI_Status> stat(n_comms);

    int i = 0;
    int j = 0;
    for (auto &fn : hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            size_t siz = fn.recv_buf_size * elements * sizeof(cmplx_t);
            if (siz >= (1ULL << 31)) {
                hila::out << "Too large MPI message in pencils! Size " << siz << " bytes\n";
                hila::terminate(1);
            }

            MPI_Irecv(rec_p[j], (int)siz, MPI_BYTE, fn.node, WRK_GATHER_TAG, lattice->mpi_comm_lat,
                      &recreq[i]);

            i++;
        }
        j++;
    }

    i = 0;
    for (auto &fn : hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            cmplx_t *p = send_buf + fn.column_offset * elements;
            int n = fn.column_number * elements * lattice->mynode.size[dir] * sizeof(cmplx_t);

            MPI_Isend(p, n, MPI_BYTE, fn.node, WRK_GATHER_TAG, lattice->mpi_comm_lat, &sendreq[i]);
            i++;
        }
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq.data(), stat.data());
        MPI_Waitall(n_comms, sendreq.data(), stat.data());
    }

    pencil_MPI_timer.stop();
}

//////////////////////////////////////////////////////////////////////////////////////
/// inverse of gather_data

template <typename cmplx_t>
void hila_fft<cmplx_t>::scatter_data() {


    extern hila::timer pencil_MPI_timer;
    pencil_MPI_timer.start();

    int n_comms = hila_pencil_comms[dir].size() - 1;

    std::vector<MPI_Request> sendreq(n_comms), recreq(n_comms);
    std::vector<MPI_Status> stat(n_comms);

    int i = 0;

    for (auto &fn : hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {
            cmplx_t *p = send_buf + fn.column_offset * elements;
            int n = fn.column_number * elements * lattice->mynode.size[dir] * sizeof(cmplx_t);

            MPI_Irecv(p, n, MPI_BYTE, fn.node, WRK_SCATTER_TAG, lattice->mpi_comm_lat, &recreq[i]);

            i++;
        }
    }

    i = 0;
    int j = 0;
    for (auto &fn : hila_pencil_comms[dir]) {
        if (fn.node != hila::myrank()) {

            MPI_Isend(rec_p[j], (int)(fn.recv_buf_size * elements * sizeof(cmplx_t)), MPI_BYTE,
                      fn.node, WRK_SCATTER_TAG, lattice->mpi_comm_lat, &sendreq[i]);

            i++;
        }
        j++;
    }

    // and wait for the send and receive to complete
    if (n_comms > 0) {
        MPI_Waitall(n_comms, recreq.data(), stat.data());
        MPI_Waitall(n_comms, sendreq.data(), stat.data());
    }

    pencil_MPI_timer.stop();
}

///////////////////////////////////////////////////////////////////////////////////
/// Separate reflect operation
/// Reflect flips the coordinates so that negative direction becomes positive,
/// and x=0 plane remains,
/// r(x) <- f(L - x)  -  note that x == 0 layer is as before
/// r(0) = f(0), r(1) = f(L-1), r(2) = f(L-2)  ...
///////////////////////////////////////////////////////////////////////////////////

template <typename cmplx_t>
inline void hila_fft<cmplx_t>::reflect() {
    extern unsigned hila_fft_my_columns[NDIM];
    extern hila::timer fft_plan_timer, fft_buffer_timer, fft_execute_timer;

    const int ncols = hila_fft_my_columns[dir] * elements;

    const int length = lattice.size(dir);

    cmplx_t *buf = (cmplx_t *)memalloc(sizeof(cmplx_t) * length);

    for (int i = 0; i < ncols; i++) {
        // collect stuff from buffers

        cmplx_t *cp = buf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(cp, rec_p[j] + i * rec_size[j], sizeof(cmplx_t) * rec_size[j]);
            cp += rec_size[j];
        }

        // reflect
        for (int j = 1; j < length / 2; j++) {
            std::swap(buf[j], buf[length - j]);
        }

        cp = buf;
        for (int j = 0; j < rec_p.size(); j++) {
            memcpy(rec_p[j] + i * rec_size[j], cp, sizeof(cmplx_t) * rec_size[j]);
            cp += rec_size[j];
        }
    }

    free(buf);
}


#endif