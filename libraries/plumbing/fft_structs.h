#ifndef HILA_FFT_STRUCTS_H
#define HILA_FFT_STRUCTS_H

namespace hila {

// hold static fft node data structures
struct pencil_struct {
    int node;               // node rank to send stuff for fft:ing
    unsigned size_to_dir;   // size of "node" to fft-dir
    unsigned column_offset; // first perp-plane column to be handled by "node"
    unsigned column_number; // and number of columns to be sent
    size_t recv_buf_size;   // size of my fft collect buffer (in units of sizeof(T)
                            // for stuff received from / returned to "node"
};


struct fftdata_struct {
    std::vector<pencil_struct> hila_pencil_comms[NDIM];
    unsigned hila_fft_my_columns[NDIM]; // how many columns does this node take care of
    size_t pencil_recv_buf_size[NDIM];
};

} // namespace hila


#endif