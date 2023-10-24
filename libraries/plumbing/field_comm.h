/** @file field_comm.h */

#ifndef FIELD_COMM_H
#define FIELD_COMM_H

#include "field.h"

/////////////////////////////////////////////////////////////////////////////////
/// This file collects implementations of field communication routines


/////////////////////////////////////////////////////////////////////////////////////////
/// Gather boundary elements for communication

template <typename T>
void Field<T>::field_struct::gather_comm_elements(
    Direction d, Parity par, T *RESTRICT buffer,
    const lattice_struct::comm_node_struct &to_node) const {
#ifndef VECTORIZED
#ifdef SPECIAL_BOUNDARY_CONDITIONS
    // note: -d in is_on_edge, because we're about to send stuff to that
    // Direction (gathering from Direction +d)
    if (boundary_condition[d] == hila::bc::ANTIPERIODIC && lattice.mynode.is_on_edge(-d)) {
        payload.gather_comm_elements(buffer, to_node, par, lattice, true);
    } else {
        payload.gather_comm_elements(buffer, to_node, par, lattice, false);
    }
#else
    payload.gather_comm_elements(buffer, to_node, par, lattice, false);
#endif

#else
    // this is vectorized branch
    bool antiperiodic = false;
#ifdef SPECIAL_BOUNDARY_CONDITIONS
    if (boundary_condition[d] == hila::bc::ANTIPERIODIC && lattice.mynode.is_on_edge(-d)) {
        antiperiodic = true;
    }
#endif

    if constexpr (hila::is_vectorizable_type<T>::value) {
        // now vectorized layout
        if (vector_lattice->is_boundary_permutation[abs(d)]) {
            // with boundary permutation need to gather elems 1-by-1
            int n;
            const unsigned *index_list = to_node.get_sitelist(par, n);
            if (!antiperiodic) {
                payload.gather_elements(buffer, index_list, n, lattice);
            } else {
                payload.gather_elements_negated(buffer, index_list, n, lattice);
            }
        } else {
            // without it, can do the full block
            payload.gather_comm_vectors(buffer, to_node, par, vector_lattice, antiperiodic);
        }
    } else {
        // not vectoizable, standard methods
        int n;
        const unsigned *index_list = to_node.get_sitelist(par, n);
        if (!antiperiodic)
            payload.gather_elements(buffer, index_list, n, lattice);
        else {
            payload.gather_elements_negated(buffer, index_list, n, lattice);
        }
    }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Place boundary elements from neighbour

template <typename T>
void Field<T>::field_struct::place_comm_elements(
    Direction d, Parity par, T *RESTRICT buffer,
    const lattice_struct::comm_node_struct &from_node) {

#ifdef VECTORIZED
    if constexpr (hila::is_vectorizable_type<T>::value) {
        // now vectorized layout, act accordingly
        if (vector_lattice->is_boundary_permutation[abs(d)]) {
            payload.place_recv_elements(buffer, d, par, vector_lattice);
        } else {
            // nothing to do here, comms directly in place
        }
    } else {
        // non-vectorized, using vanilla method, again nothing to do
    }
#else
    // this one is only for CUDA
    payload.place_comm_elements(d, par, buffer, from_node, lattice);
#endif
    // #endif
}

/////////////////////////////////////////////////////////////////////////////////////////
/// Place boundary elements from local lattice (used in vectorized version)

template <typename T>
void Field<T>::field_struct::set_local_boundary_elements(Direction dir, Parity par) {

#ifdef SPECIAL_BOUNDARY_CONDITIONS
    bool antiperiodic =
        (boundary_condition[dir] == hila::bc::ANTIPERIODIC && lattice.mynode.is_on_edge(dir));
#else
    bool antiperiodic = false;
#endif
    payload.set_local_boundary_elements(dir, par, lattice, antiperiodic);
}


/////////////////////////////////////////////////////////////////////////////////////////
/// get the receive buffer pointer for the communication.

template <typename T>
T *Field<T>::field_struct::get_receive_buffer(Direction d, Parity par,
                                              const lattice_struct::comm_node_struct &from_node) {
#if defined(VANILLA)

    return (T *)payload.get_buffer() + from_node.offset(par);

#elif defined(CUDA) || defined(HIP)

    unsigned offs = 0;
    if (par == ODD)
        offs = from_node.sites / 2;
    if (receive_buffer[d] == nullptr) {
        receive_buffer[d] = payload.allocate_mpi_buffer(from_node.sites);
    }
    return receive_buffer[d] + offs;

#elif defined(VECTORIZED)

    if constexpr (!hila::is_vectorizable_type<T>::value) {
        // use vanilla type, field laid out in std fashion
        return (T *)payload.get_buffer() + from_node.offset(par);
    } else {
        unsigned offs = 0;
        if (par == ODD)
            offs = from_node.sites / 2;

        if (vector_lattice->is_boundary_permutation[abs(d)]) {
            // extra copy operation needed
            if (receive_buffer[d] == nullptr) {
                receive_buffer[d] = payload.allocate_mpi_buffer(from_node.sites);
            }
            return receive_buffer[d] + offs;
        } else {
            // directly to halo buffer
            constexpr unsigned vector_size = hila::vector_info<T>::vector_size;
            return ((T *)payload.get_buffer() +
                    (vector_lattice->halo_offset[d] * vector_size + offs));
        }
    }
#endif
} // end of get_receive_buffer


#define NAIVE_SHIFT
#if defined(NAIVE_SHIFT)

template <typename T>
Field<T> &Field<T>::shift(const CoordinateVector &v, Field<T> &res, const Parity par) const {

    // use this to store remaining moves
    CoordinateVector rem = v;

    // check the parity of the move
    Parity par_s;

    int len = 0;
    foralldir(d) len += abs(rem[d]);

    // no move, just copy field
    if (len == 0) {
        res = *this;
        return res;
    }

    // opp_parity(ALL) == ALL
    if (len % 2 == 0)
        par_s = opp_parity(par);
    else
        par_s = par;

    // is this already gathered from one of the dirs in v?
    bool found_dir = false;
    Direction mdir;
    foralldir(d) {
        if (rem[d] > 0 && gather_status(par_s, d) != gather_status_t::NOT_DONE) {
            mdir = d;
            found_dir = true;
            break;
        } else if (rem[d] < 0 && gather_status(par_s, -d) != gather_status_t::NOT_DONE) {
            mdir = -d;
            found_dir = true;
            break;
        }
    }

    if (!found_dir) {
        // now did not find a 'ready' dir. Take the 1st available
        foralldir(d) {
            if (rem[d] > 0) {
                mdir = d;
                break;
            } else if (rem[d] < 0) {
                mdir = -d;
                break;
            }
        }
    }

    // Len 1, copy directly
    if (len == 1) {
        res[par_s] = (*this)[X + mdir];
        return res;
    }

    // now longer - need buffer
    Field<T> r1;
    Field<T> *from, *to;

    // this ensures that the final move lands on res
    if (len % 2 == 0) {
        from = &r1;
        to = &res;
    } else {
        from = &res;
        to = &r1;
    }
    // and copy initially to "from"
    (*from)[par_s] = (*this)[X + mdir];

    // and subtract remaining moves from rem
    rem = rem - mdir;
    par_s = opp_parity(par_s);

    foralldir(d) {
        if (rem[d] != 0) {
            mdir = (rem[d] > 0) ? d : -d;

            while (rem[d] != 0) {

                (*to)[par_s] = (*from)[X + mdir];

                par_s = opp_parity(par_s);
                rem = rem - mdir;
                std::swap(to, from);
            }
        }
    }

    return res;
}

#endif // NAIVE_SHIFT

/// start_gather(): Communicate the field at Parity par from Direction
/// d. Uses accessors to prevent dependency on the layout.
/// return the Direction mask bits where something is happening
template <typename T>
dir_mask_t Field<T>::start_gather(Direction d, Parity p) const {

    // get the mpi message tag right away, to ensure that we are always synchronized
    // with the mpi calls -- some nodes might not need comms, but the tags must be in
    // sync

    int tag = get_next_msg_tag();

    lattice_struct::nn_comminfo_struct &ci = lattice.nn_comminfo[d];
    lattice_struct::comm_node_struct &from_node = ci.from_node;
    lattice_struct::comm_node_struct &to_node = ci.to_node;

    // check if this is done - either gathered or no comm to be done in the 1st place

    if (is_gathered(d, p)) {
        lattice.n_gather_avoided++;
        return 0; // nothing to wait for
    }

    // No comms to do, nothing to wait for -- we'll use the is_gathered
    // status to keep track of vector boundary shuffle anyway

    if (from_node.rank == hila::myrank() && to_node.rank == hila::myrank()) {
        fs->set_local_boundary_elements(d, p);
        mark_gathered(d, p);
        return 0;
    }

    // if this parity or ALL-type gather is going on nothing to be done
    if (!gather_not_done(d, p) || !gather_not_done(d, ALL)) {
        lattice.n_gather_avoided++;
        return get_dir_mask(d); // nothing to do, but still need to wait
    }

    Parity par = p;
    // if p is ALL but ODD or EVEN is going on/done, turn off parity which is not needed
    // corresponding wait must do the same thing
    if (p == ALL) {
        if (!gather_not_done(d, EVEN) && !gather_not_done(d, ODD)) {
            // even and odd are going on or ready, nothing to be done
            lattice.n_gather_avoided++;
            return get_dir_mask(d);
        }
        if (!gather_not_done(d, EVEN))
            par = ODD;
        else if (!gather_not_done(d, ODD))
            par = EVEN;
        // if neither is the case par = ALL
    }

    mark_gather_started(d, par);

    // Communication hasn't been started yet, do it now

    int par_i = static_cast<int>(par) - 1; // index to dim-3 arrays

    constexpr size_t size = sizeof(T);

    T *receive_buffer;
    T *send_buffer;

    size_t size_type;
    MPI_Datatype mpi_type = get_MPI_number_type<T>(size_type);

    if (from_node.rank != hila::myrank() && boundary_need_to_communicate(d)) {

        // HANDLE RECEIVES: get node which will send here

        // buffer can be separate or in Field buffer
        receive_buffer = fs->get_receive_buffer(d, par, from_node);

        size_t n = from_node.n_sites(par) * size / size_type;

        if (n >= (1ULL << 31)) {
            hila::out << "Too large MPI message!  Size " << n << '\n';
            hila::terminate(1);
        }

        post_receive_timer.start();

        // c++ version does not return errors
        MPI_Irecv(receive_buffer, (int)n, mpi_type, from_node.rank, tag, lattice.mpi_comm_lat,
                  &fs->receive_request[par_i][d]);

        post_receive_timer.stop();
    }

    if (to_node.rank != hila::myrank() && boundary_need_to_communicate(-d)) {
        // HANDLE SENDS: Copy Field elements on the boundary to a send buffer and send

        unsigned sites = to_node.n_sites(par);

        if (fs->send_buffer[d] == nullptr)
            fs->send_buffer[d] = fs->payload.allocate_mpi_buffer(to_node.sites);

        send_buffer = fs->send_buffer[d] + to_node.offset(par);

        fs->gather_comm_elements(d, par, send_buffer, to_node);

        size_t n = sites * size / size_type;
#ifdef GPU_AWARE_MPI
        gpuStreamSynchronize(0);
        // gpuDeviceSynchronize();
#endif

        start_send_timer.start();

        MPI_Isend(send_buffer, (int)n, mpi_type, to_node.rank, tag, lattice.mpi_comm_lat,
                  &fs->send_request[par_i][d]);

        start_send_timer.stop();
    }

    // and do the boundary shuffle here, after MPI has started
    // NOTE: there should be no danger of MPI and shuffle overwriting, MPI writes
    // to halo buffers only if no permutation is needed.  With a permutation MPI
    // uses special receive buffer
    fs->set_local_boundary_elements(d, par);

    return get_dir_mask(d);
}

/// @internal
///  wait_gather(): Wait for communication at parity par from
///  Direction d completes the communication in the function.
///  If the communication has not started yet, also calls
///  start_gather()
///
///  NOTE: This will be called even if the field is marked const.
///  Therefore this function is const, even though it does change
///  the internal content of the field, the halo. From the point
///  of view of the user, the value of the field does not change.
template <typename T>
void Field<T>::wait_gather(Direction d, Parity p) const {

    lattice_struct::nn_comminfo_struct &ci = lattice.nn_comminfo[d];
    lattice_struct::comm_node_struct &from_node = ci.from_node;
    lattice_struct::comm_node_struct &to_node = ci.to_node;

    // check if this is done - either gathered or no comm to be done in the 1st place
    if (is_gathered(d, p))
        return;

    // this is the branch if no comms -- shuffle was done in start_gather
    if (from_node.rank == hila::myrank() && to_node.rank == hila::myrank())
        return;

    // if (!is_gather_started(d,p)) {
    //   hila::out0 << "Wait gather error - wait_gather without corresponding
    //   start_gather\n"; exit(1);
    // }

    // Note: the move can be Parity p OR ALL -- need to wait for it in any case
    // set par to be the "sum" over both parities
    // There never should be ongoing ALL and other parity gather -- start_gather takes
    // care

    // check here consistency, this should never happen
    if (p != ALL && is_gather_started(d, p) && is_gather_started(d, ALL)) {
        exit(1);
    }

    Parity par;
    int n_wait = 1;
    // what par to wait for?
    if (is_gather_started(d, p))
        par = p; // standard match
    else if (p != ALL) {
        if (is_gather_started(d, ALL))
            par = ALL; // if all is running wait for it
        else {
            exit(1);
        }
    } else {
        // now p == ALL and ALL is not running
        if (is_gathered(d, EVEN) && is_gather_started(d, ODD))
            par = ODD;
        else if (is_gathered(d, ODD) && is_gather_started(d, EVEN))
            par = EVEN;
        else if (is_gather_started(d, EVEN) && is_gather_started(d, ODD)) {
            n_wait = 2; // need to wait for both!
            par = ALL;
        } else {
            exit(1);
        }
    }

    if (n_wait == 2)
        par = EVEN; // we'll flip both

    for (int wait_i = 0; wait_i < n_wait; ++wait_i) {

        int par_i = (int)par - 1;

        if (from_node.rank != hila::myrank() && boundary_need_to_communicate(d)) {
            wait_receive_timer.start();

            MPI_Status status;
            MPI_Wait(&fs->receive_request[par_i][d], &status);

            wait_receive_timer.stop();

#ifndef VANILLA
            fs->place_comm_elements(d, par, fs->get_receive_buffer(d, par, from_node), from_node);
#endif
        }

        // then wait for the sends
        if (to_node.rank != hila::myrank() && boundary_need_to_communicate(-d)) {
            wait_send_timer.start();
            MPI_Status status;
            MPI_Wait(&fs->send_request[par_i][d], &status);
            wait_send_timer.stop();
        }

        // Mark the parity gathered from Direction dir
        mark_gathered(d, par);

        // Keep count of communications
        lattice.n_gather_done += 1;

        par = opp_parity(par); // flip if 2 loops
    }
}


/// Gather a list of elements to a single node
/// coord_list must be same on all nodes, buffer is needed only on "root"
template <typename T>
void Field<T>::field_struct::gather_elements(T *RESTRICT buffer,
                                             const std::vector<CoordinateVector> &coord_list,
                                             int root) const {

    std::vector<unsigned> index_list;
    std::vector<int> sites_on_rank(lattice.n_nodes());
    std::vector<int> reshuffle_list(coord_list.size());

    std::fill(sites_on_rank.begin(), sites_on_rank.end(), 0);

    int nranks = 0;

    int i = 0;
    for (const CoordinateVector &c : coord_list) {
        int rank = lattice.node_rank(c);
        if (hila::myrank() == rank) {
            index_list.push_back(lattice.site_index(c));
        }

        if (sites_on_rank[rank] == 0 && rank != root)
            nranks++;
        sites_on_rank[rank]++;
        reshuffle_list[i++] = rank;
    }

    std::vector<T> send_buffer(index_list.size());
    payload.gather_elements((T *)send_buffer.data(), index_list.data(), send_buffer.size(),
                            lattice);
    if (hila::myrank() != root && sites_on_rank[hila::myrank()] > 0) {
        MPI_Send((char *)send_buffer.data(), sites_on_rank[hila::myrank()] * sizeof(T), MPI_BYTE,
                 root, hila::myrank(), lattice.mpi_comm_lat);
    }
    if (hila::myrank() == root) {

        // allocate buffer for receiving data
        T *b;
        std::vector<T> pb(coord_list.size() - sites_on_rank[root]);
        b = pb.data();
        // vector for node ptrs -- point to stuff from nodes
        std::vector<T *> nptr(lattice.n_nodes());

        std::vector<MPI_Request> mpi_req(nranks);
        int nreqs = 0;
        for (int n = 0; n < sites_on_rank.size(); n++) {
            if (sites_on_rank[n] > 0) {
                if (n != root) {
                    MPI_Status status;
                    MPI_Irecv(b, (int)(sites_on_rank[n] * sizeof(T)), MPI_BYTE, n, n,
                              lattice.mpi_comm_lat, &mpi_req[nreqs++]);

                    nptr[n] = b;
                    b += sites_on_rank[n];

                } else {

                    nptr[n] = send_buffer.data();
                }
            }
        }

        if (nreqs > 0) {
            std::vector<MPI_Status> stat_arr(nreqs);
            MPI_Waitall(nreqs, mpi_req.data(), stat_arr.data());
        }

        // copy the data from bp to buffer, reordering
        for (int i = 0; i < coord_list.size(); i++) {
            buffer[i] = *nptr[reshuffle_list[i]];
            nptr[reshuffle_list[i]]++;
        }
    }
}

/// Send elements from a single node to a list of coordinates
/// coord_list must be the same on all nodes, but buffer is needed only on "root"!

template <typename T>
void Field<T>::field_struct::scatter_elements(T *RESTRICT buffer,
                                              const std::vector<CoordinateVector> &coord_list,
                                              int root) {

    std::vector<unsigned> index_list;
    std::vector<int> sites_on_rank(lattice.n_nodes());
    std::vector<int> reshuffle_list(coord_list.size());
    std::fill(sites_on_rank.begin(), sites_on_rank.end(), 0);

    int nranks = 0;
    int i = 0;
    for (CoordinateVector c : coord_list) {
        int rank = lattice.node_rank(c);
        if (hila::myrank() == rank) {
            index_list.push_back(lattice.site_index(c));
        }

        if (sites_on_rank[rank] == 0 && rank != root)
            nranks++;
        sites_on_rank[rank]++;
        reshuffle_list[i++] = rank;
    }

    // payload.gather_elements((T *)recv_buffer.data(), index_list.data(),
    //                         recv_buffer.size(), lattice);

    if (hila::myrank() != root && sites_on_rank[hila::myrank()] > 0) {
        std::vector<T> recv_buffer(index_list.size());
        MPI_Status status;

        MPI_Recv((char *)recv_buffer.data(), sites_on_rank[hila::myrank()] * sizeof(T), MPI_BYTE,
                 root, hila::myrank(), lattice.mpi_comm_lat, &status);

        payload.place_elements((T *)recv_buffer.data(), index_list.data(), recv_buffer.size(),
                               lattice);
    }
    if (hila::myrank() == root) {
        // reordering buffers
        std::vector<T> pb(coord_list.size());
        // vector for node counters -- point to stuff from nodes
        std::vector<unsigned> nloc(lattice.n_nodes());
        std::vector<unsigned> ncount(lattice.n_nodes());
        nloc[0] = ncount[0] = 0;

        for (int n = 1; n < lattice.n_nodes(); n++) {
            nloc[n] = nloc[n - 1] + sites_on_rank[n - 1];
            ncount[n] = 0;
        }
        for (int i = 0; i < coord_list.size(); i++) {
            int node = reshuffle_list[i];
            pb[nloc[node] + ncount[node]] = buffer[i];
            ncount[node]++;
        }

        std::vector<MPI_Request> mpi_req(nranks);
        int nreqs = 0;
        for (int n = 0; n < sites_on_rank.size(); n++) {
            if (sites_on_rank[n] > 0) {
                if (n != root) {
                    MPI_Isend(pb.data() + nloc[n], (int)(sites_on_rank[n] * sizeof(T)), MPI_BYTE, n,
                              n, lattice.mpi_comm_lat, &mpi_req[nreqs++]);
                }
            }
        }

        payload.place_elements(pb.data() + nloc[root], index_list.data(), index_list.size(),
                               lattice);

        if (nreqs > 0) {
            std::vector<MPI_Status> stat_arr(nreqs);
            MPI_Waitall(nreqs, mpi_req.data(), stat_arr.data());
        }
    }
}

template <typename T>
void Field<T>::set_elements(const std::vector<T> &elements,
                            const std::vector<CoordinateVector> &coord_list) {
    assert(elements.size() == coord_list.size() && "vector size mismatch in set_elments");
    std::vector<unsigned> my_indexes;
    std::vector<T> my_elements;
    for (int i = 0; i < coord_list.size(); i++) {
        CoordinateVector c = coord_list[i];
        if (lattice.is_on_mynode(c)) {
            my_indexes.push_back(lattice.site_index(c));
            my_elements.push_back(elements[i]);
        }
    }
    fs->payload.place_elements(my_elements.data(), my_indexes.data(), my_indexes.size(), lattice);
    mark_changed(ALL);
}


/// Get a list of elements and store them into a vector
template <typename T>
std::vector<T> Field<T>::get_elements(const std::vector<CoordinateVector> &coord_list,
                                      bool bcast) const {

    std::vector<T> res;
    if (hila::myrank() == 0)
        res.resize(coord_list.size());

    fs->gather_elements(res.data(), coord_list);
    if (bcast)
        hila::broadcast(res);

    return res;
}


/// get a subvolume of the field elements to all nodes
template <typename T>
std::vector<T> Field<T>::get_subvolume(const CoordinateVector &cmin, const CoordinateVector &cmax,
                                       bool bcast) const {

    size_t vol = 1;
    foralldir(d) {
        vol *= cmax[d] - cmin[d] + 1;
        assert(cmax[d] >= cmin[d] && cmin[d] >= 0 && cmax[d] < lattice.size(d));
    }
    std::vector<CoordinateVector> clist(vol);
    CoordinateVector c;

    size_t i = 0;
    forcoordinaterange(c, cmin, cmax) {
        clist[i++] = c;
    }
    return get_elements(clist, bcast);
}


/// and get a slice (subvolume)
template <typename T>
std::vector<T> Field<T>::get_slice(const CoordinateVector &c, bool bcast) const {
    CoordinateVector cmin, cmax;
    foralldir(d) if (c[d] < 0) {
        cmin[d] = 0;
        cmax[d] = lattice.size(d) - 1;
    }
    else {
        cmin[d] = cmax[d] = c[d];
    }
    return get_subvolume(cmin, cmax, bcast);
}


//////////////////////////////////////////////////////////////////////////////////
/// @internal
/// Copy the local (mpi process) data to a "logical array"
/// on gpu code, copies to host

template <typename T>
void Field<T>::copy_local_data(std::vector<T> &buffer) const {

    // copy to local variables to avoid lattice ptr
    CoordinateVector nmin = lattice.mynode.min;
    Vector<NDIM, unsigned> nmul = lattice.mynode.size_factor;

    buffer.resize(lattice.mynode.volume());
#if defined(CUDA) || defined(HIP)
    // d_malloc mallocs from device if needed
    T *data = (T *)d_malloc(sizeof(T) * lattice.mynode.volume());
#else
    T *data = buffer.data();
#endif

#pragma hila novector direct_access(data)
    onsites(ALL) {
        Vector<NDIM, unsigned> nodec;
        nodec = X.coordinates() - nmin;

        unsigned i = nodec.dot(nmul);
        data[i] = (*this)[X];
    }

#if defined(CUDA) || defined(HIP)
    gpuMemcpy(buffer.data(), data, sizeof(T) * lattice.mynode.volume(), gpuMemcpyDeviceToHost);
    d_free(data);
#endif
}

////////////////////////////////////////////////////////////////////////////////////
/// @internal
/// set the local data from an array

template <typename T>
void Field<T>::set_local_data(const std::vector<T> &buffer) {

    // copy to local variables to avoid lattice ptr
    CoordinateVector nmin = lattice.mynode.min;
    Vector<NDIM, unsigned> nmul = lattice.mynode.size_factor;

    assert(buffer.size() >= lattice.mynode.volume());

#if defined(CUDA) || defined(HIP)
    // d_malloc mallocs from device if needed
    T *data = (T *)d_malloc(sizeof(T) * lattice.mynode.volume());
    gpuMemcpy(data, buffer.data(), sizeof(T) * lattice.mynode.volume(), gpuMemcpyHostToDevice);
#else
    T *data = buffer.data();
#endif

#pragma hila novector direct_access(data)
    onsites(ALL) {
        Vector<NDIM, unsigned> nodec;
        nodec = X.coordinates() - nmin;

        unsigned i = nodec.dot(nmul);
        (*this)[X] = data[i];
    }

#if defined(CUDA) || defined(HIP)
    d_free(data);
#endif

    this->mark_changed(ALL);
}


////////////////////////////////////////////////////////////////////////////////////
/// Copy local data with halo - useful for visualization

template <typename T>
inline void collect_field_halo_data_(T *data, const Field<T> &src, Field<T> &dest,
                                     const Vector<NDIM, int> &dirs, int ndir) {

    // get the coords of the min point of the halo array
    CoordinateVector nmin = lattice.mynode.min;
    nmin.asArray() -= 1;

    // construct the mult vector to access the data
    Vector<NDIM, unsigned> nmul;
    nmul.e(0) = 1;
    for (int i = 1; i < NDIM; i++)
        nmul.e(i) = nmul.e(i - 1) * (lattice.mynode.size[i - 1] + 2);


    Vector<NDIM, int> node_min;
    Vector<NDIM, int> node_max;
    foralldir(d) {
        node_min[d] = lattice.mynode.min[d];
        node_max[d] = lattice.mynode.min[d] + lattice.mynode.size[d] - 1;
    }

#pragma hila novector direct_access(data)
    onsites(ALL) {
        Vector<NDIM, unsigned> nodec;
        CoordinateVector c = X.coordinates();
        bool gotit = true;

        for (int i = 0; i < ndir; i++) {
            Direction d = (Direction)dirs[i];
            if (c.e(d) == node_min[d]) {
                c.e(d) -= 1;
            } else if (c.e(d) == node_max[d]) {
                c.e(d) += 1;
            } else
                gotit = false;
        }

        Direction d = (Direction)dirs[ndir];
        if (gotit && c.e(d) == node_min[d]) {
            c.e(d) -= 1;
            nodec = c - nmin;
            data[nodec.dot(nmul)] = src[X - d];
            dest[X] = src[X - d];

        } else if (gotit && c.e(d) == node_max[d]) {
            c.e(d) += 1;
            nodec = c - nmin;
            data[nodec.dot(nmul)] = src[X + d];
            dest[X] = src[X + d];
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void Field<T>::copy_local_data_with_halo(std::vector<T> &buffer) const {

    // get the coords of the min point of the halo array
    CoordinateVector nmin = lattice.mynode.min;
    nmin.asArray() -= 1;

    // construct the mult vector to access the data
    Vector<NDIM, unsigned> nmul;
    nmul.e(0) = 1;
    for (int i = 1; i < NDIM; i++)
        nmul.e(i) = nmul.e(i - 1) * (lattice.mynode.size[i - 1] + 2);

    // full size of the buffer
    size_t siz = 1;
    foralldir(d) siz *= (lattice.mynode.size[d] + 2);

    buffer.resize(siz);
#if defined(CUDA) || defined(HIP)
    // d_malloc mallocs from device if needed
    T *data = (T *)d_malloc(sizeof(T) * siz);
#else
    T *data = buffer.data();
#endif

    // now collect bulk
#pragma hila novector direct_access(data)
    onsites(ALL) {
        Vector<NDIM, unsigned> nodec;
        nodec = X.coordinates() - nmin;

        unsigned i = nodec.dot(nmul);
        data[i] = (*this)[X];
    }

    // collect nn-halos

    Field<T> corners = 0;
    Field<T> corner2 = 0;
#if NDIM > 2
    Field<T> corner3 = 0;
#if NDIM > 3
    Field<T> corner4 = 0;
#endif
#endif

    Vector<NDIM, int> dirs;
    foralldir(d1) {
        dirs[0] = d1;
        // gather d1 halo
        collect_field_halo_data_(data, (*this), corners, dirs, 0);

        for (int d2 = d1 + 1; d2 < NDIM; ++d2) {
            dirs[1] = d2;
            collect_field_halo_data_(data, corners, corner2, dirs, 1);
#if NDIM > 2
            for (int d3 = d2 + 1; d3 < NDIM; ++d3) {
                dirs[2] = d3;
                collect_field_halo_data_(data, corner2, corner3, dirs, 2);
#if NDIM > 3
                for (int d4 = d3 + 1; d4 < NDIM; ++d4) {
                    dirs[3] = d4;
                    collect_field_halo_data_(data, corner3, corner4, dirs, 3);
                }
#endif
            }
#endif
        }
    }


#if defined(CUDA) || defined(HIP)
    gpuMemcpy(buffer.data(), data, sizeof(T) * siz, gpuMemcpyDeviceToHost);
    d_free(data);
#endif
}


#endif
