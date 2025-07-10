#ifndef CLUSTERS_H_
#define CLUSTERS_H_

#include "hila.h"

#include <algorithm>

#include "gpucub.h"

/**
 * @file clusters.h
 * @brief provides cluster finding tools
 * @details clusters are grown with the help of auxiliary Field variable of type
 * Field<uint8_t>. Elements of this accept values 0 .. 254, allowing for 255 cluster classes.
 * Special value hila::clusters::background indicates neutral, background sites.
 *
 * Example:
 * ----------------------
 * #include "clusters.h"
 *
 * Field<uint8_t> cltype;
 * // initialize field cltype to values 0, 1 or background
 * onsites(ALL) {
 *     if (<some condition 1>) cltype[X] = 0;
 *     else if (condition 2>) cltype[X] = 1;
 *     else cltype[X] = hila::clusters::background;
 * }
 *
 * // create connected clusters of sites with cltype values 0 and 1
 * hila::clusters cl(cltype);
 *
 * // above is equivalent to "hila::clusters cl;  cl.find(cltype);"
 *
 * if (hila::myrank() == 0) {
 *     hila::out << "Got " << cl.number() << " clusters\n";
 *     for (int i = 0; i < cl.number(); i++) {
 *         hila::out << "Cluster " << i << " type " << cl.type(i) << " size " << cl.size(i) << '\n';
 *     }
 * }
 * ----------------------
 *
 * Functions:
 *
 * Constructor: initialize with constructor hila::clusters(const Field<uint8_t> & clustertype)
 * Example:
 *      hila::clusters cluster(cltype);   // initialize var cluster
 *
 * where cltype is of type Field<uint8_t>.  This constructor builds connected clusters of the
 * sites with the same cltype value. Special value of hila::clusters::background are ignored.
 *
 * Above is equivalent with
 *      hila::clusters cluster;
 *      cluster.find(cltype);
 *
 * Note: initalization is a relatively expensive operation
 *
 * size_t hila::clusters::number() - return the total number of clusters.
 * Example:
 *      hila::out0 << "Found " << cluster.number() << " clusters\n";
 *
 * int64_t hila::clusters::size(size_t cl_number) - return the size of cluster number cl_number
 * Example:
 *      hila::out0 << "Size of cluster 0 is " << cluster.size(0) << '\n';
 *
 * uint8_t hila::clusters::type(size_t cl_number) - return the cluster type
 *      hila::out0 << "Type of cluster 0 is " << cluster.size(0) << '\n';
 *
 * std::vector<SiteIndex> hila::clusters::sites(size_t cl_number) - return the sites of cluster
 * cl_number Example: print the coordinates of cluster number 0: auto clsites = cluster.sites(0);
 *      for (auto & r : clsites) {
 *          hila::out0 << r.coordinates() << '\n;
 *      }
 * Note: .sites() must be called by all MPI ranks, otherwise deadlock occurs
 *
 * int64_t hila::clusters::area(size_t cl_number) - return the area of the cluster
 * Area is defined by the number of links where one end belongs to the cluster, another does not.
 * Note: .area() must be called by all MPI ranks
 */


namespace hila {

#define CLUSTER_BACKGROUND_ 0xFF

inline uint64_t set_cl_label(const CoordinateVector &cv, const uint8_t type) {
    uint64_t r;
    // we'll give the same label for all background sites
    if (type == CLUSTER_BACKGROUND_)
        r = static_cast<uint64_t>(CLUSTER_BACKGROUND_) << 54;
    else
        r = SiteIndex(cv).value | (static_cast<uint64_t>(type) << 54);
    return r;
}

// return the top 8 bits of the encoded site+type
inline uint8_t get_cl_label_type(const uint64_t v) {
    return 0xff & (v >> 54);
}


class clusters {

  private:
    struct cl_struct {
        uint64_t label;
        int64_t size, area;
    };

    // clist vector contains information for all clusters. NOTE: the
    // content is valid only on rank == 0.
    std::vector<cl_struct> clist;

    // We'll use uint64_t to encode the site and the type of the site.
    // Set the top 8 bits to type, and the bottom 54 bits to SiteIndex.
    Field<uint64_t> labels;

    inline void make_local_clist();

    void assert_cl_index(size_t i) const {
        if (i >= clist.size()) {
            hila::out0 << "Too large cluster index " << i << ", there are only " << clist.size()
                       << " clusters\n";
            exit(0);
        }
    }


  public:
    /// hila::clusters::background is special value to mark "non-interesting" sites.
    /// Using this can make finding clusters faster
    static constexpr uint8_t background = CLUSTER_BACKGROUND_; // 8 ones here = 255

    clusters() = default;
    ~clusters() = default;
    clusters(const Field<uint8_t> &type) {
        find(type);
    }

    /// @brief find nearest-neighbour -connected clusters which have the same type
    /// @param type field which contains the type of the site, possible values 0-254.
    /// special value hila::clusters::background indicates site does not belong to any cluster
    void find(const Field<uint8_t> &type) {
        make_labels(type);
        classify();
    }

    /// @brief number of clusters found

    size_t number() const {
        return clist.size();
    }

    /// @brief return size of cluster i
    /// @param i cluster number

    int64_t size(size_t i) const {
        assert_cl_index(i);
        return clist[i].size;
    }

    /// @brief return type of the cluster
    /// @param i cluster number

    uint8_t type(size_t i) const {
        assert_cl_index(i);
        return get_cl_label_type(clist[i].label);
    }

    /// @brief return the label of cluster i
    /// @param i cluster number
    /// Each cluster has unique 64-bit label

    uint64_t label(size_t i) const {
        assert_cl_index(i);
        return clist[i].label;
    }


    /// @brief returns the area of the cluster number i
    /// @param i cluster number
    /// First time call involves reduction, after that value buffered
    /// Must be called from all MPI ranks

    int64_t area(size_t i) {
        assert_cl_index(i);
        // check if this is already computed
        if (clist[i].area > 0)
            return clist[i].area;

        uint64_t lbl = label(i);
        uint64_t cl_area = 0;

        onsites(ALL) {
            if (labels[X] == lbl) {
                for (Direction d = e_x; d < NDIRS; ++d) {
                    if (labels[X + d] != lbl)
                        cl_area += 1;
                }
            }
        }

        clist[i].area = cl_area;
        return cl_area;
    }

    /// @brief returns std::vector of SiteIndex for cluster number i
    /// @param i cluster number
    /// Expensive operation, can possibly overflow the memory of the node 0
    /// Must be called from all MPI ranks

    std::vector<SiteIndex> sites(size_t i) const {
        SiteSelect sites;
        uint64_t la = label(i);

        onsites(ALL) {
            if (labels[X] == la)
                sites.select(X);
        }
        hila::out << "Node " << hila::myrank() << " got " << sites.size() << " sites\n";
        return sites.move_sites();
    }


    /// @brief make only the Field representation of cluster labels, without counting the clusters
    /// @param type - input field classifying the sites
    /// @return a const reference to label Field
    /// @details Usually this call is not needed, use hila::clusters::find() or constructor
    void make_labels(const Field<uint8_t> &type) {

        // mark every site with site index on 54 low, leaving 8 bits at the top for the type
        labels[ALL] = set_cl_label(X.coordinates(), type[X]);

        Reduction<int64_t> changed;
        changed.delayed();
        do {
            changed = 0;
            for (Parity par : {EVEN, ODD}) {
                // take away annoying warning
#pragma hila safe_access(labels)
                onsites(par) {
                    auto type_0 = get_cl_label_type(labels[X]);
                    if (type_0 != background) {
                        for (Direction d = e_x; d < NDIRS; ++d) {
                            auto label_1 = labels[X + d];
                            auto type_1 = get_cl_label_type(label_1);
                            if (type_0 == type_1 && labels[X] > label_1) {
                                labels[X] = label_1;
                                changed += 1;
                            }
                        }
                    }
                }
            }

            // hila::out0 << "g " << changed.value() << '\n';

        } while (changed.value() > 0);
    }

    /// @brief obtain const refence to cluster label Field var
    /// clusters::find() must have been called before this

    const auto &get_labels() const {
        return labels;
    }

    /// @brief classify clusters from existing label field (make_labels() called before)
    /// @details Constructs vector with cluster info, which can be queried
    /// Normally not needed, hila::clusters::find() already does the search
    void classify() {

        if (!labels.is_initialized(ALL)) {
            hila::error("hila::clusters::classify() called without preceding "
                        "hila::clusters::make_labels()");
        }

        clist.clear();

        make_local_clist();

        // Now merge the clist across mpi ranks
        // communicate and merge in node pairs

        int nn = hila::number_of_nodes();
        int myrank = hila::myrank();

        for (int step = 1; step < nn; step *= 2) {
            if (myrank % (2 * step) == step) {

                // we're "odd" pair, send to "even"
                hila::send_to(myrank - step, clist);
                clist.clear(); // in order to save space

            } else if (myrank % (2 * step) == 0 && myrank + step < nn) {

                // Now "even", receive from "odd" and merge data
                std::vector<cl_struct> clist_n, merged;
                hila::receive_from(myrank + step, clist_n);

                int i = 0, n = 0;
                while (i < clist.size() && n < clist_n.size()) {
                    if (clist[i].label < clist_n[n].label) {
                        merged.push_back(clist[i++]);

                    } else if (clist[i].label == clist_n[n].label) {
                        merged.push_back(clist[i++]);
                        merged.back().size += clist_n[n++].size;
                    } else {
                        merged.push_back(clist_n[n++]);
                    }
                }
                // only 1 of the following is done (if either)
                while (i < clist.size()) {
                    merged.push_back(clist[i++]);
                }
                while (n < clist_n.size()) {
                    merged.push_back(clist_n[n++]);
                }
                clist_n.clear();
                clist = std::move(merged);
            }
        }

        hila::barrier();

        // deliver clist to all ranks after all
        hila::broadcast(clist);

    } // void classify()


}; // class clusters


#if (defined(CUDA) || defined(HIP)) && !defined(HILAPP)

inline void clusters::make_local_clist() {

    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    uint64_t *buf;

    gpuMalloc(&buf, lattice.mynode.volume() * sizeof(uint64_t));
    GPU_CHECK(gpucub::DeviceRadixSort::SortKeys(
        d_temp_storage, temp_storage_bytes, labels.field_buffer(), buf, lattice.mynode.volume()));

    // Allocate temporary storage
    gpuMalloc(&d_temp_storage, temp_storage_bytes);

    GPU_CHECK(gpucub::DeviceRadixSort::SortKeys(
        d_temp_storage, temp_storage_bytes, labels.field_buffer(), buf, lattice.mynode.volume()));

    gpuFree(d_temp_storage);
    // now buf contains sorted labels

    d_temp_storage = nullptr;
    temp_storage_bytes = 0;
    uint64_t *d_labels;
    int64_t *d_sizes, *d_num_clust;
    gpuMalloc(&d_labels, lattice.mynode.volume() * sizeof(uint64_t));
    gpuMalloc(&d_sizes, (lattice.mynode.volume() + 1) * sizeof(int64_t));
    d_num_clust = d_sizes + lattice.mynode.volume(); // the last element

    GPU_CHECK(gpucub::DeviceRunLengthEncode::Encode(d_temp_storage, temp_storage_bytes, buf,
                                                    d_labels, d_sizes, d_num_clust,
                                                    lattice.mynode.volume()));

    // Allocate temporary storage
    gpuMalloc(&d_temp_storage, temp_storage_bytes);

    // Run encoding
    GPU_CHECK(gpucub::DeviceRunLengthEncode::Encode(d_temp_storage, temp_storage_bytes, buf,
                                                    d_labels, d_sizes, d_num_clust,
                                                    lattice.mynode.volume()));

    gpuFree(d_temp_storage);

    int64_t nclusters;
    gpuMemcpy(&nclusters, d_num_clust, sizeof(int64_t), gpuMemcpyDeviceToHost);

    std::vector<uint64_t> labels(nclusters);
    std::vector<int64_t> sizes(nclusters);
    if (nclusters > 0) {
        gpuMemcpy(labels.data(), d_labels, sizeof(uint64_t) * nclusters, gpuMemcpyDeviceToHost);
        gpuMemcpy(sizes.data(), d_sizes, sizeof(int64_t) * nclusters, gpuMemcpyDeviceToHost);
    }

    gpuFree(d_sizes);
    gpuFree(d_labels);

    clist.resize(nclusters);
    for (int i = 0; i < nclusters; i++) {
        clist[i] = {labels[i], sizes[i], 0};
    }
    if (clist.size() > 0 && get_cl_label_type(clist.back().label) == hila::clusters::background)
        clist.pop_back();
}

#else

inline void clusters::make_local_clist() {

    // this changes the ordering of labels - copy it
    auto lb = labels;
    auto *buf = lb.field_buffer();
    std::sort(buf, buf + lattice.mynode.volume());
    // background labels are last after sort, stop handling when these are met
    int64_t i = 0;
    while (i < lattice.mynode.volume() && get_cl_label_type(buf[i]) != background) {
        int64_t istart = i;
        auto label = buf[i];
        for (++i; i < lattice.mynode.volume() && buf[i] == label; ++i)
            ;
        clist.push_back({label, i - istart, 0});
    }
}

#endif


} // namespace hila


#endif
