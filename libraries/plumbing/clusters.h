#ifndef HILA_CLUSTERS_H_
#define HILA_CLUSTERS_H_

#include "hila.h"

#include <algorithm>

#include "gpucub.h"

/**
 * @file clusters.h
 * @brief provides cluster finding tools
 * @details clusters are grown with the help of auxiliary Field variable of type
 * Field<uint8_t> or VectorField<uint8_t>. Elements of the first case accept values
 * 0 .. 254, allowing for 255 cluster "colours".
 *
 * Special value hila::Clusters::background indicates neutral, background sites.
 *
 * For the VectorField ("link") type non-zero value of the link[d][X] indicates that
 * points X and X+d are connected. In this case there are no classes.
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
 *     else cltype[X] = hila::Clusters::background;
 * }
 *
 * // create connected clusters of sites with cltype values 0 and 1
 * hila::Clusters cl(cltype);
 *
 * // above is equivalent to "hila::Clusters cl;  cl.find(cltype);"
 *
 * if (hila::myrank() == 0) {
 *     hila::out << "Got " << cl.number() << " clusters\n";
 *     for (int i = 0; i < cl.number(); i++) {
 *         hila::out << "Cluster " << i << " type " << cl[i].type() << " size " << cl[i].size() <<
 * '\n';
 *     }
 * }
 * ----------------------
 * Alternatively, using links:
 * ----------------------
 * VectorField<uint8_t> links{hila::Clusters::background};
 *
 * foralldir(d) {
 *     onsites(ALL) {
 *         if (<cond>) link[d][X] = 1;
 *     }
 * }
 * hila::Clusters::cl(links);
 * ...
 * // same functions as in first case
 * ----------------------
 *
 * Functions:
 *
 * Constructor: initialize with constructor hila::Clusters(const Field<uint8_t> & clustertype)
 * Example:
 *      hila::Clusters cluster(cltype);   // initialize var cluster
 *
 * where cltype is of type Field<uint8_t>.  This constructor builds connected clusters of the
 * sites with the same cltype value. Special value of hila::Clusters::background are ignored.
 *
 * Above is equivalent with
 *      hila::Clusters cluster;
 *      cluster.find(cltype);
 *
 * Note: initalization is an expensive operation
 *
 * size_t hila::Clusters::number() - return the total number of clusters.
 * Example:
 *      hila::out0 << "Found " << cluster.number() << " clusters\n";
 *
 * The properties of the cluster number i are accessed with array index syntax []:
 * hila::Clusters::cluster_ref  hila::Clusters::operator[](int64_t)
 * i.e. `cluster[i]` returns a const reference to cluster number i
 *
 * Properties can be queried by functions
 * int64_t hila::Clusters::cluster_ref::size() - the size of cluster
 * int64_t hila::Clusters::cluster_ref::area() - the surface area of cluster
 * std::vector<SiteIndex> hila::Clusters::cluster_ref::sites() - vector of cluster sites
 * uint8_t hila::Clusters::cluster_ref::type() - the type of the cluster, 0..254
 *
 * area() and sites() are expensive operators and must be called by all MPI ranks
 * Area is defined by the number of links where one end belongs to the cluster, another does not.
 *
 * Example:
 *      Field<uint8_t> cltypes:
 *      ... // fill cltypes
 *
 *      hila::Clusters cl(cltypes);   // create clusters
 *
 *      hila::out0 << "Number of clusters found: " << cl.number() << '\n';
 *      for (int64_t i = 0; i < cl.number(); ++i) {
 *           hila::out0 << "Cluster " << i << " size: " << cl[i].size()
 *                      << " surface area: " << cl[i].area()
 *                      << " type " << cl[i].type() << '\n';
 *      }
 *
 *      hila::out0 << "Cluster 0 coordinates:\n";
 *      auto sitelist = cl[0].sites();
 *      for (auto s : sitelist) {
 *           hila::out0 << s.coordinates() << '\n';
 *      }
 *
 *
 */


namespace hila {

#define CLUSTER_BACKGROUND_ 0xFF

template <typename inttype, std::enable_if_t<std::is_integral<inttype>::value, int> = 0>
inline uint64_t set_cl_label(const CoordinateVector &cv, const inttype type) {
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


class Clusters {

  private:
    // We'll use uint64_t to encode the site and the type of the site.
    // Set the top 8 bits to type, and the bottom 54 bits to SiteIndex.
    Field<uint64_t> labels;

    struct cl_struct {
        uint64_t label;
        int64_t size;
    };

    // clist vector contains information for all clusters. NOTE: the
    // content is valid only on rank == 0.
    std::vector<cl_struct> clist;

    class cluster_ref {
        const Clusters &clusters;
        size_t index;

      public:
        cluster_ref() = delete;
        cluster_ref(const Clusters &clust, const size_t i) : clusters(clust), index(i) {}
        ~cluster_ref() = default;

        /// @brief return size of cluster
        auto size() const {
            return clusters.clist[index].size;
        }

        /// @brief return type of the cluster
        uint8_t type() const {
            return get_cl_label_type(clusters.clist[index].label);
        }

        /// @brief return the label of cluster i
        /// Each cluster has unique 64-bit label
        uint64_t label() const {
            return clusters.clist[index].label;
        }

        /// @brief returns the area of the cluster
        /// Must be called from all MPI ranks - involves reduction
        int64_t area() const {

            uint64_t lbl = clusters.clist[index].label;
            int64_t cl_area = 0;

            onsites (ALL) {
                if (clusters.labels[X] == lbl) {
                    for (Direction d = e_x; d < NDIRS; ++d) {
                        if (clusters.labels[X + d] != lbl)
                            cl_area += 1;
                    }
                }
            }

            return cl_area;
        }

        /// @brief returns std::vector of SiteIndex for cluster
        /// Expensive operation, can possibly overflow the memory of the node 0
        /// Must be called from all MPI ranks

        std::vector<SiteIndex> sites() const {
            SiteSelect sites;
            uint64_t la = clusters.clist[index].label;

            onsites (ALL) {
                if (clusters.labels[X] == la)
                    sites.select(X);
            }
            // hila::out << "Node " << hila::myrank() << " got " << sites.size() << " sites\n";
            return sites.move_sites();
        }

        ///////////////////////////////////////////////////////////////////
        /// Make this iterator by the rest of the definitions
        /// Now we can use 
        /// for (auto r : cls) hila::out0 << r.size();

        using iterator_category = std::input_iterator_tag;
        using difference_type = std::ptrdiff_t;

        /// deref operator returns the item, not a pointer to it!
        /// this tomfoolery makes the iterator work
        auto operator*() const {
            return *this;
        }

        auto operator++() {
            index++;
            return *this;
        }
        auto operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }
        friend bool operator!=(const hila::Clusters::cluster_ref &a,
                               const hila::Clusters::cluster_ref &b) {
            return a.index != b.index;
        }
        friend bool operator==(const hila::Clusters::cluster_ref &a,
                               const hila::Clusters::cluster_ref &b) {
            return a.index == b.index;
        }

        /// iterator definitions end
        /////////////////////////////////////////////////////////////////

    };


    inline void make_local_clist();

    void assert_cl_index(size_t i) const {
        if (i >= clist.size()) {
            hila::out0 << "Incorrect cluster index " << i << ", there are only " << clist.size()
                       << " clusters\n";
            exit(0);
        }
    }


  public:
    /// @brief hila::Clusters::background is special value to mark "non-interesting" sites.
    /// Using this can make finding clusters faster
    static constexpr uint8_t background = CLUSTER_BACKGROUND_; // 8 ones here = 255

    Clusters() = default;
    ~Clusters() = default;

    template <typename inttype, std::enable_if_t<std::is_integral<inttype>::value, int> = 0>
    Clusters(const Field<inttype> &type) {
        find(type);
    }

    template <typename inttype, std::enable_if_t<std::is_integral<inttype>::value, int> = 0>
    Clusters(const VectorField<inttype> &vfield) {
        find(vfield);
    }

    /// @brief find nearest-neighbour -connected clusters which have the same type
    /// @param type field which contains the type of the site, possible values 0-254.
    /// special value hila::Clusters::background indicates site does not belong to any cluster
    template <typename inttype, std::enable_if_t<std::is_integral<inttype>::value, int> = 0>
    void find(const Field<inttype> &type) {
        make_labels(type);
        classify();
    }

    /// @brief nearest-neighbour -connected clusters which connecting link
    /// @param link[d][X] is non-zero if sites X and X + d are to be connected
    template <typename inttype, std::enable_if_t<std::is_integral<inttype>::value, int> = 0>
    void find(const VectorField<inttype> &link) {
        Field<uint8_t> f{hila::Clusters::background};
        foralldir (d) {
            onsites (ALL) {
                if (link[d][X] != 0 || link[d][X - d] != 0)
                    f[X] = 1;
            }
        }
        find(f);
    }


    /// @brief number of clusters found
    size_t number() const {
        return clist.size();
    }

    /// @brief number of clusters found
    size_t size() const {
        return clist.size();
    }

    /// @brief access cluster number i
    const auto operator[](size_t i) const {
        assert_cl_index(i);
        return cluster_ref(*this, i);
    }


    cluster_ref begin() {
        return cluster_ref(*this, 0);
    }
    cluster_ref end() {
        return cluster_ref(*this, clist.size());
    }


    /// @brief make only the Field representation of cluster labels, without counting the clusters
    /// @param type - input field classifying the sites
    /// @return a const reference to label Field
    /// @details Usually this call is not needed, use hila::Clusters::find() or constructor
    template <typename inttype, std::enable_if_t<std::is_integral<inttype>::value, int> = 0>
    void make_labels(const Field<inttype> &type) {

        constexpr unsigned comm_interval = 1;

        // mark every site with site index on 54 low, leaving 8 bits at the top for the type
        labels[ALL] = set_cl_label(X.coordinates(), type[X]);

        Reduction<int64_t> changed;
        changed.delayed();
        unsigned i = 0;
        do {
            if (i % comm_interval != 0) {
                for (Parity par : {EVEN, ODD}) {
                    for (Direction d = e_x; d < NDIRS; ++d)
                        labels.mark_gathered(d, ALL);
                    // take away annoying warning with pragma
                    #pragma hila safe_access(labels)
                    onsites (par) {
                        auto type_0 = get_cl_label_type(labels[X]);
                        if (type_0 != hila::Clusters::background) {
                            for (Direction d = e_x; d < NDIRS; ++d) {
                                auto label_1 = labels[X + d];
                                auto type_1 = get_cl_label_type(label_1);
                                if (type_0 == type_1 && labels[X] > label_1) {
                                    labels[X] = label_1;
                                }
                            }
                        }
                    }
                }
            } else {
                changed = 0;
                for (Parity par : {EVEN, ODD}) {
                    // take away annoying warning with pragma
                    #pragma hila safe_access(labels)
                    onsites (par) {
                        auto type_0 = get_cl_label_type(labels[X]);
                        if (type_0 != hila::Clusters::background) {
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
            }
            i++;
            // hila::out0 << "g " << changed.value() << '\n';
        } while ((i - 1) % comm_interval != 0 || changed.value() > 0);

        // hila::out0 << "cluster labels: " << i << " iterations\n";
    }

    /// @brief obtain const refence to cluster label Field var
    /// Clusters::find() must have been called before this

    const auto &get_labels() const {
        return labels;
    }

    /// @brief classify clusters from existing label field (make_labels() called before)
    /// @details Constructs vector with cluster info, which can be queried
    /// Normally not needed, hila::Clusters::find() already does the search
    void classify() {

        if (!labels.is_initialized(ALL)) {
            hila::error("hila::Clusters::classify() called without preceding "
                        "hila::Clusters::make_labels()");
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

inline void Clusters::make_local_clist() {

    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    uint64_t *buf;

    gpuMalloc(&buf, lattice->mynode.volume * sizeof(uint64_t));
    GPU_CHECK(gpucub::DeviceRadixSort::SortKeys(
        d_temp_storage, temp_storage_bytes, labels.field_buffer(), buf, lattice->mynode.volume));

    // Allocate temporary storage
    gpuMalloc(&d_temp_storage, temp_storage_bytes);

    GPU_CHECK(gpucub::DeviceRadixSort::SortKeys(
        d_temp_storage, temp_storage_bytes, labels.field_buffer(), buf, lattice->mynode.volume));

    gpuFree(d_temp_storage);
    // now buf contains sorted labels

    d_temp_storage = nullptr;
    temp_storage_bytes = 0;
    uint64_t *d_labels;
    int64_t *d_sizes, *d_num_clust;
    gpuMalloc(&d_labels, lattice->mynode.volume * sizeof(uint64_t));
    gpuMalloc(&d_sizes, (lattice->mynode.volume + 1) * sizeof(int64_t));
    d_num_clust = d_sizes + lattice->mynode.volume; // the last element

    GPU_CHECK(gpucub::DeviceRunLengthEncode::Encode(d_temp_storage, temp_storage_bytes, buf,
                                                    d_labels, d_sizes, d_num_clust,
                                                    lattice->mynode.volume));

    // Allocate temporary storage
    gpuMalloc(&d_temp_storage, temp_storage_bytes);

    // Run encoding
    GPU_CHECK(gpucub::DeviceRunLengthEncode::Encode(d_temp_storage, temp_storage_bytes, buf,
                                                    d_labels, d_sizes, d_num_clust,
                                                    lattice->mynode.volume));

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
        clist[i] = {labels[i], sizes[i]};
    }
    if (clist.size() > 0 && get_cl_label_type(clist.back().label) == hila::Clusters::background)
        clist.pop_back();
}

#else

inline void Clusters::make_local_clist() {

    // this changes the ordering of labels - copy it
    auto lb = labels;
    lb.will_change();
    auto *buf = lb.field_buffer();
    std::sort(buf, buf + lattice->mynode.volume);
    // background labels are last after sort, stop handling when these are met
    int64_t i = 0;
    while (i < lattice->mynode.volume && get_cl_label_type(buf[i]) != hila::Clusters::background) {
        int64_t istart = i;
        auto label = buf[i];
        for (++i; i < lattice->mynode.volume && buf[i] == label; ++i)
            ;
        clist.push_back({label, i - istart});
    }
}

#endif


} // namespace hila


#endif
