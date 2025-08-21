#ifndef HILA_SITE_SELECT_H_
#define HILA_SITE_SELECT_H_

// We insert the GPU code in the same file too
// hilapp should not read in .cuh, because it does not understand it


#include "hila.h"

#include "gpucub.h"

//////////////////////////////////////////////////////////////////////////////////
/// Site selection: special vector to accumulate chosen sites or sites + variable
///
/// SiteSelect s;
/// SiteValueSelect<T> sv;
///
/// To be used within site loops as
///   onsites(ALL ) {
///       if ( condition1 )
///           s.select(X);
///       if ( condition2 )
///           sv.select(X, A[X]);
///   }
///
///
///

// just an empty class used to flag select operations
class site_select_type_ {};

class SiteSelect {

  protected:
    std::vector<SiteIndex> sites;


    /// status variables of reduction
    bool auto_join = true;
    bool joined = false;

    // max number of elements to collect - default volume
    size_t nmax = lattice.volume();

    size_t current_index = 0;
    size_t previous_site = SIZE_MAX;
    size_t n_overflow = 0;

  public:
    /// Initialize to zero by default (? exception to other variables)
    /// allreduce = true by default
    explicit SiteSelect() {
        auto_join = true;
        joined = false;
        nmax = lattice.volume();
        current_index = 0;
        previous_site = SIZE_MAX;
        n_overflow = 0;
    }

    SiteSelect(const SiteSelect &a) = default;

    /// Destructor cleans up communications if they are in progress
    ~SiteSelect() = default;

    /// Selection - use only inside loops

    site_select_type_ select(const X_index_type x) {
        return site_select_type_();
        // filled in by hilapp
    }

    // this makes sense only for cpu targets
    void select_site(const SiteIndex s) {
        if (s.value == previous_site) {
            sites[current_index - 1] = s;
        } else {
            sites[current_index] = s;
            previous_site = s.value;
            current_index++;
        }
    }

    SiteSelect &no_join() {
        auto_join = false;
        return *this;
    }

    SiteSelect &max_size(size_t _max) {
        nmax = _max;
        return *this;
    }

    void setup() {
        sites.resize(lattice->mynode.volume);
        current_index = 0;
        previous_site = SIZE_MAX;
        n_overflow = 0;
        joined = false;
    }

    void clear() {
        sites.clear();
        current_index = 0;
        previous_site = SIZE_MAX;
    }

    size_t size() const {
        return sites.size();
    }

    const CoordinateVector coordinates(size_t i) const {
        return sites.at(i).coordinates();
    }

    const SiteIndex site_index(size_t i) const {
        return sites.at(i);
    }

    // Don't even implement assignments

    /// @brief  std::move SiteIndex vector of selected sites, invalidating this variable
    std::vector<SiteIndex> move_sites() {
        return std::move(sites);
    }

    void join() {
        if (!joined) {
            std::vector<std::nullptr_t> v;
            join_data_vectors(v);
            joined = true;
        }
    }

    /// For delayed collect, joining starts or completes the reduction operation
    template <typename T>
    void join_data_vectors(std::vector<T> &dp) {
        if (hila::myrank() == 0) {
            for (int n = 1; n < hila::number_of_nodes(); n++) {
                size_t nsend = nmax - sites.size();
                hila::send_to(n, nsend);

                if (nsend > 0) {
                    std::vector<SiteIndex> s;
                    hila::receive_from(n, s);

                    // last element of s contains the overflow number
                    n_overflow += s.back().value;
                    s.pop_back();

                    sites.reserve(sites.size() + s.size());
                    sites.insert(sites.end(), s.begin(), s.end());

                    if constexpr (!std::is_same<T, std::nullptr_t>::value) {
                        std::vector<T> recvdata;
                        hila::receive_from(n, recvdata);
                        dp.reserve(sites.size());
                        dp.insert(dp.end(), recvdata.begin(), recvdata.end());
                    }
                } else {
                    // get the overflow number in any case
                    size_t over;
                    hila::receive_from(n, over);
                    n_overflow += over;
                }
            }

        } else {
            // now rank /= 0
            // wait for the number to be sent
            size_t nsend;
            hila::receive_from(0, nsend);
            if (nsend > 0) {
                if (nsend < sites.size()) {
                    n_overflow += sites.size() - nsend;
                    sites.resize(nsend);
                }

                // append overflow info
                sites.push_back(n_overflow);
                hila::send_to(0, sites);

                if constexpr (!std::is_same<T, std::nullptr_t>::value) {
                    dp.resize(sites.size() - 1);
                    hila::send_to(0, dp);
                }

            } else {
                // send overflow
                hila::send_to(0, sites.size() + n_overflow);
            }
            // empty data to release space
            clear();
        }
    }

    size_t overflow() {
        return n_overflow;
    }

#if !(defined(CUDA) || defined(HIP)) || defined(HILAPP)

    void endloop_action() {
        if (current_index > nmax) {
            // too many elements, trunc
            n_overflow = current_index - nmax;
            current_index = nmax;
        }
        sites.resize(current_index);
        if (auto_join)
            join();
    }

#else

    // this is GPU version of endloop_action
    // skip this for hilapp
    template <typename T>
    void copy_data_to_host_vector(std::vector<T> &dvec, const char *flag, const T *d_data) {
        void *d_temp_storage = nullptr;
        size_t temp_storage_bytes = 0;

        T *out;
        gpuMalloc(&out, lattice->mynode.volume * sizeof(T));

        int *num_selected_d;
        gpuMalloc(&num_selected_d, sizeof(int));


        GPU_CHECK(gpucub::DeviceSelect::Flagged(d_temp_storage, temp_storage_bytes, d_data, flag,
                                                out, num_selected_d, lattice->mynode.volume));

        gpuMalloc(&d_temp_storage, temp_storage_bytes);

        GPU_CHECK(gpucub::DeviceSelect::Flagged(d_temp_storage, temp_storage_bytes, d_data, flag,
                                                out, num_selected_d, lattice->mynode.volume));

        gpuFree(d_temp_storage);

        int num_selected;
        gpuMemcpy(&num_selected, num_selected_d, sizeof(int), gpuMemcpyDeviceToHost);
        gpuFree(num_selected_d);

        if (num_selected > nmax) {
            n_overflow = num_selected - nmax;
            num_selected = nmax;
        }
        dvec.resize(num_selected);

        gpuMemcpy(dvec.data(), out, sizeof(T) * num_selected, gpuMemcpyDeviceToHost);
        gpuFree(out);
    }

    // endloop action for this
    void endloop_action(const char *flag, const SiteIndex *d_sites) {

        copy_data_to_host_vector(sites, flag, d_sites);

        if (auto_join)
            join();
    }

#endif // GPU
};

class site_value_select_type_ {};

template <typename T>
class SiteValueSelect : public SiteSelect {
  protected:
    std::vector<T> values;

  public:
    explicit SiteValueSelect() : SiteSelect() {
        values.clear();
    }
    ~SiteValueSelect() = default;
    SiteValueSelect(const SiteValueSelect &v) = default;

    void setup() {
        SiteSelect::setup();
        values.resize(lattice->mynode.volume);
    }

    void clear() {
        SiteSelect::clear();
        values.clear();
    }

    site_value_select_type_ select(const X_index_type x, const T &val) {
        return site_value_select_type_();
    }

    void select_site_value(const SiteIndex s, const T &val) {
        values[current_index] = val;
        SiteSelect::select_site(s);
    }


    T value(size_t i) {
        return values.at(i);
    }

    void join() {
        if (!joined)
            join_data_vectors(values);
        joined = true;
    }

#if !(defined(CUDA) || defined(HIP)) || defined(HILAPP)

    void endloop_action() {
        bool save = auto_join;
        auto_join = false;
        SiteSelect::endloop_action();
        values.resize(current_index);
        auto_join = save;
        if (auto_join)
            join();
    }

#else
    // skip this for hilapp
    void endloop_action(const char *flag, const SiteIndex *d_sites, const T *d_values) {
        copy_data_to_host_vector(sites, flag, d_sites);
        copy_data_to_host_vector(values, flag, d_values);

        if (auto_join)
            join();
    }

#endif // GPU
};


#ifdef HILAPP

// Make hilapp generate __device__ versions of SiteIndex function - this is removed in final program

inline void dummy_func_2() {
    onsites(ALL) {
        auto s = SiteIndex(X.coordinates());
    }
}

#endif


#endif
