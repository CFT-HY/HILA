#ifndef SITE_SELECT_H_
#define SITE_SELECT_H_

#include "hila.h"

// MPI is needed here


//////////////////////////////////////////////////////////////////////////////////
/// Site selection: special vector to accumulate chosen sites or sites + variable
///
/// SiteSelect<> s;
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
    // Define iterators using std::vector iterators
    using iterator = typename std::vector<SiteIndex>::iterator;
    using const_iterator = typename std::vector<SiteIndex>::const_iterator;

    iterator begin() {
        return sites.begin();
    }
    iterator end() {
        return sites.end();
    }
    const_iterator begin() const {
        return sites.begin();
    }
    const_iterator end() const {
        return sites.end();
    }

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
        sites.resize(lattice.mynode.volume());
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
        values.resize(lattice.mynode.volume());
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

    void endloop_action() {
        bool save = auto_join;
        auto_join = false;
        SiteSelect::endloop_action();
        values.resize(current_index);
        auto_join = save;
        if (auto_join)
            join();
    }

    void join() {
        if (!joined)
            join_data_vectors(values);
        joined = true;
    }
};


#endif
