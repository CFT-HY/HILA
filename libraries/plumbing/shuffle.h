
#ifndef HILA_SHUFFLE_H_
#define HILA_SHUFFLE_H_

#include "hila.h"

//

namespace hila {
struct dir_and_parity {
    Direction direction;
    Parity parity;
};


//////////////////////////////////////////////////////////////////////////
/// Shuffle an existing std::array or std::vector to random order
/// The type T should be trivially copyable, and "lightweight" - does extra copies
///

template <typename T,
          std::enable_if_t<hila::is_std_array<T>::value || hila::is_std_vector<T>::value, int> = 0>
void shuffle(T &arr) {
    // go backwards in array, swap with random earlier element incl. itself
    for (int j = arr.size() - 1; j > 0; j--) {
        int i = (j + 1) * hila::random(); // rand from 0 to j inclusive
        if (i != j)
            std::swap(arr[i], arr[j]);
    }
}

///////////////////////////////////////////////////////////////////////////////
/// Function returning all directions and parities as shuffled to random order.
/// Return type is std::array<dir_and_parity,2*NDIM>,
///
/// Use case is in Gauge heatbath/overrelax updates, to randomize order.
/// Synchronizes with all MPI nodes, i.e. all ranks will have the same content.
/// Function marked inline in order to avoid ODR
///
/// Typical use: no need to declare std::array
///
///    for (const auto & s : hila::shuffle_directions_and_parities()) {
///        update_parity_dir(U, s.parity, s.direction);
///    }
///
/// Can also be assigned, i.e.
///    auto shuffled = hila::shuffle_directions_and_parities();

inline auto shuffle_directions_and_parities() {
    std::array<struct dir_and_parity, 2 * NDIM> arr;

    int i = 0;
    foralldir(d) {
        for (Parity p : {EVEN, ODD}) {
            arr[i].direction = d;
            arr[i].parity = p;
            i++;
        }
    }
    shuffle(arr);
    hila::broadcast(arr); // sync with all nodes
    return arr;
}

///////////////////////////////////////////////////////////////////////////////
/// Shuffling directions and parities separately with functions
///
/// std::array<Direction,NDIM> hila::shuffle_directions()
/// std::array<Parity,2> hila::shuffle_parities()
///
/// Thus, (almost) the same operation as above with shuffle_directions_and_parities() 
/// can be done with
///
///    for (const Direction d : hila::shuffle_directions()) {
///        for (const Parity p : hila::shuffle_parities()) {
///            update_parity_dir(U, p, d);
///        }
///    }
///
/// This does more MPI synchronization than the shuffle_directions_and_parities()
///

inline auto shuffle_directions() {
    std::array<Direction, NDIM> arr;

    int i = 0;
    foralldir(d) {
        arr[i] = d;
        i++;
    }

    shuffle(arr);
    hila::broadcast(arr); // sync with all nodes
    return arr;
}

inline auto shuffle_parities() {
    std::array<Parity, 2> arr;

    arr[0] = EVEN;
    arr[1] = ODD;

    shuffle(arr);
    hila::broadcast(arr); // sync with all nodes
    return arr;
}


} // namespace hila

#endif
