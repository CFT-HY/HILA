#ifndef VECTOR_TYPES_H
#define VECTOR_TYPES_H

namespace hila {

/// is_vectorizable_type<T>::value  is always false if the target is not vectorizable
template <typename T, typename A = void> struct is_vectorizable_type {
    static constexpr bool value = false;
};

#ifdef VECTORIZED
///  specialize is_vectorizable_type<T>::value to true if the base_type_struct<T>::type
///  exists
template <typename T>
struct is_vectorizable_type<T, typename ::std::enable_if_t<::std::is_arithmetic<
                                   typename hila::base_type_struct<T>::type>::value>> {
    static constexpr bool value = true;
};

/// do forward definition here, enables inclusion
template <int vector_size> struct vectorized_lattice_struct;

/// Utility for selecting a vector type by base type and length
template <typename T, int vector_len> struct vector_base_type {};

template <> struct vector_base_type<double, 4> { using type = Vec4d; };

template <> struct vector_base_type<double, 8> { using type = Vec8d; };

template <> struct vector_base_type<float, 8> { using type = Vec8f; };

template <> struct vector_base_type<float, 16> { using type = Vec16f; };

template <> struct vector_base_type<int, 4> { using type = Vec4i; };

template <> struct vector_base_type<int, 8> { using type = Vec8i; };

template <> struct vector_base_type<int, 16> { using type = Vec16i; };

template <> struct vector_base_type<int64_t, 4> { using type = Vec4q; };

template <> struct vector_base_type<int64_t, 8> { using type = Vec8q; };

// template<>
// struct vector_base_type<CoordinateVector, 4> {
//   using type = Vec4i;
// };

// template<>
// struct vector_base_type<CoordinateVector, 8> {
//   using type = Vec8i;
// };

#endif // VECTORIZED

/// Construct the vector info for the type.
/// first, if the type is not vectorizable
template <typename T, typename A = void> struct vector_info {
    static constexpr bool is_vectorizable = false;
    using base_type = T;
    // Find vector length
    static constexpr int vector_size = 1;
    // Find the vector type from above
    using type = void;
    // Number of elements in the full type
    static constexpr int elements = 1;
    // Size of the base type
    static constexpr int base_type_size = sizeof(base_type);
};

#ifdef VECTORIZED

/// and specializre the same for vectorizable type
template <typename T>
struct vector_info<T, typename ::std::enable_if_t<::std::is_arithmetic<
                          typename hila::base_type_struct<T>::type>::value>> {
    static constexpr bool is_vectorizable = true;
    // Get base type first
    using base_type = hila::number_type<T>;
    // Find vector length
    static constexpr int vector_size = VECTOR_SIZE / sizeof(base_type);
    // Find the vector type from above
    using type = typename vector_base_type<base_type, vector_size>::type;
    // Number of elements in the full type
    static constexpr int elements = sizeof(T) / sizeof(base_type);
    // Size of the base type
    static constexpr int base_type_size = sizeof(base_type);
};

#endif // VECTORIZED

} // namespace hila

#endif