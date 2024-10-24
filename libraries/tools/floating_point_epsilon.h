/** @file floating_point_epsilon.h */

#ifndef FLOATING_POINT_EPSILON_H_
#define FLOATING_POINT_EPSILON_H_

template <typename T>
struct fp {
    static constexpr T epsilon=0;
};

template <>
struct fp<long double> {
    static constexpr long double epsilon=1.0842e-19;
};

template <>
struct fp<double> {
    static constexpr double epsilon=2.22045e-16;
};

template <>
struct fp<float> {
    static constexpr float epsilon=1.19209e-07f;
};

#endif