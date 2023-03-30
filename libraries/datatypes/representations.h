#ifndef REPRESENTATIONS_H_
#define REPRESENTATIONS_H_

#include "sun.h"

/// A matrix in the adjoint representation of the SU(N) group
///
/// Members
/// adjoint.represent(sun m): projects the sU(N) matrix to the
///     adjoint representation and replaces this
///
/// Class functions:
/// adjoint::generator(int i): returns generators of SU(N)
///     in the fundamental representation
/// adjoint::represented_generator_I(int i): returns the
///     generator (times I) in the adjoint representation
/// adjoint::project_force(squarematrix): projects a square
///     matrix the size of an adjoint matrix to the SU(N)
///     algebra. This is used to calculate the force of an
///     adjoint action term to a derivative of the underlying
///     su(N) group
///
template <int N, typename radix>
class adjointRep : public SquareMatrix<N * N - 1, radix> {
  public:
    static_assert(hila::is_arithmetic<radix>::value, "adjointRep<type>: type has to be real");

    /// The underlying arithmetic type of the matrix
    using base_type = hila::scalar_type<radix>;
    using argument_type = radix;
    /// The SU(N) type the adjoint matrix is constructed of
    using sun = SU<N, radix>;

    /// Matrix size
    constexpr static int size = N * N - 1;
    /// Info on group generators
    static sun generator(int i) { return sun::generator(i); }

    /// std ops required for triviality
    adjointRep() = default;
    /// std ops required for triviality
    ~adjointRep() = default;
    /// std ops required for triviality
    adjointRep(const adjointRep &a) = default;

    /// Use square matrix constructor from radix
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    adjointRep(const scalart m) : SquareMatrix<size, radix>(m) {}

    /// Copy constructor
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    adjointRep(const adjointRep<N, scalart> m) : SquareMatrix<size, scalart>(m) {}

    /// Automatic conversion from SquareMatrix is needed!
    adjointRep(const SquareMatrix<size,radix> m) : SquareMatrix<size,radix>(m) {}


    /// Return a SU(N) generator in the adjoint representation,
    /// multiplied by I
    static adjointRep represented_generator_I(int i) {
        static bool initialize = true;
        static adjointRep r_generators[size];
        if (initialize)
            for (int g = 0; g < size; g++) {
                r_generators[g] = 0;
                sun tg = sun::generator(g);
                for (int j = 0; j < size; j++) {
                    sun tj = generator(j);
                    for (int k = 0; k < size; k++) {
                        sun tk = generator(k);
                        Complex<radix> tr1 = (tg * tj * tk).trace();
                        Complex<radix> tr2 = (tj * tg * tk).trace();
                        r_generators[g].e(j, k) = 2 * (tr1.im - tr2.im);
                    }
                }
                initialize = false;
            }
        return r_generators[i];
    }

    /// Project a matrix into the adjoint representation
    void represent(sun &m) {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                (*this).e(i, j) =
                    2 * (m.adjoint() * generator(i) * m * generator(j)).trace().re;
            }
    }

    /// Project a complex adjoint matrix into the algebra and
    /// represent as a complex NxN (momentum) matrix
    static SquareMatrix<N, Complex<radix>>
    project_force(SquareMatrix<size, Complex<radix>> rforce) {
        SquareMatrix<N, Complex<radix>> fforce = 0;
        for (int g = 0; g < size; g++) {
            adjointRep rg = represented_generator_I(g);
            radix C = (rg.transpose() * rforce).trace().re;
            fforce += C * sun::generator(g);
        }
        Complex<radix> ct(0, 2.0);
        fforce = fforce * ct;
        project_antihermitean(fforce);
        return fforce;
    }
};

/// A matrix in the antisymmetric representation of the SU(N) group
///
/// Members
/// antisymmetric.represent(sun m): projects the sU(N) matrix to the
///     antisymmetric representation and replaces this
///
/// Class functions:
/// antisymmetric::generator(int i): returns antisymmetric matrices
///     in the fundamental representation
/// antisymmetric::represented_generator_I(int i): returns
///     antisymmetric SU(N) matrices (times I) in the antisymmetric
///     representation
/// antisymmetric::project_force(squarematrix): projects a square
///     matrix the size of an antisymmetric matrix to the SU(N)
///     algebra. This is used to calculate the force of an
///     antisymmetric action term to a derivative of the underlying
///     su(N) group
///
template <int N, typename radix>
class antisymmetric : public SquareMatrix<N *(N - 1) / 2, Complex<radix>> {
  public:
    /// The underlying arithmetic type of the matrix
    using base_type = hila::scalar_type<radix>;
    using argument_type = radix;
    /// The SU(N) type the adjoint matrix is constructed of
    using sun = SU<N, radix>;

    /// Matrix size
    constexpr static int size = N * (N - 1) / 2;

    /// Use square matrix constructors
    using SquareMatrix<size, Complex<radix>>::SquareMatrix;
    /// Use square matrix constructors
    using SquareMatrix<size, Complex<radix>>::operator=;

    /// default constructor
    antisymmetric() = default;

    /// Square matrix constructor from scalar
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    antisymmetric(const scalart m) : SquareMatrix<size, Complex<radix>>(m) {}

    /// Copy constructor
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    antisymmetric(const antisymmetric<N, scalart> m) {
        for (int j = 0; j < size; j++)
            for (int i = 0; i < size; i++) {
                this->e(i, j) = m.e(i, j);
            }
    }

    /// automatic conversion from SquareMatrix -- needed for methods!
    antisymmetric(const SquareMatrix<size,Complex<radix>> m) : SquareMatrix<size,Complex<radix>>(m) {}

    /// Needs assignment as well
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    inline antisymmetric &operator=(const antisymmetric<N, scalart> m) {
        for (int j = 0; j < size; j++)
            for (int i = 0; i < size; i++) {
                this->e(i, j) = m.e(i, j);
            }
        return *this;
    }

    /// The antisymmetric generator
    static sun generator(int ng) {
        static bool initialize = true;
        static sun generators[size];
        if (initialize)
            for (int g = 0; g < size; g++) {
                generators[g] = 0;
                int k = 0;
                for (int m1 = 0; m1 < N; m1++)
                    for (int m2 = m1 + 1; m2 < N; m2++) {
                        if (g == k) {
                            generators[g].e(m1, m2).re = 0.5;
                            generators[g].e(m2, m1).re = -0.5;
                        }
                        k++;
                    }
                initialize = false;
            }
        return generators[ng];
    }

    /// Return a SU(N) generator (times I) in the antisymmetric representation
    static antisymmetric represented_generator_I(int i) {
        static bool initialize = true;
        static antisymmetric r_generators[size];
        if (initialize)
            for (int g = 0; g < size; g++) {
                r_generators[g] = 0;
                sun tg = sun::generator(g);
                for (int j = 0; j < N * (N - 1) / 2; j++) {
                    sun tj = generator(j);
                    for (int k = 0; k < N * (N - 1) / 2; k++) {
                        sun tk = generator(k);

                        Complex<radix> tr = (tj * tg * tk).trace();
                        r_generators[g].e(j, k) = Complex<radix>(0, 4) * tr;
                    }
                }
                initialize = false;
            }
        return r_generators[i];
    }

    /// Project a matrix into the antisymmetric representation
    void represent(sun &m) {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                (*this).e(i, j) =
                    2 *
                    (generator(i) * m * generator(j).adjoint() * m.transpose()).trace();
            }
    }

    /// Project a complex antisymmetric matrix into the algebra and
    /// represent as a complex NxN (momentum) matrix
    static SquareMatrix<N, Complex<radix>>
    project_force(SquareMatrix<size, Complex<radix>> rforce) {
        SquareMatrix<N, Complex<radix>> fforce = 0;
        for (int g = 0; g < size; g++) {
            antisymmetric rg = represented_generator_I(g);
            radix C = (rg.adjoint() * rforce).trace().re;
            fforce += C * sun::generator(g);
        }
        Complex<radix> ct(0, -2.0);
        fforce = fforce * ct;
        project_antihermitean(fforce);
        return fforce;
    }
};

/// A matrix in the symmetric representation of the SU(N) group
///
/// Members
/// symmetric.represent(sun m): projects the sU(N) matrix to the
///     symmetric representation and replaces this
///
/// Class functions:
/// symmetric::generator(int i): returns symmetric matrices
///     in the fundamental representation
/// symmetric::represented_generator_I(int i): returns
///     symmetric SU(N) matrices (times I) in the symmetric
///     representation
/// symmetric::project_force(squarematrix): projects a square
///     matrix the size of an symmetric matrix to the SU(N)
///     algebra. This is used to calculate the force of an
///     symmetric action term to a derivative of the underlying
///     su(N) group
///
template <int N, typename radix>
class symmetric : public SquareMatrix<N *(N + 1) / 2, Complex<radix>> {
  public:
    /// The underlying arithmetic type of the matrix
    using base_type = hila::scalar_type<radix>;
    using argument_type = radix;
    /// The SU(N) type the adjoint matrix is constructed of
    using sun = SU<N, radix>;

    /// Use square matrix constructors
    using SquareMatrix<N *(N + 1) / 2, Complex<radix>>::SquareMatrix;
    /// Use square matrix constructors
    using SquareMatrix<N *(N + 1) / 2, Complex<radix>>::operator=;

    /// Matrix size
    constexpr static int size = N * (N + 1) / 2;

    /// Use default constructor
    symmetric() = default;

    /// Constructor from scalar
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    symmetric(const scalart m) {
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                this->e(i, j) = 0;
            }
            this->c[j][j] = m;
        }
    }

    /// Constructor from a symmetric matrix
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    symmetric(const symmetric<N, scalart> m) {
        for (int j = 0; j < size; j++)
            for (int i = 0; i < size; i++) {
                this->e(i, j) = m.e(i, j);
            }
    }

    /// Needs assignment as well
    template <typename scalart, std::enable_if_t<hila::is_arithmetic<scalart>::value, int> = 0>
    inline symmetric &operator=(const symmetric<N, scalart> m) {
        for (int j = 0; j < size; j++)
            for (int i = 0; i < size; i++) {
                this->e(i, j) = m.e(i, j);
            }
        return *this;
    }

    /// Symmetric generators as SU(N) matrices
    static sun generator(int ng) {
        static bool initialize = true;
        static sun generators[size];
        if (initialize)
            for (int g = 0; g < size; g++) {
                generators[g] = 0;
                if (g < N) {
                    generators[g].e(g, g).re = sqrt(0.5);
                }
                int k = N;
                for (int m1 = 0; m1 < N; m1++)
                    for (int m2 = m1 + 1; m2 < N; m2++) {
                        if (g == k) {
                            generators[g].e(m1, m2).re = 0.5;
                            generators[g].e(m2, m1).re = 0.5;
                        }
                        k++;
                    }
                initialize = false;
            }
        return generators[ng];
    }

    /// Return a symmetric generators (times I) in the symmetric representation
    static symmetric represented_generator_I(int i) {
        static bool initialize = true;
        static symmetric r_generators[size];
        if (initialize)
            for (int g = 0; g < size; g++) {
                r_generators[g] = 0;
                sun tg = sun::generator(g);
                for (int j = 0; j < N * (N + 1) / 2; j++) {
                    sun tj = generator(j);
                    for (int k = 0; k < N * (N + 1) / 2; k++) {
                        sun tk = generator(k);

                        Complex<radix> tr = (tj * tg * tk).trace();
                        r_generators[g].e(j, k) = Complex<radix>(0, 4) * tr;
                    }
                }
                initialize = false;
            }
        return r_generators[i];
    }

    /// Project a matrix into the symmetric representation
    void represent(sun &m) {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++) {
                (*this).e(i, j) =
                    2 * (generator(i) * m * generator(j) * m.transpose()).trace();
            }
    }

    /// Project a complex symmetric matrix into the algebra and
    /// represent as a complex NxN (momentum) matrix
    static SquareMatrix<N, Complex<radix>>
    project_force(SquareMatrix<size, Complex<radix>> rforce) {
        SquareMatrix<N, Complex<radix>> fforce = 0;
        for (int g = 0; g < size; g++) {
            symmetric rg = represented_generator_I(g);
            radix C = (rg * rforce).trace().re;
            fforce += C * sun::generator(g);
        }
        Complex<radix> ct(0, -2.0);
        fforce = fforce * ct;
        project_antihermitean(fforce);
        return fforce;
    }
};

#endif