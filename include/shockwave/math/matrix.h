#pragma once

#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <cmath>
#include <utility>

#include "shockwave/utility/types.h"
#include "shockwave/utility/repeat.h"


namespace sw
{

////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \brief Matrix data, stored row-major
///
/// Two methods for random access are provided,
///     For matrices,  mat(i, j) is recommended
///     For vectors,   vec[i]    is recommended
///
/// Direct access to the data member is allowed, but for matrices
///     it is preferable to use operator() to avoid errors
////////////////////////////////////////////////////////////
template<class T, u32 n_rows_, u32 m_cols_>
struct MatrixData
{
    typedef T        value_type;
    typedef T*       iterator;
    typedef const T* const_iterator;

    static constexpr u32 n_rows = n_rows_;
    static constexpr u32 m_cols = m_cols_;

    constexpr MatrixData() = default;

    template<class... Args>
    constexpr MatrixData(Args&&... values)
        : MatrixData({static_cast<T>(std::forward<Args>(values))...}) {}

    constexpr MatrixData(std::initializer_list<T> values) {
        for ( u32 i = 0; i < size(); i++ ) {
            data[i] = i < values.size() ? *(values.begin() + i) : T{};
        }
    }

    template<class U>
    explicit constexpr MatrixData(const MatrixData<U, n_rows, m_cols>& other) {
        repeat<size()>([&](u32 i) { data[i] = static_cast<T>(other.data[i]); });
    }

    constexpr static u32 size() { return n_rows * m_cols; }

    T&       operator()(u32 row, u32 col) { return data[row * m_cols + col]; }
    const T& operator()(u32 row, u32 col) const {
        return data[row * m_cols + col];
    }

    T&       operator[](u32 index) { return data[index]; }
    const T& operator[](u32 index) const { return data[index]; }


    iterator begin() { return data; }
    iterator end() { return data + size(); }

    const_iterator begin() const { return data; }
    const_iterator end() const { return data + size(); }

    T data[size()];
};


template<class, u32, u32>
struct Matrix;

////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \brief Vector class
///
/// Vectors are always column vectors,
///     use transpose(vec) to obtain a row vector when needed.
////////////////////////////////////////////////////////////
template<class T, u32 dim>
using Vector = Matrix<T, dim, 1>;

template<class T, u32 m, u32 n, bool isSquare = (m == n)>
struct SpecialConstructors
{
};

template<class T, u32 dim>
struct SpecialConstructors<T, dim, dim, true>
{
    static inline Matrix<T, dim, dim> identity();
    static inline Matrix<T, dim, dim> diagonal(const Vector<T, dim>& elements);
};

template<class T, u32 dim>
struct SpecialConstructors<T, dim, 1, false>
{
    static inline Vector<T, dim> i();
    static inline Vector<T, dim> j();
    static inline Vector<T, dim> k();
    static inline Vector<T, dim> w();
};

////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \brief Matrix class for small matrices
///
/// Most operations are provided as global functions
/// \see operations
////////////////////////////////////////////////////////////
template<class T, u32 m, u32 n>
struct Matrix : public MatrixData<T, m, n>, public SpecialConstructors<T, m, n>
{
    Matrix() = default;

    template<class U>
    explicit inline Matrix(const Matrix<U, m, n>& other);

    explicit inline Matrix(const MatrixData<T, m, n>& data);

    explicit inline Matrix(const std::initializer_list<T>& init);

    static inline Matrix filled(T val);

    template<class U>
    static inline Matrix fromRows(
        const std::initializer_list<Vector<U, n>>& rows);

    template<class U>
    static inline Matrix fromRows(const Vector<Vector<U, n>, m>& rows);

    template<class U>
    static inline Matrix fromCols(
        const std::initializer_list<Vector<U, m>>& cols);

    template<class U>
    static inline Matrix fromCols(const Vector<Vector<U, m>, n>& cols);

    inline Vector<Vector<T, n>, m> asRows() const;
    inline Vector<Vector<T, m>, n> asCols() const;


    inline Matrix operator+() const;
    inline Matrix operator-() const;

    template<class U>
    inline Matrix& operator+=(const Matrix<U, m, n>& other);

    template<class U>
    inline Matrix& operator-=(const Matrix<U, m, n>& other);

    template<class U>
    inline Matrix& operator*=(U scalar);

    template<class U>
    inline Matrix& operator/=(U scalar);
};

////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \defgroup operations
/// \{
////////////////////////////////////////////////////////////

template<class T, class U>
using Promoted = typename std::common_type<T, U>::type;


template<class T, class U, u32 m, u32 n>
inline Matrix<Promoted<T, U>, m, n> operator+(const Matrix<T, m, n>& lhs,
                                              const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline Matrix<Promoted<T, U>, m, n> operator-(const Matrix<T, m, n>& lhs,
                                              const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline Matrix<Promoted<T, U>, m, n> operator*(U                      scalar,
                                              const Matrix<T, m, n>& mat);

template<class T, class U, u32 m, u32 n>
inline Matrix<Promoted<T, U>, m, n> operator*(const Matrix<T, m, n>& mat,
                                              U                      scalar);

template<class T, class U, u32 m, u32 n>
inline Matrix<Promoted<T, U>, m, n> operator/(const Matrix<T, m, n>& mat,
                                              U                      scalar);

template<class T, class U, u32 mLeft, u32 nLeft, u32 nRight>
inline Matrix<Promoted<T, U>, mLeft, nRight> operator*(
    const Matrix<T, mLeft, nLeft>&  lhs,
    const Matrix<U, nLeft, nRight>& rhs);


template<class T, u32 m, u32 n>
inline Matrix<T, n, m> transpose(const Matrix<T, m, n>& original);


////////////////////////////////////////////////////////////
/// \brief Return x such that Ax = b
///
////////////////////////////////////////////////////////////
template<class T, class U, u32 n>
inline Vector<Promoted<T, U>, n> solve(const Matrix<T, n, n>& A,
                                       const Vector<U, n>&    b);

template<class T, u32 n>
class Solver
{
public:

    inline Solver(const Matrix<T, n, n>& A);

    template<class U>
    inline Vector<Promoted<T, U>, n> solve(const Vector<U, n>& b) const;

    Matrix<T, n, n> original() const { return transpose(QT_) * R_; }

    bool isValid() const { return isValid_; }

private:

    Matrix<T, n, n> QT_;
    Matrix<T, n, n> R_;
    bool            isValid_ = true;
};


template<class T, u32 n>
inline Matrix<T, n, n> invert(const Matrix<T, n, n>& mat);


////////////////////////////////////////////////////////////
/// \brief Return x such that Ax >= b with bounds on x
///
/// proj must project x onto its valid bounds (ie clamp its values). With the
/// following signature:
///     T proj(const Vector<T, n>& x, u32 index)
/// where x[index] must be clamped and returned.
///
/// when the absolute change of components of x drops below epsilon, iteration
/// terminates.
///
/// Uses projected Gauss-Seidel to solve the MLCP,
/// See A. Enzenhofer's master thesis (McGill): Numerical Solutions of MLCP
////////////////////////////////////////////////////////////
template<class T, u32 n, std::invocable<Vector<T, n>, u32> Proj>
inline Vector<T, n> solveInequalities(const Matrix<T, n, n>& A,
                                      Vector<T, n>           b,
                                      Proj                   proj,
                                      Vector<T, n>           initialGuess
                                      = Vector<T, n>{},
                                      float epsilon = sw::EPSILON);

////////////////////////////////////////////////////////////
/// \brief Return x such that Ax >= b with bounds on x
///
/// proj must project x onto its valid bounds (ie clamp its values).
///
/// when the absolute change of components of x drops below epsilon, iteration
/// terminates.
///
/// Uses projected Gauss-Seidel to solve the MLCP,
/// See A. Enzenhofer's master thesis (McGill): Numerical Solutions of MLCP
////////////////////////////////////////////////////////////
template<class T, u32 n, std::invocable<Vector<T, n>> Proj>
inline Vector<T, n> solveInequalities(const Matrix<T, n, n>& A,
                                      Vector<T, n>           b,
                                      Proj                   proj,
                                      Vector<T, n>           initialGuess
                                      = Vector<T, n>{},
                                      float epsilon = sw::EPSILON);


////////////////////////////////////////////////////////////
/// \brief Solves the LCP Ax >= b with x >= 0
///
/// This gives Ax - b = w with the residuals w >= 0,
///     the complementarity condition is dot(x, w) = 0
///     (xi = 0 or wi = 0 for each index i)
///
/// Uses total enumeration, taken from box2d's contact solver.
/// If no value is returned, the LCP had no solution.
///
////////////////////////////////////////////////////////////
template<class T>
inline Vector<T, 2> solveLcp(const Matrix<T, 2, 2>& A, const Vector<T, 2>& b);

template<class T>
class LcpSolver;


////////////////////////////////////////////////////////////
// Vector operations
////////////////////////////////////////////////////////////

template<class T, class U, u32 dim>
inline Vector<Promoted<T, U>, dim> elementWiseMul(const Vector<T, dim>& lhs,
                                                  const Vector<U, dim>& rhs) {
    Vector<Promoted<T, U>, dim> res;
    for ( u32 i = 0; i < dim; ++i )
        res[i] = lhs[i] * rhs[i];

    return res;
}

template<class T, class U, u32 dim>
inline Promoted<T, U> dot(const Vector<T, dim>& lhs, const Vector<U, dim>& rhs);

template<class T, class U>
inline Vector<Promoted<T, U>, 3> cross(const Vector<T, 3>& lhs,
                                       const Vector<U, 3>& rhs);

template<class T, class U>
inline Promoted<T, U> cross(const Vector<T, 2>& lhs, const Vector<U, 2>& rhs);

template<class T, u32 dim>
inline T normSquared(const Vector<T, dim>& v);

template<class T, u32 dim>
inline T norm(const Vector<T, dim>& v);

template<class T, u32 dim>
inline Vector<T, dim> normalized(const Vector<T, dim>& v);

template<class T, class U, u32 dim>
inline Vector<Promoted<T, U>, dim> projection(const Vector<T, dim>& ofThis,
                                              const Vector<U, dim>& onThat);


////////////////////////////////////////////////////////////
/// \brief rotates v by 90 degrees
///
/// If clockwise is false, then this is (k x v)
/// If clockwise is true, then this is  (v x k)
///
////////////////////////////////////////////////////////////
template<class T>
inline Vector<T, 2> perp(const Vector<T, 2>& v, bool clockwise = false);

/// \}


////////////////////////////////////////////////////////////
// Aliases
////////////////////////////////////////////////////////////

typedef Vector<float, 2> Vec2;
typedef Vector<float, 3> Vec3;
typedef Vector<float, 4> Vec4;

typedef Vector<i32, 2> Vec2i;
typedef Vector<i32, 3> Vec3i;
typedef Vector<i32, 4> Vec4i;

typedef Matrix<float, 2, 2> Mat2;
typedef Matrix<float, 3, 3> Mat3;
typedef Matrix<float, 4, 4> Mat4;

typedef Matrix<i32, 2, 2> Mat2i;
typedef Matrix<i32, 3, 3> Mat3i;
typedef Matrix<i32, 4, 4> Mat4i;


////////////////////////////////////////////////////////////
// Comparison facilities
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \defgroup comparison
/// \{
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
/// \brief Comparison Matrix returned by comparisons of matrices
///
/// This is not convertible to bool, use the functions all() and any()
///     to evaluate to a bool.
///
/// All matrix comparisons are done element-wise.
///
/// element-wise boolean operators are also provided for ComparisonMatrix
////////////////////////////////////////////////////////////
template<u32 m, u32 n>
using ComparisonMatrix = MatrixData<bool, m, n>;


template<class T, class U, u32 m, u32 n>
inline ComparisonMatrix<m, n> operator==(const Matrix<T, m, n>& lhs,
                                         const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline ComparisonMatrix<m, n> operator!=(const Matrix<T, m, n>& lhs,
                                         const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline ComparisonMatrix<m, n> operator<(const Matrix<T, m, n>& lhs,
                                        const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline ComparisonMatrix<m, n> operator<=(const Matrix<T, m, n>& lhs,
                                         const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline ComparisonMatrix<m, n> operator>(const Matrix<T, m, n>& lhs,
                                        const Matrix<U, m, n>& rhs);

template<class T, class U, u32 m, u32 n>
inline ComparisonMatrix<m, n> operator>=(const Matrix<T, m, n>& lhs,
                                         const Matrix<U, m, n>& rhs);


////////////////////////////////////////////////////////////
/// \brief true if all elements are true
///
////////////////////////////////////////////////////////////
template<u32 m, u32 n>
inline bool all(const ComparisonMatrix<m, n>& comp);

////////////////////////////////////////////////////////////
/// \brief true if any element is true
///
////////////////////////////////////////////////////////////
template<u32 m, u32 n>
inline bool any(const ComparisonMatrix<m, n>& comp);

template<u32 m, u32 n>
inline ComparisonMatrix<m, n> operator&&(const ComparisonMatrix<m, n>& lhs,
                                         const ComparisonMatrix<m, n>& rhs);

template<u32 m, u32 n>
inline ComparisonMatrix<m, n> operator||(const ComparisonMatrix<m, n>& lhs,
                                         const ComparisonMatrix<m, n>& rhs);

template<u32 m, u32 n>
inline ComparisonMatrix<m, n> operator!(const ComparisonMatrix<m, n>& unary);

/// \}

}  // namespace sw


////////////////////////////////////////////////////////////
// std overloads
////////////////////////////////////////////////////////////

namespace std
{

////////////////////////////////////////////////////////////
/// \ingroup operations
/// \{
////////////////////////////////////////////////////////////

template<class T, class U, sw::u32 m, sw::u32 n>
struct common_type<sw::Matrix<T, m, n>, sw::Matrix<U, m, n>>
{
    typedef sw::Matrix<typename std::common_type<T, U>::type, m, n> type;
};

template<class T, sw::u32 m, sw::u32 n>
inline sw::Matrix<T, m, n> abs(const sw::Matrix<T, m, n>& mat);

template<class T, sw::u32 m, sw::u32 n>
inline sw::Matrix<T, m, n> round(const sw::Matrix<T, m, n>& mat);

template<class T, sw::u32 m, sw::u32 n>
inline sw::Matrix<T, m, n> min(const sw::Matrix<T, m, n>& lhs,
                               const sw::Matrix<T, m, n>& rhs);

template<class T, sw::u32 m, sw::u32 n>
inline sw::Matrix<T, m, n> max(const sw::Matrix<T, m, n>& lhs,
                               const sw::Matrix<T, m, n>& rhs);

/// \}

}  // namespace std
