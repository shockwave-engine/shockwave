#pragma once

#include <algorithm>
#include <compare>
#include <functional>
#include <initializer_list>
#include <type_traits>
#include <cmath>
#include <utility>

#include "shockwave/utility/types.h"


namespace sw
{

template<class, u32, u32>
struct Matrix;

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
template<class T,
         u32 m_rows_,
         u32 n_cols_,
         template<class, u32, u32>
         class Derived>
struct MatrixData
{
    typedef T        value_type;
    typedef T*       iterator;
    typedef const T* const_iterator;

    static constexpr u32 m_rows = m_rows_;
    static constexpr u32 n_cols = n_cols_;

    template<class U = T, u32 m = m_rows, u32 n = n_cols>
    using Result = Derived<std::common_type_t<T, U>, m, n>;

    template<class U = T, u32 m = m_rows, u32 n = n_cols>
    using Mat = Derived<U, m, n>;

    using Self = Mat<>;

    using Bool = Mat<bool>;


    ////////////////////////////////////////////////////////////
    // Constructors
    ////////////////////////////////////////////////////////////

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
    explicit constexpr MatrixData(const Mat<U>& other) {
        for ( u32 i = 0; i < size(); i++ ) {
            data[i] = static_cast<T>(other.data[i]);
        }
    }


    ////////////////////////////////////////////////////////////
    // Indexing
    ////////////////////////////////////////////////////////////

    constexpr static u32 size() { return m_rows * n_cols; }

    T&       operator()(u32 row, u32 col) { return data[row * n_cols + col]; }
    const T& operator()(u32 row, u32 col) const {
        return data[row * n_cols + col];
    }

    T&       operator[](u32 index) { return data[index]; }
    const T& operator[](u32 index) const { return data[index]; }

    iterator begin() { return data; }
    iterator end() { return data + size(); }

    const_iterator begin() const { return data; }
    const_iterator end() const { return data + size(); }

    ////////////////////////////////////////////////////////////
    // Operations
    ////////////////////////////////////////////////////////////

    template<std::invocable<const T&> Fn>
    static Mat<std::result_of<Fn(const T&)>> unary_element_wise(const Self& mat,
                                                                Fn&& fn) {
        Mat<std::result_of<Fn(const T&)>> res{mat};
        for ( T& value : res ) {
            value = fn(value);
        }

        return res;
    }

    template<class U, std::invocable<const T&> Fn>
    static Mat<std::result_of<Fn(const T&, const U&)>>
    element_wise(const Self& lhs, const Mat<U>& rhs, Fn&& fn) {
        Mat<std::result_of<Fn(const T&, const U&)>> res;

        for ( u32 i = 0; i < lhs.size(); i++ ) {
            res[i] = fn(lhs[i], rhs[i]);
        }

        return res;
    }


    T data[size()];
};

template<class T,
         u32 m_rows_,
         u32 n_cols_,
         template<class, u32, u32>
         class Derived>
struct MatrixOperations : public MatrixData<T, m_rows_, n_cols_, Derived>
{
    using Base = MatrixData<T, m_rows_, n_cols_, Derived>;
    using Base::Base;

    template<class U = T, u32 m = Base::m_rows, u32 n = Base::n_cols>
    using Result = Derived<std::common_type_t<T, U>, m, n>;

    template<class U = T, u32 m = Base::m_rows, u32 n = Base::n_cols>
    using Mat = Derived<U, m, n>;

    using Self = Mat<>;

    using Bool = Mat<bool>;

    ////////////////////////////////////////////////////////////
    // Unary operators
    ////////////////////////////////////////////////////////////

    Result<> operator+() const { return *this; }
    Result<> operator-() const {
        return Base::unary_element_wise(*this, std::negate<>());
    }

    ////////////////////////////////////////////////////////////
    // Binary operators
    ////////////////////////////////////////////////////////////

    template<class U>
    friend Result<U> operator+(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, std::plus<>());
    }

    template<class U>
    friend Result<U> operator-(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, std::minus<>());
    }

    template<class U>
    friend Result<U> operator*(U scalar, const Self& mat) {
        return Base::unary_element_wise(mat,
                                        [=](const T& x) { return x * scalar; });
    }

    template<class U>
    friend Result<U> operator*(const Self& mat, U scalar) {
        return scalar * mat;
    }

    template<class U>
    friend Result<U> operator/(const Self& mat, U scalar) {
        return Base::unary_element_wise(mat,
                                        [=](const T& x) { return x / scalar; });
    }

    ////////////////////////////////////////////////////////////
    // Compound assignement operators
    ////////////////////////////////////////////////////////////

    template<class U>
    Self& operator+=(const Mat<U>& other) {
        return *this = *this + other;
    }

    template<class U>
    Self& operator-=(const Mat<U>& other) {
        return *this = *this - other;
    }

    template<class U>
    Self& operator*=(U scalar) {
        return *this = *this * scalar;
    }

    template<class U>
    Self& operator/=(U scalar) {
        return *this = *this / scalar;
    }

    ////////////////////////////////////////////////////////////
    // Matrix operations
    ////////////////////////////////////////////////////////////

    template<class U, u32 p_cols>
    friend Result<U, m_rows_, p_cols> operator*(
        const Self&                    lhs,
        const Mat<U, n_cols_, p_cols>& rhs) {
        Result<U, m_rows_, p_cols> res{};

        for ( u32 row = 0; row < m_rows_; ++row ) {
            for ( u32 col = 0; col < p_cols; ++col ) {
                for ( u32 k = 0; k < n_cols_; ++k ) {
                    res(row, col) += lhs(row, k) * rhs(k, col);
                }
            }
        }

        return res;
    }

    Mat<T, n_cols_, m_rows_> transposed() const {
        Matrix<T, n_cols_, m_rows_> res;
        for ( u32 row = 0; row < m_rows_; ++row ) {
            for ( u32 col = 0; col < n_cols_; ++col ) {
                res(col, row) = (*this)(row, col);
            }
        }

        return res;
    }

    ////////////////////////////////////////////////////////////
    // STD operations
    ////////////////////////////////////////////////////////////

    friend Self abs(const Self& mat) {
        return Base::unary_element_wise(mat, [](const T& x) {
            using namespace std;
            return abs(x);
        });
    }

    friend Self round(const Self& mat) {
        return Base::unary_element_wise(mat, [](const T& x) {
            using namespace std;
            return round(x);
        });
    }

    template<class U>
    friend Result<U> min(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, [](const T& x) {
            using namespace std;
            return min(x);
        });
    }

    template<class U>
    friend Result<U> max(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, [](const T& x) {
            using namespace std;
            return max(x);
        });
    }


    ////////////////////////////////////////////////////////////
    // Comparison operators
    ////////////////////////////////////////////////////////////

    template<class U>
    friend Mat<std::strong_ordering> operator<=>(const Self&   lhs,
                                                 const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, std::compare_three_way());
    }

    template<class U>
    friend Bool operator==(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, std::equal_to<>());
    }

    template<class U>
    friend Bool operator!=(const Self& lhs, const Mat<U>& rhs) {
        return !(lhs == rhs);
    }

    template<class U>
    friend Bool operator<(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) == std::strong_ordering::less;
    }

    template<class U>
    friend Bool operator<=(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) != std::strong_ordering::greater;
    }

    template<class U>
    friend Bool operator>(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) == std::strong_ordering::greater;
    }

    template<class U>
    friend Bool operator>=(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) != std::strong_ordering::less;
    }
};


template<u32 m_rows_, u32 n_cols_, template<class, u32, u32> class Derived>
struct MatrixOperations<bool, m_rows_, n_cols_, Derived>
    : public MatrixData<bool, m_rows_, n_cols_, Derived>
{
    using Base = MatrixData<bool, m_rows_, n_cols_, Derived>;
    using Base::Base;

    using Self = Base::Self;
    using Bool = Self;

    ////////////////////////////////////////////////////////////
    // Unary operators
    ////////////////////////////////////////////////////////////

    Bool operator!() const {
        return Base::unary_element_wise(*this, std::logical_not<>());
    }

    ////////////////////////////////////////////////////////////
    // Binary operators
    ////////////////////////////////////////////////////////////

    friend Bool operator&&(const Bool& lhs, const Bool& rhs) {
        return Base::element_wise(lhs, rhs, std::logical_and<>());
    }

    friend Bool operator||(const Bool& lhs, const Bool& rhs) {
        return Base::element_wise(lhs, rhs, std::logical_or<>());
    }

    ////////////////////////////////////////////////////////////
    // Reductions
    ////////////////////////////////////////////////////////////

    bool all() { std::all_of(this->begin(), this->end(), std::identity()); }
    bool any() { std::any_of(this->begin(), this->end(), std::identity()); }
    bool none() { std::none_of(this->begin(), this->end(), std::identity()); }
};


////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \brief Vector class
///
/// Vectors are always column vectors,
///     use transpose(vec) to obtain a row vector when needed.
////////////////////////////////////////////////////////////
template<class T, u32 dim>
using Vector = Matrix<T, dim, 1>;


template<class T,
         u32 m_rows_,
         u32 n_cols_,
         template<class, u32, u32>
         class Derived>
struct VectorOperations : public MatrixData<T, m_rows_, n_cols_, Derived>
{
};

template<u32 m_rows_, u32 n_cols_, template<class, u32, u32> class Derived>
struct VectorOperations<bool, m_rows_, n_cols_, Derived>
    : public MatrixData<bool, m_rows_, n_cols_, Derived>
{
};

template<class T, u32 m_rows_, template<class, u32, u32> class Derived>
struct VectorOperations<T, m_rows_, 1, Derived>
    : public MatrixData<T, m_rows_, 1, Derived>
{
    using Base = MatrixData<T, m_rows_, 1, Derived>;
    using Self = Base::Self;

    template<class U>
    auto dot(const Vector<U, m_rows_>& rhs) const {
        return (this->transposed() * rhs)[0];
    }

    T norm_squared() const { return this->dot(*this); }
    T norm() const { return std::sqrt(norm_squared()); }

    Self normalized() const { return *this / norm(); }

    template<class U>
    auto projected_on(const Vector<U, m_rows_>& that) const {
        return that * this->dot(that) / that.norm_squared();
    }
};

template<class T,
         u32 m_rows_,
         u32 n_cols_,
         template<class, u32, u32>
         class Derived>
struct Vec2Operations : public VectorOperations<T, m_rows_, n_cols_, Derived>
{
};

template<class T, template<class, u32, u32> class Derived>
struct Vec2Operations<T, 2, 1, Derived>
    : public VectorOperations<T, 2, 1, Derived>
{
    constexpr auto perp(bool clockwise) const {
        return clockwise ? -this->perp(false)
                         : Derived{-this->data[1], this->data[0]};
    }

    template<class U>
    constexpr auto cross(const Vector<U, 2>& rhs) const {
        return this->data[0] * rhs[1] - this->data[1] * rhs[0];
    }
};

template<class T,
         u32 m_rows_,
         u32 n_cols_,
         template<class, u32, u32>
         class Derived>
struct Vec3Operations : public VectorOperations<T, m_rows_, n_cols_, Derived>
{
};

template<class T, template<class, u32, u32> class Derived>
struct Vec3Operations<T, 3, 1, Derived>
    : public VectorOperations<T, 3, 1, Derived>
{
    template<class U>
    using Result = Derived<std::common_type_t<T, U>, 3, 1>;

    template<class U>
    constexpr Result<U> cross(const Vector<U, 3>& rhs) const {
        using VecT = Vector<T, 2>;
        using VecU = Vector<U, 2>;
        return Result<U>{
            VecT(this->data[1], this->data[2]).cross(VecU{rhs[1], rhs[2]}),
            -VecT(this->data[0], this->data[2]).cross(VecU{rhs[0], rhs[2]}),
            VecT(this->data[0], this->data[1]).cross(VecU{rhs[0], rhs[1]})};
    }
};


template<class T, u32 m, u32 n, template<class, u32, u32> class Derived>
struct StaticConstructors : public MatrixOperations<T, 1, 1, Derived>
{
};

template<class T, u32 dim, template<class, u32, u32> class Derived>
struct StaticConstructors<T, dim, dim, Derived>
    : public MatrixOperations<T, dim, dim, Derived>
{
    static inline Matrix<T, dim, dim> identity();
    static inline Matrix<T, dim, dim> diagonal(const Vector<T, dim>& elements);
};

template<class T, u32 dim, template<class, u32, u32> class Derived>
struct StaticConstructors<T, dim, 1, Derived>
    : public MatrixOperations<T, dim, 1, Derived>
{
    static inline Vector<T, dim> i();
    static inline Vector<T, dim> j();
    static inline Vector<T, dim> k();
    static inline Vector<T, dim> w();
};

template<class T, template<class, u32, u32> class Derived>
struct StaticConstructors<T, 1, 1, Derived>
    : public MatrixOperations<T, 1, 1, Derived>
{
    explicit operator T() const { return this->data[0]; }
};

////////////////////////////////////////////////////////////
/// \ingroup LinearAlgebra
/// \brief Matrix class for small matrices
///
/// Most operations are provided as global functions
/// \see operations
////////////////////////////////////////////////////////////
template<class T, u32 m, u32 n>
struct Matrix : public StaticConstructors<T, m, n, Matrix>
{
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
};


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
/// proj must project x onto its valid bounds (ie clamp its values). With
/// the following signature:
///     T proj(const Vector<T, n>& x, u32 index)
/// where x[index] must be clamped and returned.
///
/// when the absolute change of components of x drops below epsilon,
/// iteration terminates.
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
/// when the absolute change of components of x drops below epsilon,
/// iteration terminates.
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

}  // namespace sw


namespace std
{

template<class T, class U, sw::u32 m, sw::u32 n>
struct common_type<sw::Matrix<T, m, n>, sw::Matrix<U, m, n>>
{
    typedef sw::Matrix<typename std::common_type<T, U>::type, m, n> type;
};

}  // namespace std
