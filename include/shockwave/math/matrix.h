#pragma once

#include <algorithm>
#include <compare>
#include <functional>
#include <initializer_list>
#include <ostream>
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

private:

    template<std::invocable<std::size_t> Generator, std::size_t... Is>
    constexpr MatrixData(Generator&& generate, std::index_sequence<Is...>)
        : data{generate(Is)...} {}

public:

    constexpr MatrixData() = default;

    template<std::invocable<u32> Fn>
    constexpr MatrixData(Fn&& generator)
        : MatrixData(std::forward<Fn>(generator),
                     std::make_index_sequence<size()>()) {}

    ////////////////////////////////////////////////////////////
    // Indexing
    ////////////////////////////////////////////////////////////

    constexpr static u32 size() { return m_rows * n_cols; }

    constexpr static u32 to_index(u32 row, u32 col) {
        return row * n_cols + col;
    }

    constexpr u32 static to_row(u32 index) { return index / n_cols_; }
    constexpr u32 static to_col(u32 index) { return index % n_cols_; }

    T& operator()(u32 row, u32 col) { return this->data[to_index(row, col)]; }
    const T& operator()(u32 row, u32 col) const {
        return this->data[to_index(row, col)];
    }

    T&       operator[](u32 index) { return this->data[index]; }
    const T& operator[](u32 index) const { return this->data[index]; }

    iterator begin() { return this->data; }
    iterator end() { return this->data + size(); }

    const_iterator begin() const { return this->data; }
    const_iterator end() const { return this->data + size(); }

    ////////////////////////////////////////////////////////////
    // Operations
    ////////////////////////////////////////////////////////////

    template<std::invocable<const T&> Fn>
    static Mat<std::result_of_t<Fn(const T&)>> unary_element_wise(
        const Self& mat,
        Fn&&        fn) {
        return Mat<std::result_of_t<Fn(const T&)>>(
            [&](u32 i) { return fn(mat[i]); });
    }

    template<class U, std::invocable<const T&, const U&> Fn>
    static Mat<std::result_of_t<Fn(const T&, const U&)>>
    element_wise(const Self& lhs, const Mat<U>& rhs, Fn&& fn) {
        return Mat<std::result_of_t<Fn(const T&, const U&)>>(
            [&](u32 i) { return fn(lhs[i], rhs[i]); });
    }

    ////////////////////////////////////////////////////////////
    // Debug
    ////////////////////////////////////////////////////////////

    friend std::ostream& operator<<(std::ostream& os, const Self& mat) {
        os << '{';
        for ( u32 i = 0; i < m_rows_; i++ ) {
            os << '{';
            for ( u32 j = 0; j < n_cols_; j++ ) {
                os << mat(i, j) << (j == n_cols_ - 1 ? "}" : ", ");
            }
        }

        return os << '}';
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
    using Result = Base::template Result<U, m, n>;

    template<class U = T, u32 m = Base::m_rows, u32 n = Base::n_cols>
    using Mat = Base::template Mat<U, m, n>;

    using Self = Mat<>;

    using Bool = Mat<bool>;

    template<class U>
    using Cmp = Mat<std::compare_three_way_result_t<T, U>>;

    ////////////////////////////////////////////////////////////
    // Unary operators
    ////////////////////////////////////////////////////////////

    Result<> operator+() const { return static_cast<const Self&>(*this); }
    Result<> operator-() const {
        return Base::unary_element_wise(static_cast<const Self&>(*this),
                                        std::negate<>());
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

    template<class U, std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
    friend Result<U> operator*(U scalar, const Self& mat) {
        return Base::unary_element_wise(mat,
                                        [=](const T& x) { return x * scalar; });
    }

    template<class U, std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
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
        return static_cast<Self&>(*this
                                  = static_cast<const Self&>(*this) + other);
    }

    template<class U>
    Self& operator-=(const Mat<U>& other) {
        return static_cast<Self&>(*this
                                  = static_cast<const Self&>(*this) - other);
    }

    template<class U>
    Self& operator*=(U scalar) {
        return static_cast<Self&>(*this
                                  = static_cast<const Self&>(*this) * scalar);
    }

    template<class U>
    Self& operator/=(U scalar) {
        return static_cast<Self&>(*this
                                  = static_cast<const Self&>(*this) / scalar);
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
    friend auto min(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, [](const T& x) {
            using namespace std;
            return min(x);
        });
    }

    template<class U>
    friend auto max(const Self& lhs, const Mat<U>& rhs) {
        return Base::element_wise(lhs, rhs, [](const T& x) {
            using namespace std;
            return max(x);
        });
    }


    ////////////////////////////////////////////////////////////
    // Comparison operators
    ////////////////////////////////////////////////////////////

    template<class U>
    friend Cmp<U> operator<=>(const Self& lhs, const Mat<U>& rhs) {
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
        return lhs <=> rhs == Cmp<U>::filled(std::strong_ordering::less);
    }

    template<class U>
    friend Bool operator<=(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) != Cmp<U>::filled(std::strong_ordering::greater);
    }

    template<class U>
    friend Bool operator>(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) == Cmp<U>::filled(std::strong_ordering::greater);
    }

    template<class U>
    friend Bool operator>=(const Self& lhs, const Mat<U>& rhs) {
        return (lhs <=> rhs) != Cmp<U>::filled(std::strong_ordering::less);
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
        return Base::unary_element_wise(static_cast<const Self&>(*this),
                                        std::logical_not<>());
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

    bool all() {
        return std::all_of(this->begin(), this->end(), std::identity());
    }

    bool any() {
        return std::any_of(this->begin(), this->end(), std::identity());
    }

    bool none() {
        return std::none_of(this->begin(), this->end(), std::identity());
    }
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
         class Derived,
         class = void>
struct VectorOperations : public MatrixOperations<T, m_rows_, n_cols_, Derived>
{
    using Base = MatrixOperations<T, m_rows_, n_cols_, Derived>;
    using Base::Base;
};

template<u32 m_rows_, u32 n_cols_, template<class, u32, u32> class Derived>
struct VectorOperations<bool, m_rows_, n_cols_, Derived, void>
    : public MatrixOperations<bool, m_rows_, n_cols_, Derived>
{
    using Base = MatrixOperations<bool, m_rows_, n_cols_, Derived>;
    using Base::Base;
};

template<class T, u32 m_rows_, template<class, u32, u32> class Derived>
struct VectorOperations<T,
                        m_rows_,
                        1,
                        Derived,
                        std::enable_if_t<!std::is_same_v<bool, T>>>
    : public MatrixOperations<T, m_rows_, 1, Derived>
{
    using Base = MatrixOperations<T, m_rows_, 1, Derived>;
    using Base::Base;

    using Self = Base::Self;

    template<class U>
    auto dot(const Vector<U, m_rows_>& rhs) const {
        return (this->transposed() * rhs)[0];
    }

    T norm_squared() const { return this->dot(*this); }
    T norm() const { return std::sqrt(norm_squared()); }

    Self normalized() const { return static_cast<const Self&>(*this) / norm(); }

    template<class U>
    auto projected_on(const Vector<U, m_rows_>& that) const {
        return that * this->dot(that) / that.norm_squared();
    }

    template<class U>
    auto element_wise_mul(const Vector<U, m_rows_>& rhs) const {
        return Base::element_wise(*this, rhs, std::multiplies<>());
    }
};

template<class T,
         u32 m_rows_,
         u32 n_cols_,
         template<class, u32, u32>
         class Derived>
struct Vec2Operations : public VectorOperations<T, m_rows_, n_cols_, Derived>
{
    using Base = VectorOperations<T, m_rows_, n_cols_, Derived>;
    using Base::Base;
};

template<class T, template<class, u32, u32> class Derived>
struct Vec2Operations<T, 2, 1, Derived>
    : public VectorOperations<T, 2, 1, Derived>
{
    using Base = VectorOperations<T, 2, 1, Derived>;
    using Base::Base;

    constexpr auto perp(bool clockwise = false) const {
        return clockwise ? Derived<T, 2, 1>{this->data[1], -this->data[0]}
                         : Derived<T, 2, 1>{-this->data[1], this->data[0]};
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
struct Vec3Operations : public Vec2Operations<T, m_rows_, n_cols_, Derived>
{
    using Base = Vec2Operations<T, m_rows_, n_cols_, Derived>;
    using Base::Base;
};

template<class T, template<class, u32, u32> class Derived>
struct Vec3Operations<T, 3, 1, Derived>
    : public Vec2Operations<T, 3, 1, Derived>
{
    using Base = Vec2Operations<T, 3, 1, Derived>;
    using Base::Base;

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
struct StaticConstructors : public Vec3Operations<T, m, n, Derived>
{
    using Base = Vec3Operations<T, m, n, Derived>;
    using Base::Base;
};

template<class T, u32 dim, template<class, u32, u32> class Derived>
struct StaticConstructors<T, dim, dim, Derived>
    : public Vec3Operations<T, dim, dim, Derived>
{
    using Base = Vec3Operations<T, dim, dim, Derived>;
    using Base::Base;

    using Mat = Derived<T, dim, dim>;

    static inline Mat identity() {
        return diagonal(Vector<T, dim>::filled(T{1}));
    }

    static inline Mat diagonal(const Vector<T, dim>& elements) {
        return Mat([&](u32 i) constexpr {
            return Base::to_row(i) == Base::to_col(i)
                     ? elements[Base::to_row(i)]
                     : T{};
        });
    }
};

template<class T, u32 dim, template<class, u32, u32> class Derived>
struct StaticConstructors<T, dim, 1, Derived>
    : public Vec3Operations<T, dim, 1, Derived>
{
    using Base = Vec3Operations<T, dim, 1, Derived>;
    using Base::Base;

    using Vec = Derived<T, dim, 1>;

    static inline Vec i() {
        return Vec([](u32 i) { return i == 0 ? T{1} : T{}; });
    }
    static inline Vec j() {
        return Vec([](u32 i) { return i == 1 ? T{1} : T{}; });
    }
    static inline Vec k() {
        return Vec([](u32 i) { return i == 2 ? T{1} : T{}; });
    }
    static inline Vec w() {
        return Vec([](u32 i) { return i == 3 ? T{1} : T{}; });
    }
};

template<class T, template<class, u32, u32> class Derived>
struct StaticConstructors<T, 1, 1, Derived>
    : public Vec3Operations<T, 1, 1, Derived>
{
    using Base = Vec3Operations<T, 1, 1, Derived>;
    using Base::Base;

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
    using Base = StaticConstructors<T, m, n, Matrix>;

    ////////////////////////////////////////////////////////////
    // Constructors
    ////////////////////////////////////////////////////////////

    constexpr Matrix() = default;

    template<std::invocable<u32> Fn>
    constexpr Matrix(Fn&& generator) : Base(std::forward<Fn>(generator)) {}

    template<class... Args,
             std::enable_if_t<(sizeof...(Args) > 1), bool> = true>
    constexpr Matrix(Args&&... values)
        : Matrix({static_cast<T>(std::forward<Args>(values))...}) {}

    constexpr Matrix(std::initializer_list<T> values) {
        for ( u32 i = 0; i < this->size(); i++ ) {
            this->data[i] = i < values.size() ? *(values.begin() + i) : T{};
        }
    }

    template<class U>
    explicit constexpr Matrix(const Matrix<U, m, n>& other) {
        for ( u32 i = 0; i < this->size(); i++ ) {
            this->data[i] = static_cast<T>(other.data[i]);
        }
    }

    ////////////////////////////////////////////////////////////
    // Static Constructors
    ////////////////////////////////////////////////////////////

    static inline Matrix filled(T val) {
        return Matrix([=](auto) { return val; });
    }

    template<class... Args, std::enable_if_t<sizeof...(Args) == m, bool> = true>
    static inline Matrix from_rows(Args&&... rows) {
        return from_rows({std::forward<Args>(rows)...});
    }

    template<class U>
    static inline Matrix from_rows(
        const std::initializer_list<Vector<U, n>>& rows) {
        return from_rows(Vector<Vector<U, n>, m>(rows));
    }

    template<class U>
    static inline Matrix from_rows(const Vector<Vector<U, n>, m>& rows) {
        Matrix mat{};

        auto it = mat.begin();
        for ( const auto& row : rows )
            for ( const auto& val : row )
                *it++ = val;

        return mat;
    }

    template<class... Args, std::enable_if_t<sizeof...(Args) == n, bool> = true>
    static inline Matrix from_cols(Args&&... cols) {
        return from_cols({std::forward<Args>(cols)...});
    }

    template<class U>
    static inline Matrix from_cols(
        const std::initializer_list<Vector<U, m>>& cols) {
        return Matrix<T, n, m>::from_rows(cols).transposed();
    }

    template<class U>
    static inline Matrix from_cols(const Vector<Vector<U, m>, n>& cols) {
        return Matrix<T, n, m>::from_rows(cols).transposed();
    }

    ////////////////////////////////////////////////////////////
    // As row / column vectors
    ////////////////////////////////////////////////////////////

    inline Vector<Vector<T, n>, m> as_rows() const {
        Vector<Vector<T, n>, m> rows;

        for ( u32 i = 0; i < m; i++ ) {
            for ( u32 j = 0; j < n; j++ ) {
                rows[i][j] = (*this)(i, j);
            }
        }

        return rows;
    }

    inline Vector<Vector<T, m>, n> as_cols() const {
        return this->transposed().as_rows();
    }
};


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
