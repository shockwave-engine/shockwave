#pragma once

#include "shockwave/math/matrix.h"

namespace sw
{

template<class T, u32 m, u32 n>
template<class U>
Matrix<T, m, n> Matrix<T, m, n>::fromRows(
    const std::initializer_list<Vector<U, n>>& rows) {
    return fromRows(Vector<Vector<U, n>, m>{rows});
}

template<class T, u32 m, u32 n>
template<class U>
Matrix<T, m, n> Matrix<T, m, n>::fromRows(const Vector<Vector<U, n>, m>& rows) {
    Matrix mat{};

    auto it = mat.begin();
    for ( const auto& row : rows )
        for ( const auto& val : row )
            *it++ = val;

    return mat;
}

template<class T, u32 m, u32 n>
template<class U>
Matrix<T, m, n> Matrix<T, m, n>::fromCols(
    const std::initializer_list<Vector<U, m>>& cols) {
    return transpose(Matrix<T, n, m>::fromRows(cols));
}

template<class T, u32 m, u32 n>
template<class U>
Matrix<T, m, n> Matrix<T, m, n>::fromCols(const Vector<Vector<U, m>, n>& cols) {
    return transpose(Matrix<T, n, m>::fromRows(cols));
}

template<class T, u32 m, u32 n>
Matrix<T, m, n> Matrix<T, m, n>::filled(T val) {
    Matrix<T, m, n> mat{};
    for ( auto& elem : mat )
        elem = val;

    return mat;
}

template<class T, u32 m, u32 n>
Vector<Vector<T, n>, m> Matrix<T, m, n>::asRows() const {
    Vector<Vector<T, n>, m> rows{};

    auto it = this->begin();
    for ( Vector<T, n>& row : rows )
        for ( T& val : row )
            val = *it++;

    return rows;
}

template<class T, u32 m, u32 n>
Vector<Vector<T, m>, n> Matrix<T, m, n>::asCols() const {
    return transpose(*this).asRows();
}


template<class T>
class Solver<T, 2>
{
public:

    Solver(const Matrix<T, 2, 2>& A) : A_{A} {
        invDet_  = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
        isValid_ = (invDet_ != 0.f);
        if ( isValid_ )
            invDet_ = 1.f / invDet_;
    }

    template<class U>
    Vector<Promoted<T, U>, 2> solve(const Vector<U, 2>& b) const {
        SIMU_ASSERT(isValid_,
                    "Solver invalid, ensure the original matrix has full rank");

        // use cramer's rule
        return invDet_
             * Vector<Promoted<T, U>, 2>{b[0] * A_(1, 1) - A_(0, 1) * b[1],
                                         A_(0, 0) * b[1] - b[0] * A_(1, 0)};
    }

    Matrix<T, 2, 2> original() const { return A_; }

    bool isValid() const { return isValid_; }

private:

    template<class U>
    friend class LcpSolver;

    Matrix<T, 2, 2> A_;
    T               invDet_;
    bool            isValid_;
};

template<class T, class U, u32 n>
Vector<Promoted<T, U>, n> solve(const Matrix<T, n, n>& A,
                                const Vector<U, n>&    b) {
    return Solver{A}.solve(b);
}

template<class T, u32 n>
Solver<T, n>::Solver(const Matrix<T, n, n>& A) : R_{} {
    // modified Gram-Schmidt for QR decomposition
    // https://www.math.uci.edu/~ttrogdon/105A/html/Lecture23.html

    auto Q = A.asCols();

    for ( u32 col = 0; col < n; ++col ) {
        R_(col, col) = norm(Q[col]);

        isValid_ = isValid_ && R_(col, col) != 0.f;
        if ( !isValid_ )
            return;

        Q[col] /= R_(col, col);

        for ( u32 nextCol = col + 1; nextCol < n; ++nextCol ) {
            R_(col, nextCol)  = dot(Q[col], Q[nextCol]);
            Q[nextCol]       -= R_(col, nextCol) * Q[col];
        }
    }

    QT_ = transpose(Matrix<T, n, n>::fromCols(Q));
}

template<class T, u32 n>
template<class U>
Vector<Promoted<T, U>, n> Solver<T, n>::solve(const Vector<U, n>& b) const {
    SIMU_ASSERT(isValid_,
                "Solver invalid, ensure the original matrix has full rank");

    Vector<Promoted<T, U>, n> c = QT_ * b;
    Vector<Promoted<T, U>, n> x{};

    for ( u32 row = n; row > 0; --row ) {
        for ( u32 col = row; col < n; ++col ) {
            c[row - 1] -= R_(row - 1, col) * x[col];
        }

        x[row - 1] = c[row - 1] / R_(row - 1, row - 1);
    }

    return x;
}


template<class T, u32 n>
Matrix<T, n, n> invert(const Matrix<T, n, n>& mat) {
    Solver<T, n> solver{mat};

    auto inverse = Matrix<T, n, n>::identity().asCols();
    for ( auto& col : inverse )
        col = solver.solve(col);

    return Matrix<T, n, n>::fromCols(inverse);
}


template<class T, u32 n, std::invocable<Vector<T, n>, u32> Proj>
Vector<T, n> solveInequalities(const Matrix<T, n, n>& A,
                               Vector<T, n>           b,
                               Proj                   proj,
                               Vector<T, n>           initialGuess,
                               float                  epsilon) {
    b = -b;  // Ax - b >= 0

    Vector<T, n> x = initialGuess;

    float eps = 1.f + epsilon;
    while ( eps > epsilon ) {
        eps = 0.f;
        for ( u32 row = 0; row < n; ++row ) {
            float delX = 0.f;

            for ( u32 col = 0; col < row; ++col )
                delX += A(row, col) * x[col];

            for ( u32 col = row + 1; col < n; ++col )
                delX += A(row, col) * x[col];

            delX = -(delX + b[row]) / A(row, row);

            Vector<T, n> newX{x};
            newX[row] = delX;
            delX      = proj(newX, row);

            eps    = std::max(eps, std::abs(delX - x[row]));
            x[row] = delX;
        }
    }

    return x;
}


template<class T, u32 n, std::invocable<Vector<T, n>> Proj>
Vector<T, n> solveInequalities(const Matrix<T, n, n>& A,
                               Vector<T, n>           b,
                               Proj                   proj,
                               Vector<T, n>           initialGuess,
                               float                  epsilon) {
    return solveInequalities(
        A,
        b,
        [=](const Vector<T, n>& x, u32 i) { return proj(x)[i]; },
        initialGuess,
        epsilon);
}


template<class T>
class LcpSolver
{
public:

    LcpSolver(const Matrix<T, 2, 2>& A)
        : solver{A},
          invA11{1.f / A(0, 0)},
          invA22{1.f / A(1, 1)} {
        solver.isValid_ = solver.isValid() && A(0, 0) != 0.f && A(1, 1) != 0.f;
    }

    template<class U>
    Vector<Promoted<T, U>, 2> solve(const Vector<U, 2>& b) const {
        Vector<T, 2> x = solver.solve(b);
        if ( all(x >= Vector<T, 2>::filled(0.f)) )
            return x;

        x = Vector<T, 2>{b[0] * invA11, 0.f};
        if ( (x[0] >= 0.f) && (x[0] * solver.A_(1, 0) - b[1] >= 0) )
            return x;

        x = Vector<T, 2>{0.f, b[1] * invA22};
        if ( (x[1] >= 0.f) && (x[1] * solver.A_(0, 1) - b[0] >= 0) )
            return x;

        x = Vector<T, 2>{};
        if ( all(-b >= Vector<T, 2>::filled(0.f)) )
            return x;

        return x;  // all null
    }

    Matrix<T, 2, 2> original() const { return solver.original(); }
    bool            isValid() const { return solver.isValid(); }

private:

    Solver<T, 2> solver;
    float        invA11;
    float        invA22;
};

template<class T>
Vector<T, 2> solveLcp(const Matrix<T, 2, 2>& A, const Vector<T, 2>& b) {
    return LcpSolver{A}.solve(b);
}

}  // namespace sw
