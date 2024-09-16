#include <iostream>
#include <numbers>

#include "catch2/catch_test_macros.hpp"

#include "shockwave/math/matrix.h"

using namespace sw;

// TODO: Everything here should be constexpr


TEST_CASE("Matrix") {
    SECTION("Matrix Data") {
        Mat2 null{};
        REQUIRE((null == Mat2{0, 0, 0, 0}).all());
    }

    SECTION("Construct from vectors") {
        Matrix<int, 2, 4> ref{1, 2, 3, 4, 5, 6, 7, 8};

        REQUIRE((ref
                 == Matrix<int, 2, 4>::from_rows(Vec4i{1, 2, 3, 4},
                                                 Vec4i{5, 6, 7, 8}))
                    .all());

        REQUIRE((ref
                 == Matrix<int, 2, 4>::from_cols(Vec2i{1, 5},
                                                 Vec2i{2, 6},
                                                 Vec2i{3, 7},
                                                 Vec2i{4, 8}))
                    .all());
    }

    SECTION("Filled") {
        Mat2i ref{1, 1, 1, 1};
        REQUIRE((ref == Mat2i::filled(1)).all());
    }

    SECTION("Special constructors") {
        // clang-format off
        REQUIRE((Mat4::identity()
                    == Mat4(1, 0, 0, 0, 
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1)).all());
        // clang-format on

        REQUIRE(
            (Mat4::identity()
             == Mat4::from_cols({Vec4::i(), Vec4::j(), Vec4::k(), Vec4::w()}))
                .all());

        REQUIRE(
            (Mat4::identity()
             == Mat4::from_rows({Vec4::i(), Vec4::j(), Vec4::k(), Vec4::w()}))
                .all());

        // may also be made to not compile
        REQUIRE((Vec3i::w() == Vec3i{}).all());
    }

    SECTION("Unary operators") {
        Mat2i ident = Mat2i::identity();
        REQUIRE((+ident == ident).all());
        REQUIRE((-ident + ident == Mat2i{}).all());
        REQUIRE((+-+-ident == ident).all());
    }

    SECTION("Arithmetic assignement") {
        Mat2 ones = Mat2::filled(1.f);

        ones *= 2;
        REQUIRE((ones == Mat2::filled(2.f)).all());

        ones /= 2;
        REQUIRE((ones == Mat2::filled(1.f)).all());

        ones += ones;
        REQUIRE((ones == Mat2::filled(2.f)).all());

        ones -= ones;
        REQUIRE((ones == Mat2{}).all());
    }

    SECTION("Linear combination") {
        // Binary operators returning by value must promote the return type
        auto combination = Vec3i::i() * 0.5f - Vec3i::i() + 0.5f * Vec3i::j()
                         + Vec3i::k() / 2.f;

        REQUIRE((combination == Vec3{-0.5f, 0.5f, 0.5f}).all());
    }

    SECTION("Matrix multiplication") {
        Mat2 ident = Mat2::identity();
        REQUIRE((ident * ident == ident).all());

        REQUIRE(
            (Vec2::i().transposed() * Vec2::j() == Vector<int, 1>{0}).all());

        REQUIRE(
            (Vec2{1, 2} * Vec2{2, 1}.transposed() == Mat2{2, 1, 4, 2}).all());

        float theta = std::numbers::pi_v<float> / 4.f;
        Mat2  rot{std::cos(theta),
                 -std::sin(theta),
                 std::sin(theta),
                 std::cos(theta)};

        /*REQUIRE(all(approx(Vec2{1, 1}.normalized(), Vec2::filled(1e-6f))*/
        /*                .contains(rot * Vec2::i())));*/
    }

    SECTION("Solver") {
        /*REQUIRE((solve(Mat3::identity(), Vec3::i()) == Vec3::i()).all());*/
        /**/
        /*float theta = std::numbers::pi_v<float> / 4.f;*/
        /*Mat2  rot{std::cos(theta),*/
        /*         -std::sin(theta),*/
        /*         std::sin(theta),*/
        /*         std::cos(theta)};*/
        /**/
        /*REQUIRE(all(approx(invert(rot) * rot, Mat2::filled(EPSILON))*/
        /*                .contains(Mat2::identity())));*/
        /**/
        /*REQUIRE(all(approx(solve(rot * rot, Vec2::j()),
         * Vec2::filled(EPSILON))*/
        /*                .contains(Vec2::i())));*/
    }

    SECTION("Solver (encountered issues)") {
        /*// clang-format off*/
        /*Mat3 A{*/
        /*    2.5f,  -0.5f,  1.5f,*/
        /*    -0.5f, 2.5f,  -1.5f,*/
        /*    1.5f,  -1.5f,  2.5f,*/
        /*};*/
        /*// clang-format on*/
        /**/
        /*Vec3 b{0.2f, 0.2f, 0.f};*/
        /**/
        /*Solver<float, 3> solver{A};*/
        /*Mat3             original = solver.original();*/
        /*REQUIRE(approx(original, Mat3::filled(EPSILON)).contains(A)).all());*/
        /**/
        /*Vec3 x = solver.solve(b);*/
        /**/
        /*REQUIRE(all(*/
        /*    approx(x, Vec3::filled(EPSILON)).contains(Vec3{0.1f, 0.1f,
         * 0.f})));*/
    }

    SECTION("Inequality Solver") {
        SECTION("Can solve equalities") {
            /*// clang-format off*/
            /*Mat3 A{*/
            /*    2.5f,  -0.5f,  1.5f,*/
            /*    -0.5f, 2.5f,  -1.5f,*/
            /*    1.5f,  -1.5f,  2.5f,*/
            /*};*/
            /*// clang-format on*/
            /**/
            /*Vec3 b{0.2f, 0.2f, 0.f};*/
            /**/
            /*Vec3 x = solveInequalities(A, b, [](Vec3 x) { return x; });*/
            /**/
            /*REQUIRE(all(approx(x, Vec3::filled(2 * simu::EPSILON))*/
            /*                .contains(Vec3{0.1f, 0.1f, 0.f})));*/
        }

        SECTION("Can solve LCP") {
            /*// clang-format off*/
            /*Mat2 A{*/
            /*    2, 1, */
            /*    1, 2    */
            /*};*/
            /*// clang-format on*/
            /**/
            /*Vec2 b{5, 6};*/
            /**/
            /*Vec2 x = solveInequalities(A, b, [](Vec2 x) {*/
            /*    return std::max(x, Vec2::filled(0.f));*/
            /*});*/
            /**/
            /*REQUIRE(x >= Vec2::filled(0.f)).all());*/
            /**/
            /*Vec2 res = A * x;*/
            /*REQUIRE(res >= b).all());*/
        }

        SECTION("Can solve MLCP") {
            /*// clang-format off*/
            /*Mat2 A{*/
            /*    2, 1, */
            /*    1, 2    */
            /*};*/
            /*// clang-format on*/
            /**/
            /*Vec2 b{-5, -6};*/
            /**/
            /*Vec2 x = solveInequalities(A, b, [](Vec2 x) {*/
            /*    return Vec2{std::max(x[0], 1.f), std::min(x[1], -1.f)};*/
            /*});*/
            /**/
            /*REQUIRE(x[0] >= 1.f);*/
            /*REQUIRE(x[1] <= -1.f);*/
            /**/
            /*Vec2 res = A * x;*/
            /*REQUIRE(res >= b).all());*/
        }
    }

    SECTION("Vector operations") {
        REQUIRE((Vec3::i().cross(Vec3::j()) == Vec3::k()).all());
        REQUIRE((Vec3::j().cross(Vec3::i()) == -Vec3::k()).all());

        REQUIRE((Vec3::j().cross(Vec3::k()) == Vec3::i()).all());
        REQUIRE((Vec3::k().cross(Vec3::j()) == -Vec3::i()).all());

        REQUIRE((Vec3::k().cross(Vec3::i()) == Vec3::j()).all());
        REQUIRE((Vec3::i().cross(Vec3::k()) == -Vec3::j()).all());

        REQUIRE((Vec3::i().cross(Vec3::i()) == Vec3{}).all());

        REQUIRE((Vec2::i().perp() == Vec2::j()).all());
        REQUIRE((Vec2::i().perp(true) == -Vec2::j()).all());
    }

    SECTION("std and comparison matrices") {
        using namespace std;  // required for ADL

        Vec3 v = Vec3::filled(-1.25f);

        REQUIRE((v == abs(v)).none());

        Vec3 absV = abs(v);
        REQUIRE((v < absV).all());
        REQUIRE((v <= absV).all());

        REQUIRE((absV > round(absV)).all());
        REQUIRE((absV >= round(absV)).all());
    }
}
