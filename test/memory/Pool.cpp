#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_template_test_macros.hpp"

#include "shockwave/memory/Pool.h"
#include <cstddef>

template<class T, std::size_t gen_bits>
using Pool = sw::Pool<T, std::allocator<T>, gen_bits>;

TEST_CASE("memory/Pool") {
    SECTION("Generational Index") {
        using Index = Pool<int, 8>::Index;

        Index first_gen{0, 0};
        Index second_gen{0, 1};
        REQUIRE(first_gen < second_gen);

        REQUIRE(first_gen == Index{0, 0});
        REQUIRE(first_gen != second_gen);

        Index second_of_first_gen{1, 0};
        REQUIRE(first_gen < second_of_first_gen);
        REQUIRE(second_of_first_gen < second_gen);
    }
}

TEMPLATE_PRODUCT_TEST_CASE_SIG("memory/Pool",
                               "[template][product][nttp]",
                               ((typename T, size_t S), T, S),
                               Pool,
                               ((int, 8), (int, 0))) {
    TestType x;
    SECTION("Pool") {}
}
