#include "catch2/catch_test_macros.hpp"
#include "catch2/catch_template_test_macros.hpp"
#include "shockwave/utility/assert.h"

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

        Index generation_overflow_wraps{0, 1 << Index::gen_bits};
        REQUIRE(first_gen == generation_overflow_wraps);

        REQUIRE_THROWS_AS(
            Index(static_cast<std::size_t>(1) << Index::index_bits, 0),
            sw::Exception);

        REQUIRE_NOTHROW(
            Index((static_cast<std::size_t>(1) << Index::index_bits) - 1, 0));
    }

    SECTION("No Generations Index") {
        using Index = Pool<int, 0>::Index;

        Index first_gen{0, 0};
        Index second_gen{0, 1};
        REQUIRE(first_gen == second_gen);  // wrapped

        Index second_of_first_gen{1, 0};
        REQUIRE(first_gen < second_of_first_gen);
    }
}

TEMPLATE_PRODUCT_TEST_CASE_SIG("memory/Pool",
                               "[template][product][nttp]",
                               ((typename T, size_t S), T, S),
                               Pool,
                               ((int, 8), (int, 0))) {
    SECTION("Pool") {
        TestType pool;
        using Index = TestType::Index;

        SECTION("Bad index returns nullopt") {
            REQUIRE(pool.get(Index{0, 0}) == std::nullopt);
        }

        SECTION("Can take objects and give them back") {
            auto index = pool.take();
            REQUIRE(index == Index{0, 0});

            REQUIRE(pool.get(index).has_value());

            int* value = *pool.get(index);
            *value     = 10;
            REQUIRE(**pool.get(index) == 10);

            auto index2 = pool.take();
            REQUIRE(index2 == Index{1, 0});

            pool.give_back(index);
            REQUIRE(pool.get(index) == std::nullopt);

            index = pool.take();
            if constexpr ( Index::gen_bits == 0 ) {
                REQUIRE(index == Index{0, 0});
            } else {
                REQUIRE(index != Index{0, 0});  // can give {2, 0} or {0, 1}
                REQUIRE(pool.get(Index{0, 0}) == std::nullopt);
            }
        }
    }
}
