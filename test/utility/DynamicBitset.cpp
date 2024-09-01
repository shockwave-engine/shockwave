#include "catch2/catch_test_macros.hpp"
#include "shockwave/utility/assert.h"

#include <ranges>

#include "shockwave/utility/DynamicBitset.h"

TEST_CASE("utility/DynamicBitset") {
    sw::DynamicBitset bitset;

    SECTION("Resizes to a byte width") {
        REQUIRE(bitset.size() == 0);

        bitset.resize(1);
        REQUIRE(bitset.size() == 8);

        bitset.resize(8);
        REQUIRE(bitset.size() == 8);

        bitset.resize(9);
        REQUIRE(bitset.size() == 16);
    }

    SECTION("Can set / unset bits and get their value") {
        bitset.resize(8);

        for ( std::size_t bit : std::views::iota(0u, bitset.size()) ) {
            REQUIRE_FALSE(bitset.get(bit));
        }

        for ( std::size_t bit : std::views::iota(0u, bitset.size()) ) {
            bitset.set(bit);
        }

        for ( std::size_t bit : std::views::iota(0u, bitset.size()) ) {
            REQUIRE(bitset.get(bit));
            bitset.unset(bit);
            REQUIRE_FALSE(bitset.get(bit));
        }
    }
}
