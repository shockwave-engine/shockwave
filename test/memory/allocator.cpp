#include "catch2/catch_test_macros.hpp"

#include "shockwave/memory/allocator.h"
#include <type_traits>

TEST_CASE("memory/allocator") {
    SECTION("rebind") {
        STATIC_REQUIRE(
            std::is_same_v<std::allocator<int>,
                           sw::rebind_to<std::allocator<float>, int>>);

        STATIC_REQUIRE_FALSE(
            std::is_same_v<std::allocator<int>,
                           sw::rebind_to<std::allocator<int>, float>>);
    }
}
