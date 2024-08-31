#include "catch2/catch_test_macros.hpp"

#include "shockwave/memory/allocator.h"
#include <type_traits>

TEST_CASE("memory/allocator") {
    SECTION("rebind") {
        STATIC_REQUIRE(
            std::is_same_v<std::allocator<int>,
                           sw::rebind<std::allocator<float>>::to<int>>);

        STATIC_REQUIRE_FALSE(
            std::is_same_v<std::allocator<int>,
                           sw::rebind<std::allocator<int>>::to<float>>);
    }
}
