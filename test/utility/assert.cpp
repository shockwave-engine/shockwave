#include "catch2/catch_test_macros.hpp"

#include "shockwave/utility/assert.h"

TEST_CASE("utility/assert") {
    SECTION("ensure") {
        REQUIRE_NOTHROW([]() { ensure(true, ""); }());

        REQUIRE_THROWS_AS([]() { ensure(false, ""); }(), sw::Exception);
    }
}
