#include "catch2/catch_test_macros.hpp"

#include "shockwave/shockwave.h"

TEST_CASE("Assertion") {
    shock();
    bla();
    REQUIRE(1 == 1);
}
