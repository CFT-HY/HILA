#include "hila.h"
#include "catch.hpp"

TEST_CASE("Size", "[default]") {
    REQUIRE(lattice->volume() == 32*32*32);
}