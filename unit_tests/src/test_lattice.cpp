#include "hila.h"
#include "catch.hpp"
#include "catch_main.hpp"

TEST_CASE("Size", "[MPI][.]") {
    REQUIRE(lattice->volume() == lattice_size*lattice_size*lattice_size);
}

TEST_CASE("Node information", "[MPI][.]") {
    if (hila::myrank)
    REQUIRE(true);
}