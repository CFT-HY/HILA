#include <gtest/gtest.h>
#include "hila.h"

static_assert(NDIM == 3, "NDIM must be 3 here");

int main(int argc, char **argv) {
    int result = 0;

    testing::InitGoogleTest(&argc, argv);
    hila::initialize(argc, argv);
    lattice->setup({32,32,32});
    result = RUN_ALL_TESTS();
    hila::finishrun();
    return result;
}