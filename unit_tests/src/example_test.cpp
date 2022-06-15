#include <gtest/gtest.h>
#include "hila.h"
static_assert(NDIM == 3, "NDIM must be 3 here");

TEST(LatticeTest, Size) {
    ASSERT_EQ(lattice->volume(), 256*256*256);
}

TEST(LatticeTest, SizeFail) {
    ASSERT_EQ(lattice->volume(), 255*256*256);
}
 
int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    hila::initialize(argc, argv);
    lattice->setup({256, 256, 256});
    return RUN_ALL_TESTS();
}