#include <gtest/gtest.h>
#include "hila.h"

TEST(LatticeTest, Size) {
    ASSERT_EQ(lattice->volume(), 32*32*32);
}