#include <gtest/gtest.h>
#include "hila.h"

TEST(LatticeTest, Size) {
    ASSERT_EQ(lattice->volume(), 128*128*128);
}