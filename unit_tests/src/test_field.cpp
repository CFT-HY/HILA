#include <gtest/gtest.h>
#include "hila.h"

using MyType = double;

struct FieldTest : public testing::Test {
public:
    Field<double> dummy_field;
};

TEST_F(FieldTest, TestFieldCreationNULL) {
    EXPECT_FALSE(dummy_field.is_allocated());
}

TEST_F(FieldTest, TestFieldAllocation) {
    dummy_field.allocate();
    EXPECT_TRUE(dummy_field.is_allocated());
}

TEST_F(FieldTest, TestFieldCopy) {
    dummy_field = 1;
    auto dummy_field_copy = dummy_field;
    onsites(ALL) {
        EXPECT_EQ(dummy_field_copy[X],dummy_field[X]);
    }
}