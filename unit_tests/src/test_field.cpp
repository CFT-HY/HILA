#include <gtest/gtest.h>
#include "hila.h"

using MyType = double;

struct FieldTest : public testing::Test {
public:
    Field<float> dummy_field;

    template <typename T>
    Field<float> temporary_field(T assign_val) {
        Field<float> temporary_field = assign_val;
        return temporary_field;
    }
};

TEST_F(FieldTest, NullptrConstructor) {
    EXPECT_FALSE(dummy_field.is_allocated());
}

TEST_F(FieldTest, Allocation) {
    dummy_field.allocate();
    EXPECT_TRUE(dummy_field.is_allocated());
}

TEST_F(FieldTest, ScalarConstructor) {
    onsites(ALL) {
        EXPECT_EQ(temporary_field(1)[X],1.0);
    }
}

TEST_F(FieldTest, ZeroConstructor) {
    onsites(ALL) {
        EXPECT_EQ(temporary_field(0)[X],0);
    }
}

TEST_F(FieldTest, CopyConstructor) {
    onsites(ALL) dummy_field[X] = hila::gaussrand();
    onsites(ALL) {
        EXPECT_EQ(temporary_field(dummy_field)[X],dummy_field[X]);
    }
}

