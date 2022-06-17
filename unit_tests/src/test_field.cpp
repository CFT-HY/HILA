#include <gtest/gtest.h>
#include "hila.h"

struct FieldTest : public testing::Test {
public:
    Field<float> dummy_field;

    template <typename T>
    Field<float> temporary_field(T assign_val) {
        Field<float> temporary_field = assign_val;
        std::cout << &temporary_field << '\n';
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

TEST_F(FieldTest, Destructor) {
    dummy_field.~Field();
    EXPECT_FALSE(dummy_field.is_allocated());
}

TEST_F(FieldTest, ScalarConstructor) {
    Field<float> temp = temporary_field(1);
    std::cout << &temp << '\n';
    onsites(ALL) {
        EXPECT_EQ(temporary_field(1)[X],1.0);
    }
}

TEST_F(FieldTest, ZeroConstructor) {
    Field<float> temp = temporary_field(0);
    std::cout << &temp << '\n';
    onsites(ALL) {
        EXPECT_EQ(temp[X],0);
    }
}

TEST_F(FieldTest, CopyConstructor) {
    dummy_field = 1;
    onsites(ALL) dummy_field = 1;
    EXPECT_TRUE(dummy_field.is_allocated());
    Field<float> temp = dummy_field;
    std::cout << &temp << '\n';
    onsites(ALL) {
        EXPECT_EQ(temp[X],dummy_field[X]);
    }
}