#include <gtest/gtest.h>
#include "hila.h"

struct FieldTest : public testing::Test {
public:
    Field<float> dummy_field;

    template <typename T>
    Field<float> temporary_field(T assign_val) {
        Field<float> temporary_field = assign_val;
        return temporary_field;
    }

    void fill_dummy_field() {
        onsites(ALL) this->dummy_field[X] = hila::random();
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
    EXPECT_EQ(temporary_field(1),1);
}

TEST_F(FieldTest, ZeroConstructor) {
    EXPECT_EQ(temporary_field(0),0);
}

TEST_F(FieldTest, CopyConstructor) {
    fill_dummy_field();
    EXPECT_EQ(temporary_field(dummy_field),dummy_field);
}