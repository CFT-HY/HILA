#include <gtest/gtest.h>
#include "hila.h"

using MyType = float;

struct FieldTest : public testing::Test {
public:
    Field<MyType> dummy_field;

    template <typename T>
    Field<MyType> temporary_field(T assign_val) {
        Field<MyType> temporary_field = assign_val;
        return temporary_field;
    }

    void fill_dummy_field() {
        onsites(ALL) this->dummy_field[X] = hila::random();
    }

    void fill_dummy_field(MyType assign_val) {
        onsites(ALL) this->dummy_field[X] = assign_val;
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
    fill_dummy_field(1);
    EXPECT_EQ(temporary_field(1),dummy_field);
}

TEST_F(FieldTest, ZeroConstructor) {
    fill_dummy_field(0);
    EXPECT_EQ(temporary_field(0), dummy_field);
}

TEST_F(FieldTest, CopyConstructor) {
    fill_dummy_field();
    EXPECT_EQ(temporary_field(dummy_field),dummy_field);
}

TEST_F(FieldTest, AssignmentFromField) {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    EXPECT_EQ(temporary_field,dummy_field);

}

TEST_F(FieldTest, AssignmentFromElement) {
    Field<MyType> temporary_field;
    fill_dummy_field(1);
    temporary_field = 1;
    EXPECT_EQ(temporary_field,dummy_field);
}

TEST_F(FieldTest, AssignmentFromZero) {
    Field<MyType> temporary_field;
    fill_dummy_field(0);
    temporary_field = 0;
    EXPECT_EQ(temporary_field,dummy_field);
}

TEST_F(FieldTest, ArithmeticWithConst) {
    fill_dummy_field(2);
    EXPECT_EQ(temporary_field(1)+=1,dummy_field);
    EXPECT_EQ(temporary_field(3)-=1,dummy_field);
    EXPECT_EQ(temporary_field(1)*=2,dummy_field);
    EXPECT_EQ(temporary_field(4)/=2,dummy_field);
    EXPECT_EQ(temporary_field(1)+1,dummy_field);
    EXPECT_EQ(temporary_field(3)-1,dummy_field);
    EXPECT_EQ(temporary_field(1)*2,dummy_field);
    EXPECT_EQ(temporary_field(4)/2,dummy_field);
}

TEST_F(FieldTest, ArithmeticWithField) {
    fill_dummy_field(2);
    EXPECT_EQ(temporary_field(1)+=temporary_field(1),dummy_field);
    EXPECT_EQ(temporary_field(3)-=temporary_field(1),dummy_field);
    EXPECT_EQ(temporary_field(1)*=temporary_field(2),dummy_field);
    EXPECT_EQ(temporary_field(4)/=temporary_field(2),dummy_field);
    EXPECT_EQ(temporary_field(1)+temporary_field(1),dummy_field);
    EXPECT_EQ(temporary_field(3)-temporary_field(1),dummy_field);
    EXPECT_EQ(temporary_field(1)*temporary_field(2),dummy_field);
    EXPECT_EQ(temporary_field(4)/temporary_field(2),dummy_field);
}

//UNARY OPERATOR?

TEST_F(FieldTest, Squarenorm) {
    fill_dummy_field(1);
    EXPECT_EQ(dummy_field.squarenorm()/lattice->volume(),1);
}