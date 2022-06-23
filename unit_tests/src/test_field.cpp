#include "hila.h"
#include "catch.hpp"

using MyType = float;

class FieldTest {

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

TEST_CASE_METHOD(FieldTest, "Test field constructor as Nullptr", "[Create]") {
    REQUIRE_FALSE(dummy_field.is_allocated());
}

TEST_CASE_METHOD(FieldTest, "Test field allocation", "[Create]") {
    dummy_field.allocate();
    REQUIRE(dummy_field.is_allocated());
}

TEST_CASE_METHOD(FieldTest, "Test field distructor", "[Create]") {
    dummy_field.~Field();
    REQUIRE_FALSE(dummy_field.is_allocated());
}

TEST_CASE_METHOD(FieldTest, "Test field scalar constructor", "[Create]") {
    fill_dummy_field(1);
    REQUIRE(temporary_field(1) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Test field zero constructor", "[Create]") {
    fill_dummy_field(0);
    REQUIRE(temporary_field(0) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Test field copy constructor", "[Create]") {
    fill_dummy_field();
    REQUIRE(temporary_field(dummy_field) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Test field copy", "[Assignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    REQUIRE(temporary_field == dummy_field);

}

TEST_CASE_METHOD(FieldTest, "Test field assignment from value", "[Assignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field(1);
    temporary_field = 1;
    REQUIRE(temporary_field == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Test field assignment from zero", "[Assignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field(0);
    temporary_field = 0;
    REQUIRE(temporary_field == dummy_field);
}




TEST_CASE_METHOD(FieldTest, "Test field arithmetic with constant", "[Arithmetic]") {
    fill_dummy_field(2);
    REQUIRE((temporary_field(1)+=1) == dummy_field);
    REQUIRE((temporary_field(3)-=1) == dummy_field);
    REQUIRE((temporary_field(1)*=2) == dummy_field);
    REQUIRE((temporary_field(4)/=2) == dummy_field);
    REQUIRE((temporary_field(1)+ 1) == dummy_field);
    REQUIRE((temporary_field(3)- 1) == dummy_field);
    REQUIRE((temporary_field(1)* 2) == dummy_field);
    REQUIRE((temporary_field(4)/ 2) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Test field arithmetic with value", "[Arithmetic]") {
    fill_dummy_field(2);
    REQUIRE((temporary_field(1) += temporary_field(1)) == dummy_field);
    REQUIRE((temporary_field(3) -= temporary_field(1)) == dummy_field);
    REQUIRE((temporary_field(1) *= temporary_field(2)) == dummy_field);
    REQUIRE((temporary_field(4) /= temporary_field(2)) == dummy_field);
    REQUIRE((temporary_field(1) +  temporary_field(1)) == dummy_field);
    REQUIRE((temporary_field(3) -  temporary_field(1)) == dummy_field);
    REQUIRE((temporary_field(1) *  temporary_field(2)) == dummy_field);
    REQUIRE((temporary_field(4) /  temporary_field(2)) == dummy_field);
}

//UNARY OPERATOR?

TEST_CASE_METHOD(FieldTest, "Test field squarenorm", "[Mathematical methods]") {
    fill_dummy_field(1);
    REQUIRE(dummy_field.squarenorm()/lattice->volume() == 1);
}
