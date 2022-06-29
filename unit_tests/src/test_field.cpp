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

TEST_CASE_METHOD(FieldTest, "Field nullptr", "[FieldCreation]") {
    GIVEN("A defined field without constructor") {
        THEN("Field will be unallocated as nullptr") {
            REQUIRE_FALSE(dummy_field.is_allocated());

        }
    }
}

TEST_CASE_METHOD(FieldTest, "Field allocation", "[FieldCreation]") {
    dummy_field.allocate();
    REQUIRE(dummy_field.is_allocated());
}

TEST_CASE_METHOD(FieldTest, "Field destructor", "[FieldCreation]") {
    dummy_field.~Field();
    REQUIRE_FALSE(dummy_field.is_allocated());
}

TEST_CASE_METHOD(FieldTest, "Field scalar constructor", "[FieldCreation]") {
    fill_dummy_field(1);
    REQUIRE(temporary_field(1) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field zero constructor", "[FieldCreation]") {
    fill_dummy_field(0);
    REQUIRE(temporary_field(0) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field copy constructor", "[FieldCreation]") {
    fill_dummy_field();
    REQUIRE(temporary_field(dummy_field) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field assignment from field", "[FieldAssignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    REQUIRE(temporary_field == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field assignment from value", "[FieldAssignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field(1);
    temporary_field = 1;
    REQUIRE(temporary_field == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field assignment from zero", "[FieldAssignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field(0);
    temporary_field = 0;
    REQUIRE(temporary_field == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field get element", "[FieldAssignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    MyType temporary_field_element = temporary_field.get_element({2,2,2});
    MyType dummy_field_element = dummy_field.get_element({2,2,2});
    REQUIRE(temporary_field_element == dummy_field_element);
}

TEST_CASE_METHOD(FieldTest, "Field set element", "[FieldAssignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    temporary_field.set_element(100,{2,2,2});
    dummy_field.set_element(100,{2,2,2});
    REQUIRE(temporary_field == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field arithmetic with constant", "[FieldMathematicalOperations]") {
    fill_dummy_field(2);
    CHECK((temporary_field(1)+=1) == dummy_field);
    CHECK((temporary_field(3)-=1) == dummy_field);
    CHECK((temporary_field(1)*=2) == dummy_field);
    CHECK((temporary_field(4)/=2) == dummy_field);
    CHECK((temporary_field(1)+ 1) == dummy_field);
    CHECK((temporary_field(3)- 1) == dummy_field);
    CHECK((temporary_field(1)* 2) == dummy_field);
    REQUIRE((temporary_field(4)/ 2) == dummy_field);
}

TEST_CASE_METHOD(FieldTest, "Field arithmetic with field", "[FieldMathematicalOperations]") {
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

TEST_CASE_METHOD(FieldTest, "Field squarenorm", "[FieldMathematicalOperations]") {
    fill_dummy_field(1);
    REQUIRE(dummy_field.squarenorm()/lattice->volume() == 1);
}

TEST_CASE_METHOD(FieldTest, "Field sum reduction", "[FieldMathematicalOperations]") {
    fill_dummy_field(1);
    REQUIRE(dummy_field.sum() == lattice->volume());
}

TEST_CASE_METHOD(FieldTest, "Field product reduction", "[FieldMathematicalOperations]") {
    fill_dummy_field(1);
    REQUIRE(dummy_field.product() == 1);
}

TEST_CASE_METHOD(FieldTest, "Field MinMax", "[FieldMathematicalOperations]") {
    fill_dummy_field(2);
    dummy_field.set_element(1.0,{2,2,1});
    dummy_field.set_element(3.0,{2,2,2});
    CoordinateVector loc_min,loc_max;
    REQUIRE(dummy_field.min() == 1.0);    
    REQUIRE(dummy_field.max() == 3.0);
    REQUIRE(dummy_field.min(ODD) == 1.0);    
    REQUIRE(dummy_field.max(EVEN) == 3.0);
    REQUIRE(dummy_field.min(EVEN,loc_min) == 2.0);    
    REQUIRE(dummy_field.max(ODD,loc_max) == 2.0);
}