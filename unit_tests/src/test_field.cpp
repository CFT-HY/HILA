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

TEST_CASE_METHOD(FieldTest, "Field nullptr, allocation and destructor", "[Field]") {
    SECTION("Test field nullptr") {
        REQUIRE_FALSE(dummy_field.is_allocated());
    }
    SECTION("Test field allocation") {
        dummy_field.allocate();
        REQUIRE(dummy_field.is_allocated());
    }
    SECTION("Test destructor") {
        dummy_field.~Field();
        REQUIRE_FALSE(dummy_field.is_allocated());
    }
}

TEST_CASE_METHOD(FieldTest, "Field constructors", "[Field]") {
    SECTION("Test field allocation") {
        dummy_field.allocate();
        REQUIRE(dummy_field.is_allocated());
    }
    SECTION("Test destructor") {
        dummy_field.~Field();
        REQUIRE_FALSE(dummy_field.is_allocated());
    }
    SECTION("Test scalar constructor") {
        fill_dummy_field(1);
        REQUIRE(temporary_field(1) == dummy_field);
    }
    SECTION("Test zero constructor") {
        fill_dummy_field(0);
        REQUIRE(temporary_field(0) == dummy_field);
    }
    SECTION("Test copy constructor") {
        fill_dummy_field();
        REQUIRE(temporary_field(dummy_field) == dummy_field);
    }
}

TEST_CASE_METHOD(FieldTest, "Field assignment", "[Field]") {
    Field<MyType> temporary_field;
    SECTION("Assignment from field") {
        fill_dummy_field();
        temporary_field = dummy_field;
        REQUIRE(temporary_field == dummy_field);
    }
    SECTION("Assignment from value") {
        fill_dummy_field(1);
        temporary_field = 1;
        REQUIRE(temporary_field == dummy_field);
    }
    SECTION("Assignment from zero") {
        fill_dummy_field(0);
        temporary_field = 0;
        REQUIRE(temporary_field == dummy_field);
    }
}

TEST_CASE_METHOD(FieldTest, "Field get and set element", "[FieldAssignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    SECTION("Field get element") {
        MyType temporary_field_element = temporary_field.get_element({2, 2, 2});
        MyType dummy_field_element = dummy_field.get_element({2, 2, 2});
        REQUIRE(temporary_field_element == dummy_field_element);
    }
    // SECTION("Field set element") {
    //     temporary_field.set_element(100,{2,2,2});
    //     dummy_field.set_element(100,{2,2,2});
    //     REQUIRE(temporary_field == dummy_field);
    // }
}

TEST_CASE_METHOD(FieldTest, "Field arithmetic", "[Field]") {
    fill_dummy_field(2);
    SECTION("Arithmetic with constant") {
        REQUIRE((temporary_field(1) += 1) == dummy_field);
        REQUIRE((temporary_field(3) -= 1) == dummy_field);
        REQUIRE((temporary_field(1) *= 2) == dummy_field);
        REQUIRE((temporary_field(4) /= 2) == dummy_field);
        REQUIRE((temporary_field(1) + 1) == dummy_field);
        REQUIRE((temporary_field(3) - 1) == dummy_field);
        REQUIRE((temporary_field(1) * 2) == dummy_field);
        REQUIRE((temporary_field(4) / 2) == dummy_field);
    }
    SECTION("Arithmetic with field") {
        REQUIRE((temporary_field(1) += temporary_field(1)) == dummy_field);
        REQUIRE((temporary_field(3) -= temporary_field(1)) == dummy_field);
        REQUIRE((temporary_field(1) *= temporary_field(2)) == dummy_field);
        REQUIRE((temporary_field(4) /= temporary_field(2)) == dummy_field);
        REQUIRE((temporary_field(1) + temporary_field(1)) == dummy_field);
        REQUIRE((temporary_field(3) - temporary_field(1)) == dummy_field);
        REQUIRE((temporary_field(1) * temporary_field(2)) == dummy_field);
        REQUIRE((temporary_field(4) / temporary_field(2)) == dummy_field);
    }
}

// UNARY OPERATOR?

TEST_CASE_METHOD(FieldTest, "Field mathematical operations", "[Field]") {
    fill_dummy_field(1);
    SECTION("Squarenorm") {
        REQUIRE(dummy_field.squarenorm() / lattice.volume() == 1);
    }
    SECTION("Sum reduction") {
        REQUIRE(dummy_field.sum() == lattice.volume());
    }
    SECTION("Product reduction") {
        REQUIRE(dummy_field.product() == 1);
    }
    SECTION("MinMax") {
        dummy_field[{2, 2, 2}] = 2.0;
        dummy_field[{2, 2, 2}] = 2.0;
        CoordinateVector loc_min, loc_max;
        REQUIRE(dummy_field.min() == 0.0);
        REQUIRE(dummy_field.max() == 2.0);
        REQUIRE(dummy_field.min(ODD) == 0.0);
        REQUIRE(dummy_field.max(EVEN) == 2.0);
        REQUIRE(dummy_field.min(EVEN, loc_min) == 1.0);
        REQUIRE(dummy_field.max(ODD, loc_max) == 1.0);
    }
}