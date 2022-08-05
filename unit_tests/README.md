# Unit tests

This is a unit test sweet for HILA software meant to test full functionalities of HILA applications.

The test sweet is built using [Catch2](https://github.com/catchorg/Catch2/tree/v2.x) a c++ unit testing framework. Catch2 v2.x works by simply adding a header to the c++ project which is in the path `./unit_tests/src/catch.hpp`. This is a very light weight approach and requires zero installation from user.

## Running tests

Simply compile and run

```bash
make -j4
./build/catch_main
```

Example output with all tests passing:

```
===============================================================================
All tests passed (36 assertions in 17 test cases)
```

For more verbose output run

    ./build/catch_main --success

## Selecting specific tests

To list all tests run:

    ./build/catch_main -l

To list all tags run:

    ./biuld/catch_main -t

If you want to select a specific test, for example testing field assignment from another field (./unit_tests/src/test_field.cpp):

```c++
...
  void fill_dummy_field() {
        onsites(ALL) this->dummy_field[X] = hila::random();
    }
...
TEST_CASE_METHOD(FieldTest, "Assignment from field", "[Assignment]") {
    Field<MyType> temporary_field;
    fill_dummy_field();
    temporary_field = dummy_field;
    REQUIRE(temporary_field == dummy_field);
}
```

One can run:

    ./build/catch_main "Assignment from field" --success


```
-------------------------------------------------------------------------------
Assignment from field
-------------------------------------------------------------------------------
build/test_field.cpt:7050
...............................................................................

build/test_field.cpt:7054: PASSED:
  REQUIRE( temporary_field == dummy_field )
with expansion:
  {?} == {?}
```

If you want to select all tests under the tag `[Assignment]` then run:

    ./build/catch_main [Assignment]

```
Filters: [Assignment]
===============================================================================
All tests passed (5 assertions in 5 test cases)
```

To list all tests under the `[Assignment]` tag run:

    ./build/catch_main [Assignment] -l

```
Matching test cases:
  Assignment from field
      [Assignment]
  Assignment from value
      [Assignment]
  Assignment from zero
      [Assignment]
  Get element
      [Assignment]
  Set element
      [Assignment]
5 matching test cases
```

## NOTES

test_lattice.cpp:

- discuss create_std_gathers() and nn_comminfo_struct
    - why is nn_comminfo_struct.from_node and nn_comminfo_struct.to_node the same