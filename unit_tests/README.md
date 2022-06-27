# Unit tests

This is a unit test sweet for HILA software meant to test full functionalities of HILA applications.

The test sweet is built using [GoogleTest](https://google.github.io/googletest/) a c++ unit testing framework. 

## Installing GoogleTest

Building GoogleTest requires cmake:

    sudo apt-get install cmake

Installing GoogleTest from source:
```bash
git clone https://github.com/google/googletest.git
cd googletest
mkdir build && cd build
cmake .. -DBUILD_SHARED_LIBS=ON -DINSTALL_GTEST=ON -DCMAKE_INSTALL_PREFIX:PATH=/usr
make -j8
sudo make install
sudo ldconfig
```

Test that the linker can find gtest libraries:

    ldconfig -p | grep gtest

Output should be:

```bash
libgtest_main.so.1.11.0 (libc6,x86-64) => /lib/x86_64-linux-gnu/libgtest_main.so.1.11.0
libgtest_main.so (libc6,x86-64) => /lib/x86_64-linux-gnu/libgtest_main.so
libgtest.so.1.11.0 (libc6,x86-64) => /lib/x86_64-linux-gnu/libgtest.so.1.11.0
libgtest.so (libc6,x86-64) => /lib/x86_64-linux-gnu/libgtest.so
```

## Running tests

Simply compile and run

```bash
make -j4 example_test
./build/example_test
```

Example output:

Note that the test is designed to fail for example purposes

```bash
[==========] Running 2 tests from 1 test suite.
[----------] Global test environment set-up.
[----------] 2 tests from LatticeTest
[ RUN      ] LatticeTest.Size
[       OK ] LatticeTest.Size (0 ms)
[ RUN      ] LatticeTest.SizeFail
build/tests.cpt:1912: Failure
Expected equality of these values:
  lattice->volume()
    Which is: 16777216
  255*256*256
    Which is: 16711680
[  FAILED  ] LatticeTest.SizeFail (0 ms)
[----------] 2 tests from LatticeTest (0 ms total)

[----------] Global test environment tear-down
[==========] 2 tests from 1 test suite ran. (0 ms total)
[  PASSED  ] 1 test.
[  FAILED  ] 1 test, listed below:
[  FAILED  ] LatticeTest.SizeFail

 1 FAILED TEST
```