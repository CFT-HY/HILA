#!/bin/bash

#We need these for a cpu build. Turn these into a library?
export BUILD_REQUIREMENTS="mersenne_inline.o lattice.o memory.o setup_layout_generic.o map_node_layout_trivial.o"

check(){
    if [ $? -eq 0 ];
    then
        echo -e "$1 $test \e[36msuccess\e[39m"
    else
        echo -e "$1 $test \e[31mfailed\e[39m"
	    exit 1
    fi
}

transform_c(){
    echo ../../build/transformer $1
    ../../build/transformer $1 2>/dev/null
    check transform
}

compile_c(){
    make ${1}.o $BUILD_REQUIREMENTS
    clang++-8 -o $1 ${1}.o $BUILD_REQUIREMENTS
    check compile
}


# Get a list from command line arguments or list all test files
if [ "$#" -gt 0 ]; then
    tests=( "$@" )
else
    tests=$(ls test_*.cpp  )
fi


make cleanall
for testfile in $tests; do
    test="${testfile%.*}"
    echo $test
    transform_c ${test}.cpp
    compile_c ${test}
    echo ./$test
    ./$test
    check
    rm $test ${test}.cpt
done

exit 0
