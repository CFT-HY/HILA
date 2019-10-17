#!/bin/bash

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
    make -s ${1}.exe 2>/dev/null
    check compile
}

run_c(){
    echo ./$1
    ./${1}.exe
    check
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
    run_c ${test}
    rm ${test}.exe ${test}.cpt
done
make clean

exit 0
