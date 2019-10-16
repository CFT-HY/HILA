#!/bin/bash

check(){
    if [ $? -eq 0 ];
    then
        echo -e "$1 $test \e[36msuccess\e[39m"
    else
        echo -e "$1 $test \e[30mfailede\e[39m"
    fi
}

transform_c(){
    echo ../../build/transformer $@ 2>/dev/null
    ../../build/transformer $@ 2>/dev/null
    check transform
}

compile_c(){
    echo make $@
    make $@
    check compile
}

tests=( "$@" )

for test in "${tests[@]}"; do
    transform_c ${test}.cpp
    compile_c ${test}.cpt
    echo ./$test
    ./$test
    check
    make clean
done

