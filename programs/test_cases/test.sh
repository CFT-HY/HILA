#!/bin/bash

tests=( "$@" )

for test in "${tests[@]}"; do
    ./$test
    if [ $? -eq 0 ];
    then
        echo -e "$test \e[36msuccess"
    else
        echo -e "$test \e[30mfailed"
    fi
done
