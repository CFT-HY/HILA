#!/bin/bash

tests=( "$@" )

for test in "${tests[@]}"; do
    ../../build/transformer $test
    if [ $? -eq 0 ];
    then
        echo -e "transformer $test \e[36msuccess"
    else
        echo -e "$test \e[30mfailed"
	exit 1
    fi
done

exit 0
