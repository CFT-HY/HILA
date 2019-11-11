#!/bin/bash

export ERRORLOG=$(mktemp /tmp/abc-script.XXXXXX)
export fails=0
export num_tests=0

MAKEFILE=Makefile

check(){
    if [ $? -eq 0 ];
    then
        echo -e "$1 $test \e[36msuccess\e[39m"
    else
        echo -e "$1 $test \e[31mfailed\e[39m"
        echo -e " * $1 $test" >> $ERRORLOG
        export fails=$((fails+1))
    fi
    export num_tests=$((num_tests+1))
}

transform_c(){
    echo ../../build/transformer $1
    ../../build/transformer $1 2>/dev/null
    check transform
}

compile_c(){
    make -f ${MAKEFILE} -s ${1}.exe 2>/dev/null
    check compile
}

run_c(){
    echo ./$1
    ./${1}.exe
    check run
}

run_mpi_c(){
    echo mpirun -n $2 ./${1}.exe
    mpirun -n $2 ./${1}.exe
    check run $2
}

# Get a list from command line arguments or list all test files
if [ "$#" -gt 0 ]; then
    tests=( "$@" )
else
    tests=$(ls test_*.cpp  )
fi

for D in 1 2 3 4 ; do
  sed -i 's/OPTS = .*/OPTS = -DNDIM='${D}'/' Makefile
  make -f ${MAKEFILE} clean
  for testfile in $tests; do
    test="${testfile%.*}"
    echo $test
    transform_c ${test}.cpp
    compile_c ${test}
    run_c ${test}
    run_mpi_c ${test} 1
    run_mpi_c ${test} 2
    run_mpi_c ${test} 4
    rm ${test}.exe ${test}.cpt
  done
  make -f ${MAKEFILE} clean
done


if [ ${fails} -eq 0 ]; then
	echo -e "\e[36m All tests passed \e[39m"
else
	echo ${fails}/${num_tests} tests failed
	cat $ERRORLOG
fi

rm $ERRORLOG
exit 0
