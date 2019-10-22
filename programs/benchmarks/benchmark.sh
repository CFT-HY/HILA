#!/bin/bash

if [ "$1" = "-f" ]; then
    shift
    export MAKEFILE=$1
    shift
else
    export MAKEFILE=Makefile
fi


# Get a list from command line arguments or list all benchmark files
if [ "$#" -gt 0 ]; then
    benchmarks=( "$@" )
else
    benchmarks=$(ls bench_*.cpp  )
fi

make cleanall
for benchfile in $benchmarks; do
    benchmark="${benchfile%.*}"
    make -f $MAKEFILE -s ${benchmark}.exe 2>/dev/null
    ./${benchmark}.exe
done

exit 0
