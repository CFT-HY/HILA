#!/bin/bash

# Get a list from command line arguments or list all benchmark files
if [ "$#" -gt 0 ]; then
    benchmarks=( "$@" )
else
    benchmarks=$(ls bench_*.cpp  )
fi


make cleanall
for benchfile in $benchmarks; do
    benchmark="${benchfile%.*}"
    make -s ${benchmark}.exe 2>/dev/null
    ./${benchmark}.exe
    rm -f ${benchmark}.exe
done
make clean

exit 0
