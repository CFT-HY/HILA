#!/bin/bash

export MAKEFILE=Makefile
export RUNNER="mpirun -n 1"

for arg in "$@"
do
  if [ "$1" = "CUDA" ]; then
    shift
    export MAKEFILE=Makefile_gpu
    export RUNNER=
  fi
  if [ "$1" = "AVX" ]; then
    shift
    export MAKEFILE=Makefile_avx
  fi
  if [ "$1" = "MPI" ]; then
    shift
    if [ "$1" = "-n" ]; then
      shift
      export RUNNER="mpirun -n $1"
      shift
    else
      export RUNNER=mpirun 
    fi
  fi
done


# Get a list from command line arguments or list all benchmark files
if [ "$#" -gt 0 ]; then
    benchmarks=( "$@" )
else
    benchmarks=$(ls bench_*.cpp  )
fi

make cleanall
for benchfile in $benchmarks; do
    benchmark="${benchfile%.*}"
    echo make -j -f $MAKEFILE -s ${benchmark}.exe
    make -j -f $MAKEFILE -s ${benchmark}.exe 2>/dev/null
    echo ${RUNNER} ./${benchmark}.exe
    ${RUNNER} ./${benchmark}.exe 
done

exit 0
