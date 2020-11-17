# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

# Define compiler
CC = /usr/local/cuda-9.2/bin/nvcc
LD = /usr/local/cuda-9.2/bin/nvcc -gencode arch=compute_35,code=sm_35 -fmad=false

# Define compilation flags
CXXFLAGS = -dc -x cu -gencode arch=compute_35,code=sm_35 -fmad=false -std=c++11 -DCUDA 
#CXXFLAGS = -g -x c++ --std=c++17 

LDLIBS = -lfftw3 -lm 

# Need to give include directory to mpi for hilapp and nvcc - here 2 common ones
MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
MPI_LIBS = -L/usr/lib/openmpi/lib -lmpi

LDLIBS = -lm $(MPI_LIBS)

# extra cuda objects here
HILA_OBJECTS += build/hila_cuda.o

################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA 
HILA_OPTS = -DCUDA -DUSE_MPI $(MPI_INCLUDE_DIRS)


