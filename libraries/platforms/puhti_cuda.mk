# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
# Before compiling, load gcc and cuda modules:
# module load gcc/9.1.0 cuda/11.1.0 openmpi/4.0.5-cuda
#

# Define compiler
CC = nvcc
LD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

# Define compilation flags
CXXFLAGS = -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA -DUSE_MPI  $(CXX_INCLUDE)
#CXXFLAGS = -g -x c++ --std=c++17 

# No need to give include directory to mpi for hilapp - here 2 common ones
MPI_INCLUDE_DIRS = 

MPI_LIBS =  -lmpi

LDLIBS = -lcufft -lm $(MPI_LIBS)

# extra cuda objects here
HILA_OBJECTS += build/hila_cuda.o

################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA -DCUDA -DUSE_MPI -DPUHTI $(MPI_INCLUDE_DIRS)
HILA_OPTS = -DUSE_MPI -DCUDA -DPUHTI


