# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

# Define compiler -- NOTE: NEED AT LEAST CUDA 11 TO COMPILE c++17

CC = /usr/local/cuda-11/bin/nvcc
# CC = /usr/bin/nvcc
LD = $(CC) -std c++17

# Define compilation flags - 61 and 52 work with fairly common geForce cards
CXXFLAGS = -O3 -dc -x cu -std c++17 -DCUDA 
CXXFLAGS += -gencode arch=compute_61,code=sm_61 -gencode arch=compute_52,code=sm_52  --use_fast_math --restrict
CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177"
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


