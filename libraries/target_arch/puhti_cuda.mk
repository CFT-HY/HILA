# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
# Before compiling, load gcc and cuda modules:
# module load gcc/9.1.0 cuda/11.1.0 openmpi/4.0.5-cuda
#

$(info ########################################################################)
$(info Target puhti_cuda:  remember to )
$(info   module load gcc/9.1.0 cuda/11.1.0 openmpi/4.0.5-cuda )
$(info ########################################################################)

# Define compiler
CC = nvcc
LD = nvcc -gencode arch=compute_70,code=sm_70 --use_fast_math --restrict

# Define compilation flags
CXXFLAGS = -dc -O3 -std=c++17 -x cu -gencode arch=compute_70,code=sm_70 --use_fast_math --restrict 
# 3162 is a warning about ignored inline in __global__ functions - it's not really ignored by nvcc,
# it allows definition of a function in multiple compilation units as required by c++ standard!!  
# Quiet it.
# Warning 177 is about unused variables 
CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=3162"

#CXXFLAGS = -g -x c++ --std=c++17 

STD_INCLUDE_DIRS := $(addprefix -I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

# No need to give include directory to mpi for hilapp - here 2 common ones
MPI_INCLUDE_DIRS = 

MPI_LIBS =  -lmpi

LDLIBS = -lcufft -lm $(MPI_LIBS)

# extra cuda objects here
HILA_OBJECTS += build/hila_cuda.o

################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA $(STD_INCLUDE_DIRS) $(MPI_INCLUDE_DIRS)
HILA_OPTS = -DUSE_MPI -DCUDA -DPUHTI


