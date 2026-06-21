# Platform specific makefile for cuda code.  Tuned for "Roihu" computer at CSC
#
# this is included from main.mk -file, which is in turn included from
# application makefile

$(info ########################################################################)
$(info Target roihu-gpu: Remember to load modules:)
$(info module load gcc/15.2.0 openmpi/5.0.8 cuda/13.1.1)
$(info ########################################################################)

# Define compiler
CC = nvcc
LD = nvcc -gencode arch=compute_90a,code=sm_90a --use_fast_math --restrict

# Define compilation flags  - arch 90 or 90a. 90a has GH extension features but is not backwards/forwards compatible
CXXFLAGS = -g -dc -std=c++17 -x cu -gencode arch=compute_90a,code=sm_90a --restrict
CXXFLAGS_NOOPT = $(CXXFLAGS) -O1
CXXFLAGS += -O3 --use_fast_math

# Suppress some warnings
# 3162 is a warning about ignored inline in __global__ functions - it's not really ignored by nvcc,
# it allows definition of a function in multiple compilation units as required by c++ standard!!
# Warning 177 is about unused variables
#CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=3162"
CXXFLAGS += -Xcudafe "--diag_suppress=177"


MPI_LIBS =  -lmpi

LDLIBS = -lcufft -lm $(MPI_LIBS)

# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o

################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA $(STD_HILAPP_INCLUDES)
HILA_OPTS = -DCUDA
