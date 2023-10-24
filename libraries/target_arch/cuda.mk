# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
# Show option in "make help"
#% ARCH=cuda options: use with make [...] options
#%    CUDA_VERSION=99.9   - Set cuda version, if installed (default=11.6)
#%    CUDA_ARCH=<version> - Compile with compute_NN and sm_NN (default=61)
#%    DEBUG=on            - Compile with -g, for debugging
#%

# Define compiler -- NOTE: NEED AT LEAST CUDA 11 TO COMPILE c++17
ifndef CUDA_VERSION
	CUDA_VERSION = 11.6
endif
CC := /usr/local/cuda-${CUDA_VERSION}/bin/nvcc
# CC = /usr/bin/nvcc
LD := $(CC) -std c++17

# Define compilation flags - 61 and 52 work with fairly common geForce cards
ifndef CUDA_ARCH
	CUDA_ARCH = 61
endif

CXXFLAGS := -dc -x cu -std c++17 -DCUDA -gencode arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH} --restrict
ifndef DEBUG
	CXXFLAGS += -O3 --use_fast_math
else
	CXXFLAGS += -g
endif

# -gencode arch=compute_52,code=sm_52

#
# 20050 is a warning about ignored inline in __global__ functions - it's not ignored though, it allows multiple
# definitions as per c++ standard!
CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=20050"
#CXXFLAGS = -g -x c++ --std=c++17 

LDLIBS := -L/usr/local/cuda-${CUDA_VERSION}/targets/x86_64-linux/lib/ -lcudart -lcufft -lm 

# Need to give include directory to mpi for hilapp and nvcc - here 2 common ones
MPI_INCLUDE_DIRS := -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
MPI_LIBS := -L/usr/lib/openmpi/lib -lmpi

LDLIBS += $(MPI_LIBS)
LDFLAGS += -gencode arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}

# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o

################

# These variables must be defined here
HILAPP_OPTS := -target:CUDA 
HILA_OPTS := -DCUDA $(MPI_INCLUDE_DIRS)

