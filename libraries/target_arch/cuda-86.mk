# Platform specific makefile for cuda code.  Tuned for CFT group server (Nvidia RTX A4000)
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

# Define compiler -- NOTE: NEED AT LEAST CUDA 11 TO COMPILE c++17
ifndef CUDA_VERSION
	CUDA_VERSION = 11.8
endif
CC := /usr/local/cuda-${CUDA_VERSION}/bin/nvcc
# CC = /usr/bin/nvcc
LD := $(CC) -std c++17

# Define compilation flags - 61 and 52 work with fairly common geForce cards
ifndef CUDA_ARCH
	CUDA_ARCH = 86
endif
CXXFLAGS := -O3 -dc -x cu -std c++17 -DCUDA 
CXXFLAGS += -gencode arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH} --use_fast_math --restrict

#
# 20050 is a warning about ignored inline in __global__ functions - it's not ignored though, it allows multiple
# definitions as per c++ standard!
CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=20050"
#CXXFLAGS = -g -x c++ --std=c++17 

LDLIBS := -L/usr/local/cuda-${CUDA_VERSION}/targets/x86_64-linux/lib/ -lcudart -lcufft -lm 

# Need to give include directory to mpi for hilapp and nvcc - here 2 common ones
MPI_INCLUDE_DIRS := -I/opt/local/cuda-aware-mpi/include
MPI_LIBS := -L/opt/local/cuda-aware-mpi/lib -lmpi

LDLIBS += $(MPI_LIBS)
LDFLAGS += -gencode arch=compute_${CUDA_ARCH},code=sm_${CUDA_ARCH}

# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o
################

# These variables must be defined here
HILAPP_OPTS := -target:CUDA 
HILA_OPTS := -DCUDA $(MPI_INCLUDE_DIRS)


