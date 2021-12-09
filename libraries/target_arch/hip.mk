# Platform specific makefile for HIP code, on linux with hipcc (ROCm) installed
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

# Define compiler -- NOTE: NEED AT LEAST CUDA 11 TO COMPILE c++17

CC := hipcc
# CC = /usr/bin/nvcc
LD := $(CC) --std=c++17 
#-gencode arch=compute_61,code=sm_61 -gencode arch=compute_52,code=sm_52

# Define compilation flags - 61 and 52 work with fairly common geForce cards
CXXFLAGS := -O3 --std=c++17 -x cu -dc --offload-arch=gfx1012
#--no-cuda-version-check -nocudalib 
# CXXFLAGS += -gencode arch=compute_61,code=sm_61 -gencode arch=compute_52,code=sm_52 
#
# 20050 is a warning about ignored inline in __global__ functions - it's not ignored though, it allows multiple
# definitions as per c++ standard!
# CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=20050"
#CXXFLAGS = -g -x c++ --std=c++17 

#LDLIBS := -L/usr/local/cuda-11.2/targets/x86_64-linux/lib/ -lcufft -lm 
LDLIBS := -lm

# Need to give include directory to mpi for hilapp and nvcc - here 2 common ones
MPI_INCLUDE_DIRS := -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
MPI_LIBS := -L/usr/lib/openmpi/lib -lmpi

LDLIBS += $(MPI_LIBS)

HIP_PATH ?= $(shell hipconfig --path)
HIP_INCLUDE_DIRS := -I$(HIP_PATH)/include -I$(HIP_PATH)/../hiprand/include -I$(HIP_PATH)/../hipfft/include

# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o

################

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -D__HIP_PLATFORM_AMD__
HILA_OPTS := -DHIP -DUSE_MPI $(MPI_INCLUDE_DIRS) $(HIP_INCLUDE_DIRS)


