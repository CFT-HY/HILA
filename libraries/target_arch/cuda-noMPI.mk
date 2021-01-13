# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

# Define compiler
CC = nvcc
LD = nvcc 

# Define compilation flags
CXXFLAGS = -dc -x cu -std c++14 -gencode arch=compute_30,code=sm_30 -DCUDA 
#CXXFLAGS = -g -x c++ --std=c++17 

LDLIBS = -lfftw3 -lm 

# No need to give include directory to mpi for hilapp - here 2 common ones

LDLIBS = -lm

# extra cuda objects here
HILA_OBJECTS += build/hila_cuda.o

################

# These variables must be defined here
#
HILAPP_OPTS = -no-mpi -target:CUDA -DCUDA 
HILA_OPTS = -DCUDA


