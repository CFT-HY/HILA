# Platform specific makefile for vanilla (linux) mpi code 
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

### Define compiler and options

# Define compiler
CC = mpic++
LD = mpic++

# Define compilation flags
CXXFLAGS = -O3 -x c++ --std=c++17 
#CXXFLAGS = -g -x c++ --std=c++17 

# No need to give include directory to mpi for hilapp - here 2 common ones
MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include

################

# Linker libraries and possible options

LDLIBS  = -lfftw3 -lm
LDFLAGS = 

# These variables must be defined here
#
HILAPP_OPTS = -target:AVX -DAVX -DUSE_MPI $(MPI_INCLUDE_DIRS)
HILA_OPTS = -DAVX -DUSE_MPI

