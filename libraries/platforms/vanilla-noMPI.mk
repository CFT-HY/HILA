# Platform specific makefile for vanilla (linux) mpi code 
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

### Define compiler and options

# Define compiler
CC := g++
LD := g++

# Define compilation flags
CXXFLAGS := -O3 -x c++ --std=c++17
#CXXFLAGS := -g -x c++ --std=c++17
LDLIBS := -lfftw3 -lm

STD_INCLUDE_DIRS := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))


# No need to give include directory to mpi for hilapp - here 2 common ones
# MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include

################

# These variables must be defined here
#
HILAPP_OPTS := -no-mpi $(STD_INCLUDE_DIRS)
HILA_OPTS :=

