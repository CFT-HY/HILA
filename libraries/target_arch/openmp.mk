# Platform specific makefile for vanilla+openmp (linux) mpi code 
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

### Define compiler and options

# Define compiler
CC := mpic++
LD := mpic++

# Define compilation flags
CXXFLAGS  := -O3 -x c++ --std=c++17 -fno-rtti -mavx2 -mfma -fopenmp
# CXXFLAGS := -g -x c++ --std=c++17


## The following incantation gives the include paths of the $(CC) compiler (if it is gcc or clang)
# It may be that this path is not necessary at all, usually not for "system installed" clang
# STD_INCLUDE_DIRS := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ /"))
STD_INCLUDE_DIRS :=


# Linker libraries and possible options

LDLIBS  := -lfftw3 -lfftw3f -lm -lgomp
LDFLAGS :=

# These variables must be defined here
#
HILAPP_OPTS := -target:openmp $(STD_INCLUDE_DIRS) $(MPI_INCLUDE_DIRS)
HILA_OPTS := -DUSE_MPI



