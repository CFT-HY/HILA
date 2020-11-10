# Platform specific makefile for puhti (cluster) mpi code
# Puhti does not have clang, so we use statically compiled hilapp
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
CXXFLAGS := -O3 -x c++ --std=c++17
#CXXFLAGS := -g -x c++ --std=c++17 


## Get the include dirs for stddef.h and stdarg.h etc. to hilapp
#STD_INCLUDE_DIRS := -I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/

STD_INCLUDE_DIRS := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

### Need to give include directory to mpi for hilapp - here 2 common ones
# MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
#
# This in general works with OpenMPI: --showme:incdirs gives the include path of the mpic++
MPI_INCLUDE_DIRS := $(addprefix -I, $(shell  $(CC) --showme:incdirs) )

# Linker libraries and possible options

LDLIBS  = -lfftw3 -lm
LDFLAGS = 

# These variables must be defined here
#
HILAPP_OPTS = $(STD_INCLUDE_DIRS) -DUSE_MPI $(MPI_INCLUDE_DIRS)
HILA_OPTS = -DUSE_MPI



