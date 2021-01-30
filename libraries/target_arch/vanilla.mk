# Platform specific makefile for vanilla (linux) mpi code 
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
CXXFLAGS  := -O3 -x c++ --std=c++17 -fno-rtti -mavx2 -mfma
#CXXFLAGS := -g -x c++ --std=c++17


## The following incantation gives the include paths of the $(CC) compiler (if it is gcc or clang)
# It may be that this path is not necessary at all, usually not for "system installed" clang
# STD_INCLUDE_DIRS := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ /"))
STD_INCLUDE_DIRS :=

### Need to give MPI include directory for hilapp - here 2 common ones
# MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
#
# This in general works with OpenMPI: --showme:incdirs gives the include path of the mpic++
MPI_INCLUDE_DIRS := $(addprefix -I, $(shell  $(CC) --showme:incdirs) )


# Linker libraries and possible options

LDLIBS  := -lfftw3 -lm
LDFLAGS :=

# These variables must be defined here
#
HILAPP_OPTS := $(STD_INCLUDE_DIRS) $(MPI_INCLUDE_DIRS)
HILA_OPTS := -DUSE_MPI



