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
ifndef DEBUG
	CXXFLAGS := -O3 -x c++ --std=c++17 -march=native -mavx2 -mfma -fabi-version=0 -fomit-frame-pointer
else
	CXXFLAGS := -g -x c++ --std=c++17 -march=native -mavx2 -mfma -fabi-version=0 -fomit-frame-pointer
endif

#CXXFLAGS := -g -x c++ --std=c++17

# Define this to use setup_layout_vector
# it works for non-AVX code too, but is necessary for AVX

LAYOUT_VECTOR := 1

## The following incantation gives the include paths of the $(CC) compiler (if it is gcc or clang)
# It may be that this path is not necessary at all, usually not for "system installed" clang
# THIS SEEMS TO CONFLICT WITH AVX DEFINITIONS; SO LEAVE OUT
#STD_INCLUDE_DIRS := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))
STD_INCLUDE_DIRS :=

################

# Linker libraries and possible options

LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS :=

# These variables must be defined here
#
HILAPP_OPTS := -target:AVX $(STD_INCLUDE_DIRS)
HILA_OPTS := -DAVX

