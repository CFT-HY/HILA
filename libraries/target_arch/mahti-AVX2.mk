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
CXXFLAGS  := -O3 -x c++ --std=c++17 -march=native -mtune=znver1 -mavx2 -mfma -m3dnow -fomit-frame-pointer
#CXXFLAGS := -g -x c++ --std=c++17

# Define this to use setup_layout_vector
# it works for non-AVX code too, but is necessary for AVX

LAYOUT_VECTOR := 1

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

### Need to give include directory to mpi for hilapp - here 2 common ones
# MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
#
# This in general works with OpenMPI: --showme:incdirs gives the include path of the mpic++
MPI_INCLUDE_DIRS := $(addprefix -I, $(shell  $(CC) --showme:incdirs) )

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
HILAPP_INCLUDE_LIST += $(MPI_INCLUDE_DIRS)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`


# Linker libraries and possible options

LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS :=

# These variables must be defined here
#
HILAPP_OPTS := -target:AVX $(HILAPP_INCLUDES)
HILA_OPTS := -DAVX -DUSE_MPI -DNODE_LAYOUT_BLOCK=128

