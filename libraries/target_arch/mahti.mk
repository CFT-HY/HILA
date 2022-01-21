# Platform specific makefile for vanilla (linux) mpi code 
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

$(info ########################################################################)
$(info Target mahti: remember to )
$(info   module load cuda gcc/9.3.0 openmpi/4.0.3 amdfftw/2.2-mpi)
$(info ########################################################################)


### Define compiler and options

# Define compiler
CC := mpic++
LD := mpic++

# Define compilation flags
CXXFLAGS  := -O3 -x c++ --std=c++17 -fno-rtti -mavx2 -mfma
#CXXFLAGS := -g -x c++ --std=c++17

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))


# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`


# Linker libraries and possible options

LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS :=

# These variables must be defined here
#
HILAPP_OPTS := $(HILAPP_INCLUDES)
HILA_OPTS := -DUSE_MPI -DNODE_LAYOUT_BLOCK=128



