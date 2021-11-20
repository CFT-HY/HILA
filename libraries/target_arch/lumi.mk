# Platform specific makefile for vanilla (linux) mpi code 
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

$(info ########################################################################)
$(info Target lumi: remember to )
$(info   module load cray-fftw craype-hugepages2M )
$(info ########################################################################)


### Define compiler and options

# Define compiler
CC := CC
LD := CC

# Define compilation flags
CXXFLAGS  := -Ofast -flto -x c++ --std=c++17 -fno-rtti
#CXXFLAGS := -g -x c++ --std=c++17

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

### Need to give include directory to mpi for hilapp - here 2 common ones
# MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include
#
# This in general works with OpenMPI: --showme:incdirs gives the include path of the mpic++
#MPI_INCLUDE_DIRS := $(addprefix -I, $(shell  $(CC) --showme:incdirs) )
MPI_INCLUDE_DIRS := -I/opt/cray/pe/mpich/8.1.5/ucx/cray/9.1/include

# stddef.h again!
HILAPP_INCLUDE_LIST += -I/opt/cray/pe/gcc/default/snos/lib/gcc/x86_64-suse-linux/default/include -I/opt/cray/pe/fftw/default/x86_64/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
HILAPP_INCLUDE_LIST += $(MPI_INCLUDE_DIRS)
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`


# Linker libraries and possible options

LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS :=

# These variables must be defined here
#
HILAPP_OPTS := $(HILAPP_INCLUDES)
HILA_OPTS := -DUSE_MPI



