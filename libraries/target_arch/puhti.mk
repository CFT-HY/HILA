# Platform specific makefile for puhti (cluster) mpi code
# Puhti does not have clang, so we use statically compiled hilapp
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#


$(info ########################################################################)
$(info Target puhti:  remember to )
$(info   module load gcc/9.1.0 fftw openmpi/4.1.1-cuda )
$(info ########################################################################)


### Define compiler and options

# Define compiler
CC := mpic++
LD := mpic++

# Define compilation flags
CXXFLAGS := -O3 -x c++ --std=c++17
#CXXFLAGS := -g -x c++ --std=c++17 

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

# Write hilapp includes to a file 0hilapp_incl_dirs
$(shell mkdir -p build; echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`



# Linker libraries and possible options

LDLIBS  = -lfftw3 -lm
LDFLAGS = 

# These variables must be defined here
#
HILAPP_OPTS = $(HILAPP_INCLUDES)
HILA_OPTS =



