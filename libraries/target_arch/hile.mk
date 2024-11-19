# Platform specific makefile for HILE standard (CPU) code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

$(info ########################################################################)
$(info Target hile: remember to )
$(info   module load cray-mpich cray-fftw )
$(info ########################################################################)


### Define compiler and options

# Define compiler - use cray CC wrapper
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


# stddef.h again!
HILAPP_INCLUDE_LIST +=  -I/opt/cray/pe/fftw/default/x86_rome/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`


# Linker libraries and possible options

LDLIBS  := -lstdc++fs -lfftw3 -lfftw3f -lm -lstdc++fs
LDFLAGS := $(CXXFLAGS)

# These variables must be defined here
#
HILAPP_OPTS := $(HILAPP_INCLUDES)
HILA_OPTS := -DNODE_LAYOUT_BLOCK=128



