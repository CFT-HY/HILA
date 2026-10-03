# Platform specific makefile for LUMI Cray Env (LUMI/yy.xx) AMD Mi250 GPU code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

$(info ########################################################################)
$(info Target lumi-hip-CC: remember to )
$(info   module --force purge && module --force unload LUMI )
$(info   module load LUMI/25.03 partition/G cpeCray/25.03 craype-accel-amd-gfx90a rocm )
$(info ########################################################################)

### Define compiler and options

# c++ standard level, can be set in makefile or command line
CPPSTD := c++17

# Define compiler - use cray CC wrapper
CC := CC
LD := CC

# Define compilation flags
CXXFLAGS := -std=$(CPPSTD) -fno-rtti -fgpu-rdc -xhip 

CXXFLAGS_NOOPT := $(CXXFLAGS) -O1
CXXFLAGS += -O3


# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

# when  module load LUMI/25.., the compiler wrapper CC
# is inface Clang, whcih know the pathes of Clang c++ 17 std include files.
# Then the good practice to get HILAPP_INCLUDE_LIST is call CC -xc++ --std=c++17 -Wp,-v - otherthan us g++

# HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))
HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | CC -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

# stddef.h again!
# HILAPP_INCLUDE_LIST += -I/opt/cray/pe/gcc/default/snos/lib/gcc/x86_64-suse-linux/default/include -I/opt/cray/pe/fftw/default/x86_64/include
HILAPP_INCLUDE_LIST += -I/opt/cray/pe/fftw/default/x86_64/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`

HILA_OBJECTS += build/hila_gpu.o

LDFLAGS := --hip-link
LDLIBS := -lhipfft

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP $(HILAPP_INCLUDES)
HILA_OPTS := -DHIP -DGPU_VECTOR_REDUCTION_THREAD_BLOCKS=64 -DGPU_RNG_THREAD_BLOCKS=64 $(HILA_INCLUDES)
