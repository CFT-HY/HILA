# Platform specific makefile for LUMI Cray Env (LUMI/yy.xx) AMD Mi250 GPU code
#
# this is included from main.mk -file, which is in turn included from
# application makefile
#

$(info ########################################################################)
$(info Target lumi-hip-CC: remember to )
$(info   module purge )
$(info   module load PrgEnv-cray craype-accel-amd-gfx90a rocm/6.3.4)
$(info ########################################################################)


### Define compiler and options

# c++ standard level, can be set in makefile or command line
CPPSTD := c++17

# Define compiler - use cray CC wrapper
CC := CC
LD := CC

# Flags for both the compiler and the linker
CXX_LD_FLAGS := -fgpu-rdc --offload-arch=gfx90a
# Compiler flags
CXXFLAGS := $(CXX_LD_FLAGS) -std=$(CPPSTD) -fno-rtti -O3 -xhip
## -fno-omit-frame-pointer -O3 -xhip

# ??? HILA seems to internally compile something with this...
CXXFLAGS_NOOPT := $(CXX_LD_FLAGS) -std=$(CPPSTD) -fno-rtti -xhip

#
# HILAPP needs stdlib includes for whatever reason... The following asks the compiler to print them
HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=$(CPPSTD) -Wp,-v - 2>&1 | grep "^ "))

# stddef.h again!
#HILAPP_INCLUDE_LIST += -I/opt/cray/pe/fftw/default/x86_64/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`

HILA_OBJECTS += build/hila_gpu.o

# Linker config
LDFLAGS := $(CXX_LD_FLAGS) --hip-link
LDFLAGS += -lhipfft
# Cray _should_ link GTL if an accel module is loaded but fails to do so here. Link it manually:
LDFLAGS += -lmpi_gtl_hsa

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP $(HILAPP_INCLUDES)
HILA_OPTS := -DHIP -DGPU_RNG_THREAD_BLOCKS=64 $(HILA_INCLUDES)
