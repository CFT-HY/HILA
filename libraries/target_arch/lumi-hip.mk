# Platform specific makefile for LUMI standard (CPU) code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

$(info ########################################################################)
$(info Target lumi-hip: remember to )
$(info module load CrayEnv PrgEnv-cray craype-accel-amd-gfx908 cray-mpich rocm )
$(info ########################################################################)


### Define compiler and options

# Define compiler - use cray CC wrapper
CC := CC
LD := CC

# Define compilation flags
#CXXFLAGS  := -Ofast -flto -x c++ --std=c++17 -fno-rtti
#CXXFLAGS := -g -x c++ --std=c++17
CXXFLAGS := -std=c++17 -fno-rtti --rocm-path=${ROCM_PATH} --offload-arch=gfx908 -x hip -fgpu-rdc
# CXXFLAGS := -std=c++17 --offload-arch=gfx908 -x c++

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

# stddef.h again!
# HILAPP_INCLUDE_LIST += -I/opt/cray/pe/gcc/default/snos/lib/gcc/x86_64-suse-linux/default/include -I/opt/cray/pe/fftw/default/x86_64/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`

HILA_OBJECTS += build/hila_gpu.o build/memory_pool2.o

HILA_INCLUDES := -I/appl/eap/opt/rocm-4.3.1/rocrand/include/ -I/appl/eap/opt/rocm-4.3.1/hiprand/include/ 
HILA_INCLUDES += -I/appl/eap/opt/rocm-4.3.1/hipfft/include/

# Linker libraries and possible options

# LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS := $(CXXFLAGS) --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP $(HILAPP_INCLUDES)
HILA_OPTS := -DUSE_MPI -DHIP $(HILA_INCLUDES)



