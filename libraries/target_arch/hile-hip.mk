# Platform specific makefile for LUMI standard (CPU) code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

$(info ########################################################################)
$(info Target hile-hip: currently requires the module )
$(info module load craype-accel-amd-gfx90a )
$(info on HILE gpu node. )
$(info ########################################################################)

### Define compiler and options

# Define compiler - use cray CC wrapper
CC := hipcc
LD := hipcc

ROCM_PATH := /opt/rocm

# Define compilation flags
#CXXFLAGS  := -Ofast -flto -x c++ --std=c++17 -fno-rtti
#CXXFLAGS := -g -x c++ --std=c++17
# CXXFLAGS := -std=c++17 -fno-rtti --rocm-path=${ROCM_PATH} --offload-arch=gfx908 -x hip -fgpu-rdc
CXXFLAGS := -std=c++17 -fno-rtti -O3 -xhip -fgpu-rdc --offload-arch=gfx90a 
CXXFLAGS += -D__HIP_PLATFORM_AMD__=1 -D__HIP_PLATFORM_HCC__=1 
CXXFLAGS += -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1
# CXXFLAGS := -std=c++17 --offload-arch=gfx908 -x c++
#
# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | hipcc -xc++ --std=c++20 -Wp,-v - 2>&1 | grep "^ "))

# stddef.h again!
# HILAPP_INCLUDE_LIST += -I/opt/cray/pe/gcc/default/snos/lib/gcc/x86_64-suse-linux/default/include -I/opt/cray/pe/fftw/default/x86_64/include
HILAPP_INCLUDE_LIST += -I/opt/cray/pe/fftw/default/x86_64/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`

HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o

# ROCM_LIBS := $(shell echo ${ROCM_PATH} | sed s/rocm/rocmlibs/)

# Currently LUMI EAP things require explicit setting of many include paths -- NOTE: MPI may be bad
HILA_INCLUDES := -I${ROCM_PATH}/include -I${ROCM_PATH}/include/hip
HILA_INCLUDES += -I${ROCM_PATH}/include/rocrand -I${ROCM_PATH}/include/hiprand 
HILA_INCLUDES += -I${ROCM_PATH}/include/hipfft -I${ROCM_PATH}/include/hipcub
HILA_INCLUDES += -I${CRAY_MPICH_DIR}/include


# Linker libraries and possible options

# LDLIBS  := -lfftw3 -lfftw3f -lm
#LDFLAGS := $(CXXFLAGS) -fgpu-rdc --hip-link --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64
LDFLAGS := $(CXXFLAGS) -fgpu-rdc --hip-link --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64
LDFLAGS += -L${ROCM_PATH}/hipfft/lib/ -lhipfft
LDFLAGS += -L${CRAY_MPICH_DIR}/lib ${PE_MPICH_GTL_DIR_amd_gfx90a} -lmpi ${PE_MPICH_GTL_LIBS_amd_gfx90a}

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP $(HILAPP_INCLUDES)
HILA_OPTS := -DHIP -DGPU_VECTOR_REDUCTION_THREAD_BLOCKS=64 -DGPU_RNG_THREAD_BLOCKS=64 $(HILA_INCLUDES)



