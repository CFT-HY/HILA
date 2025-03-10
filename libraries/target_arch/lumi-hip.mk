# Platform specific makefile for LUMI standard (CPU) code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

$(info ########################################################################)
$(info Target lumi-hip: remember to )
$(info module load CrayEnv PrgEnv-cray craype-accel-amd-gfx90a cray-mpich rocm )
$(info lumi-hip DOES NOT WORK AFTER UPDATE TO ROCM 6.0.1, USE lumi-hip-CC INSTEAD )
$(info SOMEBODY PLEASE REPAIR THIS )
$(info ########################################################################)


### Define compiler and options

# Define compiler - use cray CC wrapper
CC := hipcc
LD := hipcc

# Define compilation flags
#CXXFLAGS  := -Ofast -flto -x c++ --std=c++17 -fno-rtti
#CXXFLAGS := -g -x c++ --std=c++17
# CXXFLAGS := -std=c++17 -fno-rtti --rocm-path=${ROCM_PATH} --offload-arch=gfx908 -x hip -fgpu-rdc
CXXFLAGS := -std=c++17 -fno-rtti --rocm-path=${ROCM_PATH} --offload-arch=gfx90a -x hip -fgpu-rdc -Wno-cuda-compat
# CXXFLAGS := -std=c++17 --offload-arch=gfx908 -x c++

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | $(CC) -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

# stddef.h again!
# HILAPP_INCLUDE_LIST += -I/opt/cray/pe/gcc/default/snos/lib/gcc/x86_64-suse-linux/default/include -I/opt/cray/pe/fftw/default/x86_64/include
HILAPP_INCLUDE_LIST += -I/opt/cray/pe/fftw/default/x86_64/include


# Write hilapp inlcudes to a file 0hilapp_incl_dirs
$(shell mkdir -p build)
$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`

# use memory_pool.o, not memory_pool2.o - fails in gpu-aware mpi!
HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o

# ROCM_LIBS := $(shell echo ${ROCM_PATH} | sed s/rocm/rocmlibs/)

# Currently LUMI EAP things require explicit setting of many include paths -- NOTE: MPI may be bad
HILA_INCLUDES := -I${ROCM_PATH}/hip/include
HILA_INCLUDES += -I${ROCM_PATH}/rocrand/include/ -I${ROCM_PATH}/hiprand/include/ 
HILA_INCLUDES += -I${ROCM_PATH}/hipfft/include/
HILA_INCLUDES += -I${MPICH_DIR}/include -I/opt/cray/pe/fftw/default/x86_64/include

# Define FFTW_LIB for use in app Makefiles, if fftw is linked in (add $(FFTW_LIB) to link command
FFTW_LIB := -L/opt/cray/pe/fftw/default/x86_64/lib/ -lfftw3

# Linker libraries and possible options

# LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS := $(CXXFLAGS) -fgpu-rdc --hip-link --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib -lamdhip64 
LDFLAGS += -L${ROCM_PATH}/hipfft/lib/ -lhipfft
LDFLAGS += -L${MPICH_DIR}/lib -lmpi -L${CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP $(HILAPP_INCLUDES)
HILA_OPTS := -DHIP $(HILA_INCLUDES)



