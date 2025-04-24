# Platform specific makefile for LUMI standard (CPU) code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#  1) amd/6.0.3                 8) craype-x86-trento
#  2) craype/2.7.31.11          9) craype-accel-amd-gfx90a
#  3) cray-dsmml/0.3.0         10) libfabric/1.15.2.0
#  4) cray-mpich/8.1.29        11) craype-network-ofi
#  5) cray-libsci/24.03.0      12) xpmem/2.8.2-1.0_5.1__g84a27a5.shasta
#  6) PrgEnv-amd/8.5.0         13) partition/G                          (S)
#  7) LUMI/24.03          (S)  14) rocm/6.2.2
$(info ########################################################################)
$(info Target lumi-hip-6.2.2: remember to )
$(info   module load LUMI/24.03 partition/G  craype-accel-amd-gfx90a PrgEnv-amd rocm/6.2.2 )
$(info ########################################################################)


### Define compiler and options

# Define compiler - use cray CC wrapper
CC := hipcc
LD := hipcc

# Define compilation flags
#CXXFLAGS  := -Ofast -flto -x c++ --std=c++17 -fno-rtti
#CXXFLAGS := -g -x c++ --std=c++17
# CXXFLAGS := -std=c++17 -fno-rtti --rocm-path=${ROCM_PATH} --offload-arch=gfx908 -x hip -fgpu-rdc
CXXFLAGS := -x hip --hip-link --std=c++20 --stdlib=libc++ -fno-rtti --rocm-path=${ROCM_PATH} -O3 -fgpu-rdc --offload-arch=gfx90a -D__HIP_PLATFORM_AMD__=1
CXXFLAGS += -D__HIP_PLATFORM_AMD__=1 
CXXFLAGS += -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 
GCC_CXX_INCLUDES := \
	/usr/include/c++/13 \
	/usr/include/c++/13/x86_64-suse-linux \
	/usr/include/c++/13/backward \
	/usr/lib64/gcc/x86_64-suse-linux/13/include \
	/usr/local/include \
	/usr/lib64/gcc/x86_64-suse-linux/13/include-fixed \
	/usr/lib64/gcc/x86_64-suse-linux/13/../../../../x86_64-suse-linux/include 

# Prepend as “-isystem” so these are searched before cuda_wrappers’s own
# built-in headers
CXXFLAGS += $(foreach dir,$(GCC_CXX_INCLUDES),-isystem $(dir))
# CXXFLAGS := -std=c++17 --offload-arch=gfx908 -x c++
#
# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

# when  module load LUMI/24.03 partition/C cpeGNU/24.03 are loaded, the compiler wrapper CC
# is inface Clang 17.0.1, whcih know the pathes of Clang c++ 17 std include files.
# Then the good practice to get HILAPP_INCLUDE_LIST is call CC -xc++ --std=c++17 -Wp,-v - otherthan us g++

# HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))
#HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | CC -xc++ --std=c++20 --stdlib=libstdc++ -Wp,-v - 2>&1 | grep "^ "))

# stddef.h again!
# HILAPP_INCLUDE_LIST += -I/opt/cray/pe/gcc/default/snos/lib/gcc/x86_64-suse-linux/default/include -I/opt/cray/pe/fftw/default/x86_64/include
#HILAPP_INCLUDE_LIST += -I/opt/cray/pe/fftw/default/x86_64/include

# Write hilapp inlcudes to a file 0hilapp_incl_dirs
#$(shell mkdir -p build)
#$(shell echo "$(HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
#HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`
HIPCC_LINK_FLAGS    := --hip-link -fgpu-rdc \
                       --offload-arch=gfx90a -D__HIP_PLATFORM_AMD__=1 \
                       --rocm-path=${ROCM_PATH}
HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o

# ROCM_LIBS := $(shell echo ${ROCM_PATH} | sed s/rocm/rocmlibs/)

# Currently LUMI EAP things require explicit setting of many include paths -- NOTE: MPI may be bad
HILA_INCLUDES := -I${ROCM_PATH}/hip/include -I${ROCM_PATH}/hip/include/hip
HILA_INCLUDES += -I${ROCM_PATH}/rocrand/include/ -I${ROCM_PATH}/hiprand/include/ 
HILA_INCLUDES += -I${ROCM_PATH}/hipfft/include/ -I${ROCM_PATH}/hipcub/include/hipcub
HILA_INCLUDES += -I${MPICH_DIR}/include
#HILA_INCLUDES += -I/opt/cray/pe/cce/17.0.1/cce/x86_64/include/craylibs
# HILA_INCLUDES += -I${MPICH_DIR}/include


# Linker libraries and possible options

# LDLIBS  := -lfftw3 -lfftw3f -lm

LDFLAGS := --hip-link -fgpu-rdc --offload-arch=gfx90a -D__HIP_PLATFORM_AMD__=1 --rocm-path=${ROCM_PATH}
LDFLAGS +=  -L${MPICH_DIR}/lib -L${CRAY_MPICH_ROOTDIR}/gtl/lib -L${ROCM_PATH}/hipfft/lib  -L${ROCM_PATH}/lib 

# libraries flags of amdhip, rocm-fft, mpich
LDLIBS := -lamdhip64 -lhipfft -lmpi -lmpi_gtl_hsa

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP
HILA_OPTS := -DHIP -DGPU_VECTOR_REDUCTION_THREAD_BLOCKS=64 -DGPU_RNG_THREAD_BLOCKS=64 $(HILA_INCLUDES)