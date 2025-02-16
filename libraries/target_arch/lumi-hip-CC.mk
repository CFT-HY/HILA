# Platform specific makefile for LUMI standard (CPU) code
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

$(info ########################################################################)
$(info Target lumi-hip-CC: remember to )
$(info   module load CrayEnv PrgEnv-cray craype-accel-amd-gfx90a cray-mpich rocm )
$(info )
$(info Arcording to note of changes of  update of August-September 2024 from LUMI team,)
$(info LUMI/24.04 is the only truly-supported software stack, you may condider use )
$(info LUMI/24.03 and its toolchains to make workflow smooth. )
$(info To load LUMI/24.03 with Cray Clang compiler, AMD GPU gfx90a libs, ROCm and toolchain run )
$(info )
$(info   module --force purge && module --force unload LUMI)
$(info   module load LUMI/24.03 partition/G cpeCray/24.03 buildtools/24.03 rocm/6.0.3 )
$(info LUMI/24.03 contains EasyBuild extra-packages building system, read LUMI documentaion to learn more.)
$(info Moreover, the compiler wrapper CC has full konwledge of C++ standard include files, which are significant for hilapp functionality.)
$(info ########################################################################)


### Define compiler and options

# Define compiler - use cray CC wrapper
CC := CC
LD := CC

# Define compilation flags
#CXXFLAGS  := -Ofast -flto -x c++ --std=c++17 -fno-rtti
#CXXFLAGS := -g -x c++ --std=c++17
# CXXFLAGS := -std=c++17 -fno-rtti --rocm-path=${ROCM_PATH} --offload-arch=gfx908 -x hip -fgpu-rdc
CXXFLAGS := -std=c++17 -fno-rtti -O3 -fgpu-rdc --offload-arch=gfx90a -D__HIP_PLATFORM_AMD__=1
CXXFLAGS += -D__HIP_PLATFORM_AMD__=1 
CXXFLAGS += -D__HIP_ROCclr__ -D__HIP_ARCH_GFX90A__=1 -xhip
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
HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | CC -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))

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
HILA_INCLUDES := -I${ROCM_PATH}/hip/include -I${ROCM_PATH}/hip/include/hip
HILA_INCLUDES += -I${ROCM_PATH}/rocrand/include/ -I${ROCM_PATH}/hiprand/include/ 
HILA_INCLUDES += -I${ROCM_PATH}/hipfft/include/ -I${ROCM_PATH}/hipcub/include/hipcub
# HILA_INCLUDES += -I${MPICH_DIR}/include


# Linker libraries and possible options

# LDLIBS  := -lfftw3 -lfftw3f -lm
LDFLAGS := --offload-arch=gfx90a -D__HIP_PLATFORM_AMD__=1 --hip-link --rocm-path=${ROCM_PATH} -L${ROCM_PATH}/lib 
LDFLAGS += -L${ROCM_PATH}/hipfft/lib
LDFLAGS += -L${MPICH_DIR}/lib -L${CRAY_MPICH_ROOTDIR}/gtl/lib

# libraries flags of amdhip, rocm-fft, mpich
LDLIBS := -lamdhip64 -lhipfft -lmpi

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP $(HILAPP_INCLUDES)
HILA_OPTS := -DHIP -DGPU_VECTOR_REDUCTION_THREAD_BLOCKS=64 -DGPU_RNG_THREAD_BLOCKS=64 $(HILA_INCLUDES)



