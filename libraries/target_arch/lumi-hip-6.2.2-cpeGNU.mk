# NB: If using hipSolver with the rocm/6.2.2 module, need to first load SuiteSparse.
# Looks like only hipfft is used so I guess SuiteSparse is not needed.
$(info ########################################################################)
$(info Target lumi-hip-6.2.2-cpeGNU: remember to )
$(info   module load LUMI/24.03 partition/G cpeGNU/24.03 cray-fftw rocm/6.2.2)
$(info If using singularity container for hilapp remember to:)
$(info     export SINGULARITY_BIND="{HILA_INSTALL_PATH}")
$(info ########################################################################)

# In PrgEnv-gnu we need to use hipcc directly and set HIPCC envvars to properly capture Cray's include and lib paths
# https://docs.lumi-supercomputer.eu/development/compiling/prgenv/#compile-hip-code
CC := hipcc
LD := hipcc

HILA_CPP_STD := --std=c++20

# "Base" HILA CXX and link flags. Need -xhip because .cpt is nonstandard so the compiler doesn't know how to deal with them otherwise
CXXFLAGS := -O3 -x hip $(HILA_CPP_STD) -fno-rtti -fgpu-rdc
LDFLAGS := -fgpu-rdc

# Append LUMI specific flags to hipcc manually (see docs above).
# These should be equivalent to setting HIPCC_COMPILE_FLAGS_APPEND and HIPCC_LINK_FLAGS_APPEND
CXXFLAGS += --offload-arch=gfx90a
CXXFLAGS += $(shell CC --cray-print-opts=cflags)
LDFLAGS += --offload-arch=gfx90a

# HILAPP needs stdlib include dirs. Unfortunately compilers don't expose these directly as envvars,
# but they are printed by the preprocessor in verbose mode. So we parse them with an arcane incantation
HILAPP_INCLUDES :=  $(addprefix -I, $(shell $(CC) -E $(HILA_CPP_STD) -xc++ -Wp,-v - < /dev/null 2>&1 | grep "^ "))

$(info	Using HILAPP includes: $(HILAPP_INCLUDES))

HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o

# Library -l flags for the linker. Get these from the Cray environment, and in cpeGNU we must manually add hipfft.
LDLIBS := $(shell CC --cray-print-opts=libs) -lhipfft

# These variables must be defined here
#
HILAPP_OPTS := -target:HIP -DHIP $(HILAPP_INCLUDES)
HILA_OPTS := $(HILA_INCLUDES) -DHIP -DGPU_VECTOR_REDUCTION_THREAD_BLOCKS=64 -DGPU_RNG_THREAD_BLOCKS=64
