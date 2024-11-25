# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
# Before compiling, load gcc and cuda modules:
# module load gcc/9.1.0 cuda/11.1.0 openmpi/4.0.5-cuda
#

$(info ########################################################################)
$(info Target puhti_cuda:  remember to )
$(info   module load gcc/11.3.0 openmpi/4.1.4-cuda cuda fftw )
$(info ########################################################################)

# Define compiler
CC = nvcc
LD = nvcc -gencode arch=compute_70,code=sm_70 --use_fast_math --restrict

# Define compilation flags
CXXFLAGS = -dc -O3 -std=c++17 -x cu -gencode arch=compute_70,code=sm_70 --use_fast_math --restrict 
# 3162 is a warning about ignored inline in __global__ functions - it's not really ignored by nvcc,
# it allows definition of a function in multiple compilation units as required by c++ standard!!  
# Quiet it.
# Warning 177 is about unused variables 
CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=3162"

#CXXFLAGS = -g -x c++ --std=c++17 

# hilapp needs to know where c++ system include files are located.  This is not a problem if
# hilapp was built from system installed clang, but if hilapp was statically compiled elsewhere
# and copied here it must be told.  Instead of hunting the directories by hand, we can ask
# system installed compilers.  g++ should be present almost everywhere.  The strange incantation
# below makes g++ list the search directories.  The result is written to build/0hilapp_incl_dirs

STD_HILAPP_INCLUDE_LIST := $(addprefix -I, $(shell echo | g++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ "))
$(shell mkdir -p build)
$(shell echo "$(STD_HILAPP_INCLUDE_LIST)" > build/0hilapp_incl_dirs )
STD_HILAPP_INCLUDES := `cat build/0hilapp_incl_dirs`


MPI_LIBS =  -lmpi

LDLIBS = -lcufft -lm $(MPI_LIBS)

# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o build/memory_pool.o

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA $(STD_HILAPP_INCLUDES) 
HILA_OPTS = -DCUDA -DPUHTI


