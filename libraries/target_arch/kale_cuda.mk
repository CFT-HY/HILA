# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

$(info Remember to: module load GCC/8.2.0-2.31.1 OpenMPI/3.1.3-GCC-8.2.0-2.31.1 CUDA/10.1.105)

# Define compiler
CC = nvcc
LD = nvcc 


CUDAVER = -gencode arch=compute_70,code=sm_70

# Define compilation flags
CXXFLAGS = -dc --x cu $(CUDAVER) -std c++14 $(MPI_INCLUDE_DIRS)

#CXXFLAGS = -g -x c++ --std=c++17 

STD_INCLUDE_DIRS :=  $(addprefix -I, $(shell echo | mpic++ -xc++ --std=c++17 -Wp,-v - 2>&1 | grep "^ /"))

# STD_INCLUDE_DIRS := -I/cvmfs/fgi.csc.fi/apps/el7/aalto/spack/software/gcc/9.2.0/sepgfqb/lib/gcc/x86_64-pc-linux-gnu/9.2.0/../../../../include/c++/9.2.0

# Need to give include directory to mpi for hilapp - here 2 common ones
# MPI_INCLUDE_DIRS :=  $(addprefix -I, $(shell mpic++ -showme:incdirs))

MPI_INCLUDE_DIRS := -I/cvmfs/fgi.csc.fi/apps/el7/aalto/spack/software/openmpi/3.1.4/l567div/include

# MPI_LIBS := $(addprefix -L, $(shell mpic++ -showme:libdirs)) -lmpi
# MPI_LIBS := -L/appl/opt/OpenMPI/4.0.3-GCC-9.3.0/lib -lmpi

LDLIBS = -lm -lmpi -lcuda -lcudart -lcufft -lstdc++ ${CUDAVER} # $(MPI_LIBS)


# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o build/memory_pool2.o

################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA -DCUDA -DUSE_MPI $(STD_INCLUDE_DIRS) $(MPI_INCLUDE_DIRS)
HILA_OPTS = -DUSE_MPI -DCUDA 


