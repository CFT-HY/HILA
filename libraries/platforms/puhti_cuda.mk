# Platform specific makefile for cuda code.  Tuned for "puhti" computer at CSC
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#
#

# Define compiler
CC = nvcc
LD = nvcc -gencode arch=compute_70,code=sm_70 -fmad=false

# Define compilation flags
CXXFLAGS = -dc -x cu -gencode arch=compute_70,code=sm_70 -fmad=false -DCUDA -DUSE_MPI -I../../../cub/  $(CXX_INCLUDE)
#CXXFLAGS = -g -x c++ --std=c++17 

LDLIBS = -lfftw3 -lm 

# No need to give include directory to mpi for hilapp - here 2 common ones
MPI_INCLUDE_DIRS = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/openmpi/include

MPI_LIBS = -L/appl/spack/install-tree/gcc-8.3.0/hpcx-mpi-2.4.0-7gyvq3/lib -lmpi

LDLIBS = -lm $(MPI_LIBS)


################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA -DCUDA -DUSE_MPI -DPUHTI -DHILAPP $(MPI_INCLUDE_DIRS)
HILA_OPTS = -DAVX -DUSE_MPI 


# and add extra rule for building backend_cuda -files

build/%.o: $(LIBRARIES_DIR)/plumbing/backend_cuda/%.cpp $(ALL_DEPEND)
	mkdir -p build
	$(CC) $(CXXFLAGS) $(OPTS) $(HILA_OPTS) $< -o $@

