# Platform specific makefile for cuda code.  Tuned for "Roihu" computer at CSC
#
# this is included from main.mk -file, which is in turn included from
# application makefile

$(info ########################################################################)
$(info Target roihu-gpu: Remember to load modules:)
$(info module load gcc/15.2.0 openmpi/5.0.8 cuda/13.1.1        # vanilla / overlap)
$(info module load gcc/15.2.0 openmpi/5.0.8 cuda/13.1.1 nccl   # + GPU_CCL)
$(info GPU_SHMEM uses nvhpc's bundled nvcc/HPCX/NVSHMEM via the paths below -- no extra module)
$(info ########################################################################)

# ---------------------------------------------------------------------------
# Communication-backend toolchain selection.
#
# Driven by what the application requested in APP_OPTS (-DGPU_CCL / -DGPU_SHMEM).
# Everything below uses deferred $(if $(findstring ...,$(APP_OPTS)) ...) so it is
# evaluated when a recipe runs -- by then the application Makefile has fully
# populated APP_OPTS. Works whether the backend is set in the app Makefile
# (APP_OPTS += -DGPU_CCL) or on the command line (OPTS=-DGPU_CCL).
#
#   vanilla / overlap : system nvcc + openmpi   (nothing extra)
#   GPU_CCL           : + -lnccl                (needs: module load nccl)
#   GPU_SHMEM         : nvhpc's bundled nvcc + HPCX MPI + NVSHMEM. NVSHMEM's
#                       libnvshmem_device.a is built against nvhpc's own CUDA, so
#                       mixing it with the system cuda module hits an nvlink
#                       device-code ABI mismatch ("ABI version 8 vs 7"); hence the
#                       whole shmem build compiles+links with nvhpc's nvcc and its
#                       ABI-matched HPCX libmpi.
# Runtime (launch wrapper, not the build) must add $(NVHPC_MPI_DIR)/lib and
# $(NVSHMEM_DIR)/lib to LD_LIBRARY_PATH, set OPAL_PREFIX=$(NVHPC_MPI_DIR) and
# NVSHMEM_SYMMETRIC_SIZE=2G -- see the communication-backend-test launch scripts.
# ---------------------------------------------------------------------------
NVHPC_ROOT    = /appl/soft/manual/general/aarch64/nvhpc/Linux_aarch64/26.3
NVHPC_NVCC    = $(NVHPC_ROOT)/cuda/13.1/bin/nvcc
NVSHMEM_DIR   = $(NVHPC_ROOT)/comm_libs/13.1/nvshmem
NVHPC_MPI_DIR = $(NVHPC_ROOT)/comm_libs/13.1/hpcx/hpcx-2.25.1/ompi

_shmem = $(findstring GPU_SHMEM,$(APP_OPTS))
_ccl   = $(findstring GPU_CCL,$(APP_OPTS))

# Define compiler -- nvhpc's bundled nvcc for GPU_SHMEM, system nvcc otherwise.
NVCC_BIN = $(if $(_shmem),$(NVHPC_NVCC),nvcc)
CC = $(NVCC_BIN)
LD = $(NVCC_BIN) -gencode arch=compute_90a,code=sm_90a --use_fast_math --restrict

# Define compilation flags  - arch 90 or 90a. 90a has GH extension features but is not backwards/forwards compatible
CXXFLAGS = -g -dc -std=c++17 -x cu -gencode arch=compute_90a,code=sm_90a --restrict
CXXFLAGS_NOOPT = $(CXXFLAGS) -O1
CXXFLAGS += -O3 --use_fast_math

# Suppress some warnings
# 3162 is a warning about ignored inline in __global__ functions - it's not really ignored by nvcc,
# it allows definition of a function in multiple compilation units as required by c++ standard!!
# Warning 177 is about unused variables
#CXXFLAGS += -Xcudafe "--display_error_number --diag_suppress=177 --diag_suppress=3162"
CXXFLAGS += -Xcudafe "--diag_suppress=177"

# GPU_SHMEM: NVSHMEM + HPCX-MPI headers (must precede the system openmpi's mpi.h)
CXXFLAGS += $(if $(_shmem),-I$(NVSHMEM_DIR)/include -I$(NVHPC_MPI_DIR)/include)

# MPI: HPCX (nvhpc-bundled, ABI-matched to nvhpc nvcc) for GPU_SHMEM, else system
# openmpi. The -L must precede -lmpi so HPCX's libmpi wins over the module's.
MPI_LIBS = $(if $(_shmem),-L$(NVHPC_MPI_DIR)/lib -lmpi,-lmpi)

LDLIBS = -lcufft -lm $(MPI_LIBS)
LDLIBS += $(if $(_ccl),-lnccl)
LDLIBS += $(if $(_shmem),-L$(NVSHMEM_DIR)/lib -lnvshmem_host -lnvshmem_device)

# extra cuda objects here
HILA_OBJECTS += build/hila_gpu.o

################

# These variables must be defined here
#
HILAPP_OPTS = -target:CUDA $(STD_HILAPP_INCLUDES)
HILA_OPTS = -DCUDA
