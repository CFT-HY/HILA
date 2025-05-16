# Specialize cuda.mk to CUDA_ARCH=89
#
# this is included from main.mk -file, which is in turn included from 
# application makefile
#

CUDA_ARCH := 89
include $(HILA_DIR)/libraries/target_arch/cuda.mk
