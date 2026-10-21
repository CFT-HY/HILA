# Test of hila::global variables on GPUs, see src/test_gpu_globals.cpp
#
#   make ARCH=<arch>              user globals in __device__ memory (default)
#   make ARCH=<arch> CONSTANT=1   user globals in __constant__ memory
#
# "make clean" when switching between the two.

HILA_DIR := ../..

.DEFAULT_GOAL := test_gpu_globals

ARCH := vanilla

APP_OPTS := -DNDIM=3
ifeq ($(CONSTANT),1)
APP_OPTS += -DGPU_GLOBALS_IN_CONSTANT_MEMORY=1
endif

include $(HILA_DIR)/libraries/main.mk

APP_OBJECTS := build/test_gpu_globals.o build/file_a.o build/file_b.o

test_gpu_globals: build/test_gpu_globals ; @:

build/test_gpu_globals: Makefile $(APP_OBJECTS) $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $(APP_OBJECTS) $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
