# Benchmark f[X + k*dir] access: deep halo (HALO_DEPTH > 1) against shift-copy (HALO_DEPTH=1)
#
#%     make [..] NDIM=<d>          - Dimensions of system, default 3
#%     make [..] HALO_DEPTH=<n>    - halo depth, default 3 here; 1 gives shift-copy

HILA_DIR := ../..

.DEFAULT_GOAL := bench_halo_depth

APP_HEADERS :=

ARCH := vanilla

NDIM := 3

HALO_DEPTH ?= 3

APP_OPTS := -DNDIM=$(NDIM)

include $(HILA_DIR)/libraries/main.mk

bench_halo_depth: build/bench_halo_depth ; @:

build/bench_halo_depth: Makefile build/bench_halo_depth.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_halo_depth.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
