# Benchmark nearest-neighbour kernels with/without LAYOUT_ON_THE_FLY
HILA_DIR := ../..
.DEFAULT_GOAL := bench_layout
APP_HEADERS :=
ARCH := vanilla
NDIM := 3
APP_OPTS := -DNDIM=$(NDIM)
include $(HILA_DIR)/libraries/main.mk
bench_layout: build/bench_layout ; @:
build/bench_layout: Makefile build/bench_layout.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_layout.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
