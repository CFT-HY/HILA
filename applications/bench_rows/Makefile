# Prototype: row-structured (nested) site loops in the lexicographic layout,
# compared with hilapp-generated onsites() loops.  Needs EVEN_SITES_FIRST=0.
#
#%     make [..] NDIM=<d>          - Dimensions of system, default 3

HILA_DIR := ../..

.DEFAULT_GOAL := bench_rows

APP_HEADERS :=

ARCH := vanilla

NDIM := 3

EVEN_SITES_FIRST := 0

APP_OPTS := -DNDIM=$(NDIM)

include $(HILA_DIR)/libraries/main.mk

bench_rows: build/bench_rows ; @:

build/bench_rows: Makefile build/bench_rows.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/bench_rows.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
