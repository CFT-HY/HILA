# Test deep halo gathers (HALO_DEPTH > 1): Field::gather_depth() and
# site_layout neighbour_offset()
#
#%     make [..] NDIM=<d>          - Dimensions of system, default 3
#%     make [..] HALO_DEPTH=<n>    - halo depth, default 3 here
#%     make [..] ANTIPERIODIC=1    - also test antiperiodic b.c. (SPECIAL_BOUNDARY_CONDITIONS)

HILA_DIR := ../..

.DEFAULT_GOAL := test_halo_depth

APP_HEADERS :=

ARCH := vanilla

NDIM := 3

HALO_DEPTH ?= 3

APP_OPTS := -DNDIM=$(NDIM)

ifdef ANTIPERIODIC
ifneq ($(ANTIPERIODIC),0)
APP_OPTS += -DSPECIAL_BOUNDARY_CONDITIONS
endif
endif

include $(HILA_DIR)/libraries/main.mk

test_halo_depth: build/test_halo_depth ; @:

build/test_halo_depth: Makefile build/test_halo_depth.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_halo_depth.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
