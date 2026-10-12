# Test on-the-fly site layout (site_layout.h) against lattice arrays
#
#%     make [..] NDIM=<d>    - Dimensions of system, default 3

HILA_DIR := ../..

.DEFAULT_GOAL := test_site_layout

APP_HEADERS :=

ARCH := vanilla

NDIM := 3

APP_OPTS := -DNDIM=$(NDIM)

include $(HILA_DIR)/libraries/main.mk

test_site_layout: build/test_site_layout ; @:

build/test_site_layout: Makefile build/test_site_layout.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_site_layout.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
