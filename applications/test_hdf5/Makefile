HILA_DIR := ../..
.DEFAULT_GOAL := test_hdf5
ARCH := vanilla
NDIM := 3
APP_OPTS := -DNDIM=$(NDIM)
include $(HILA_DIR)/libraries/main.mk
test_hdf5: build/test_hdf5 ; @:
build/test_hdf5: Makefile build/test_hdf5.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_hdf5.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
