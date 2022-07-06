HILA_DIR := ..

ifndef ARCH
ARCH := vanilla
endif

.DEFAULT_GOAL := catch_main

include $(HILA_DIR)/libraries/main.mk

APP_OPTS += -DNDIM=3

TEST_OBJECTS = \
	build/test_field.o \
	build/test_coordinates.o\
	build/test_real_var_ops.o\
	build/test_array.o\
	build/test_cmplx.o

HILA_OBJECTS += $(TEST_OBJECTS)
CXXFLAGS += -g

catch_main: build/catch_main ; @:

build/%: build/%.o Makefile $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $< $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

