HILA_DIR := ..

ifndef ARCH
ARCH := vanilla
endif

.DEFAULT_GOAL := tests

include $(HILA_DIR)/libraries/main.mk

LDLIBS += -lgtest_main -lgtest -lpthread
APP_OPTS += -DNDIM=3

TEST_OBJECTS = \
	build/test_lattice.o \
	build/test_field.o

HILA_OBJECTS += $(TEST_OBJECTS)

tests: build/tests ; @:

build/%: build/%.o Makefile $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $< $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

