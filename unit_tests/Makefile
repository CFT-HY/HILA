HILA_DIR := ..

ifndef ARCH
ARCH := vanilla
endif

.DEFAULT_GOAL := tests

include $(HILA_DIR)/libraries/main.mk

LDLIBS += -lgtest_main -lgtest -lpthread
APP_OPTS += -DNDIM=3

TEST_OBJECTS = \
	build/test1.o

HILA_OBJECTS += $(TEST_OBJECTS)

tests: build/tests ; @:
example_test: build/example_test ; @:

build/%: build/%.o Makefile $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $< $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)