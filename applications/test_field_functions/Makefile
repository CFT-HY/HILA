# Test of field functions: functions which access field neighbours f[X+d],
# called from site loops.  See src/test_field_functions.cpp
#
#   make ARCH=<arch>

HILA_DIR := ../..

.DEFAULT_GOAL := test_field_functions

ARCH := vanilla

APP_OPTS := -DNDIM=3

include $(HILA_DIR)/libraries/main.mk

APP_OBJECTS := build/test_field_functions.o

test_field_functions: build/test_field_functions ; @:

build/test_field_functions: Makefile $(APP_OBJECTS) $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ $(APP_OBJECTS) $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
