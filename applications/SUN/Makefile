# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

APP_HEADERS := SUN.h

# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

SUN: build/SUN ; @:

# Now the linking step for each target executable
build/SUN: Makefile build/SUN.o build/monte.o build/KennedyPendleton.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/SUN.o build/monte.o build/KennedyPendleton.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)


