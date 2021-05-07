# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
HILA_DIR := ../..

APP_OPTS := -DNDIM=2


# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

xy: build/xy ; @:

 
# Now the linking step for each target executable
build/xy: Makefile build/xy.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/xy.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

