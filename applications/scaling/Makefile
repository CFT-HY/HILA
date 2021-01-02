# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

APP_OPTS := -DNDIM=2


# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

onsim: build/onsim ; @:

 
# Now the linking step for each target executable
build/onsim: Makefile build/onsim.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/onsim.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

