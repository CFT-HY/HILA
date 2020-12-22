# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

APP_HEADERS := staggered.h

# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

staggered: build/staggered ; @:

 
# Now the linking step for each target executable
build/staggered: Makefile build/staggered.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/staggered.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

