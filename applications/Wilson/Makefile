
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

APP_HEADERS := Wilson.h

APP_OPTS := -DNDIM=4

# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk


# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

Wilson: build/Wilson ; @:

 
# Now the linking step for each target executable
build/Wilson: Makefile build/Wilson.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/Wilson.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

