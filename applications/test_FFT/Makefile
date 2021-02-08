
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

.DEFAULT_GOAL := test_FFT

# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas

APP_OPTS += -DNDIM=3

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

test_FFT: build/test_FFT ; @:

# Now the linking step for each target executable
build/test_FFT: Makefile build/test_FFT.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/test_FFT.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)





