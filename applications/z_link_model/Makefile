
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [z_link_model]        - build z_link_model program
# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

ifndef ARCH
ARCH := vanilla
endif

# Set default goal and arch
.DEFAULT_GOAL := z_link_model

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=4

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

z_link_model: build/z_link_model_${ARCH} ; @:

# Now the linking step for each target executable
build/z_link_model_${ARCH}: Makefile build/z_link_model.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/z_link_model.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 

