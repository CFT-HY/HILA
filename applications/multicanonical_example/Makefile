# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

#SRC_DIR := ./src

# Set default goal and arch
.DEFAULT_GOAL := muca_example

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=3
APP_OBJECTS = build/multicanonical.o

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

muca_example: build/muca_example; @:

# Now the linking step for each target executable
build/muca_example: Makefile build/multicanonical_example.o $(APP_OBJECTS) $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/multicanonical_example.o $(APP_OBJECTS) $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
