# This is the main makefile of hila applications
#  - to be called in application makefiles 
#  - calls platform specific makefiles in directory "platforms" 
#

# If "make clean", don't worry about targets, platforms and options
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleanall)

## hilapp binary (HILA_DIR defined in calling Makefile)

HILAPP_DIR := $(HILA_DIR)/hilapp
HILAPP := $(HILAPP_DIR)/bin/hilapp

################

LIBRARIES_DIR := $(HILA_DIR)/libraries
ARCH_DIR := $(LIBRARIES_DIR)/target_arch
HILA_INCLUDE_DIR := $(HILA_DIR)/libraries

# ARCH needs to be defined. Check.
ifndef ARCH
  $(info ########################################################################)
  $(info Target architecture (ARCH) is not defined.)
  $(info Use "make ARCH=<target-arch>" or define ARCH in application Makefile)
  $(info For available target architectures, see the directory $(ARCH_DIR) )
  $(info ########################################################################)
  $(error )
endif

.PRECIOUS: build/%.cpt build/%.o

HILA_OBJECTS = \
	build/initialize.o \
	build/input.o \
	build/mersenne_inline.o \
	build/gaussian_random.o \
	build/lattice.o \
	build/map_node_layout_trivial.o \
	build/memalloc.o \
	build/timing.o \
	build/test_gathers.o \
	build/com_mpi.o \
	build/com_single.o \
	build/fft.o

# com_mpi / com_single could be moved to platforms, but they're protected by USE_MPI guards

# Read in the appropriate platform bits and perhaps extra objects
include $(ARCH_DIR)/$(ARCH).mk

# Define LAYOUT_VECTOR if vector (SUBNODE) layout is desired
ifdef LAYOUT_VECTOR
	HILA_OBJECTS += build/setup_layout_vector.o
	HILA_OPTS += -DSUBNODE_LAYOUT
else
	HILA_OBJECTS += build/setup_layout_generic.o
endif

# To force a full remake when changing platforms or targets
CLEANED_GOALS := $(shell echo ${MAKECMDGOALS} | sed -e 's/ /_/g' -e 's/\//+/g' | cut -c1-60)
LASTMAKE := build/.lastmake.${CLEANED_GOALS}.${ARCH}

$(LASTMAKE): $(MAKEFILE_LIST)
	@mkdir -p build
	-rm -f build/.lastmake.*
	make clean
	touch ${LASTMAKE}



# Use all headers inside libraries for dependencies
HILA_HEADERS := $(wildcard $(HILA_DIR)/libraries/*/*.h) $(wildcard $(HILA_DIR)/libraries/*/*/*.h)

ALL_DEPEND := $(LASTMAKE) $(HILA_HEADERS)

HILA_OPTS += -I$(HILA_INCLUDE_DIR) -I$(HILA_INCLUDE_DIR)/plumbing

# Add the (possible) std. includes for hilapp
HILAPP_OPTS += -I$(HILAPP_DIR)/clang_include $(CUSTOM_HILAPP_OPTS)

#
#  GIT VERSION: tricks to get correct git version and build date
#  on the file

GIT_SHA := $(shell git rev-parse --short=8 HEAD)

ifneq "$(GIT_SHA)" "" 
HILA_OPTS += -DGIT_SHA_VALUE=$(GIT_SHA)
GIT_SHA_FILE := build/.git_sha_number_$(GIT_SHA)

# Force recompilation if git number has changed

$(GIT_SHA_FILE):
	-rm -f build/.git_sha_number_*
	touch $(GIT_SHA_FILE)

ALL_DEPEND += $(GIT_SHA_FILE)

endif

# Standard rules for creating and building cpt files. These
# build .o files in the build folder by first running them
# through hilapp.
# Source can be at . , ./src,  or under libraries/plumbing
# If your project has more complicated hierarcy, add the rules to the
# project Makefile

build/%.cpt: src/%.cpp Makefile $(MAKEFILE_LIST) $(ALL_DEPEND) $(APP_HEADERS)
	@mkdir -p build
	$(HILAPP) $(HILAPP_OPTS) $(APP_OPTS) $(HILA_OPTS) $< -o $@ $(HILAPP_TRAILING_OPTS)

build/%.cpt: %.cpp Makefile $(MAKEFILE_LIST) $(ALL_DEPEND) $(APP_HEADERS)
	@mkdir -p build
	$(HILAPP) $(HILAPP_OPTS) $(APP_OPTS) $(HILA_OPTS) $< -o $@ $(HILAPP_TRAILING_OPTS)

build/%.o : build/%.cpt
	$(CC) $(CXXFLAGS) $(APP_OPTS) $(HILA_OPTS) $< -c -o $@

build/%.cpt: $(LIBRARIES_DIR)/plumbing/%.cpp $(ALL_DEPEND) $(HILA_HEADERS)
	@mkdir -p build
	$(HILAPP) $(HILAPP_OPTS) $(APP_OPTS) $(HILA_OPTS) $< -o $@ $(HILAPP_TRAILING_OPTS)


# This one triggers only for cuda targets
build/%.cpt: $(LIBRARIES_DIR)/plumbing/backend_cuda/%.cpp $(ALL_DEPEND) $(HILA_HEADERS)
	@mkdir -p build
	$(HILAPP) $(HILAPP_OPTS) $(APP_OPTS) $(HILA_OPTS) $< -o $@ $(HILAPP_TRAILING_OPTS)



endif
endif   # close the "clean" bracket

.PHONY: clean cleanall

clean:
	-rm -f build/*.o build/*.cpt build/.lastmake* build/.git_sha_number*

cleanall:
	-rm -f build/*
