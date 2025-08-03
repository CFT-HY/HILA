
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [z_link_model]        - build z_link_model program
#% suN_gauge_gf make options 
#%     PARITY=<K>              - K=0: use 4D parity for staggering
#%                             - K=1: use spatial parity for staggering
#%     PLAQ_SHIFT=<S>          - S=1: site-parity-dependent plaquette shifts on spatial plaquettes 
#%                             - S=2: site-parity-dependent plaquette shifts on temporal plaquetts
#%                             - any other value for S: all plaquette shifts are set to zero
#&                           

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Number of colors
PARITY := 1
PLAQ_SHIFT := 1


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

APP_OPTS += -DNDIM=4 -DPARITY=${PARITY} -DPLAQ_SHIFT=${PLAQ_SHIFT} -DNODE_LAYOUT_TRIVIAL

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

z_link_model: build/z_link_model_par${PARITY}_ps${PLAQ_SHIFT}_${ARCH} ; @:

# Now the linking step for each target executable
build/z_link_model_par${PARITY}_ps${PLAQ_SHIFT}_${ARCH}: Makefile build/z_link_model.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/z_link_model.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 

