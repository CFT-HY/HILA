
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [suN_gauge_gf]        - build suN_gauge_gf program
#% suN_gauge_gf make options 
#%     NCOL=<N>             - SU(N) gauge simulation program (default: 3)
#%     SUN_OVERRELAX_dFJ=1  - use deForcrand-Jahn full overrelax (default: su2 subgroups)

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Number of colors
NCOL := 3

# Gauge action for gradient flow update: 
# "WILSON": 0, "BP": 1, "LW": 2, "IWASAKI": 3, "DBW2": 4
GFLOWS := 3


ifndef ARCH
ARCH := vanilla
endif

# Set default goal and arch
.DEFAULT_GOAL := suN_gauge_gf

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=4 -DNCOLOR=${NCOL} -DGFLOWACTION=${GFLOWS}

ifdef SUN_OVERRELAX_dFJ
APP_OPTS += -DSUN_OVERRELAX_dFJ
endif

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

suN_gauge_gf: build/suN_gauge_gf_su${NCOL}_fs${GFLOWS}_${ARCH} ; @:

# Now the linking step for each target executable
build/suN_gauge_gf_su${NCOL}_fs${GFLOWS}_${ARCH}: Makefile build/suN_gauge_gf.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/suN_gauge_gf.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 








