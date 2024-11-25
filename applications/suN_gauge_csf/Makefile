
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [suN_gauge_csf]        - build suN_gauge_csf program
#% suN_gauge_csf make options 
#%     NCOL=<N>              - SU(N) gauge simulation program (default: 3)
#%     SUN_OVERRELAX_dFJ=1   - use deForcrand-Jahn full overrelax (default: su2 subgroups)
#%     GFLOWS=<K>            - Gauge action for gradient flow:
#%                                                   K=0: "WILSON"
#%                                                   K=1: "Bulk-Prevention"
#%                                                   K=2: "Luscher-Weisz"
#%                                                   K=3: "IWASAKI"         (default)
#%                                                   K=4: "DBW2"

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative
HILA_DIR := ../..

# Number of colors
NCOL := 3

# Gauge action for gradient flow update: 
# "WILSON": 0, "BP": 1, "LW": 2, "IWASAKI": 3, "DBW2": 4
GFLOWS := 3
CSFLOWS := 3


ifndef ARCH
ARCH := vanilla
endif

# Set default goal and arch
.DEFAULT_GOAL := suN_gauge_csf

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=4 -DNCOLOR=${NCOL} -DGFLOWACTION=${GFLOWS} -DCSFLOWACTION=${CSFLOWS}

ifdef SUN_OVERRELAX_dFJ
APP_OPTS += -DSUN_OVERRELAX_dFJ
endif

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

suN_gauge_csf: build/su${NCOL}_gauge_csf_fs${GFLOWS}_csfs${CSFLOWS}_${ARCH} ; @:

# Now the linking step for each target executable
build/su${NCOL}_gauge_csf_fs${GFLOWS}_csfs${CSFLOWS}_${ARCH}: Makefile build/suN_gauge_csf.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/suN_gauge_csf.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 








