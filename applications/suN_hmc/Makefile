
# Following line(s) are printed with "make help".  Use columns 8 and 30
#%     make [suN_hmc]        - build suN_hmc program
#% suN_hmc make options 
#%     NCOL=<N>             - SU(N) gauge simulation program (default: 3)
#%     HMCS=<M>             - Gauge action for HMC:  M=0: "WILSON"
#%                                                   M=1: "Bulk-Prevention" (default)
#%                                                   M=2: "Luscher-Weisz"
#%                                                   M=3: "IWASAKI"
#%                                                   M=4: "DBW2"
#%     GFLOWS=<K>            - Gauge action for gradient flow:
#%                                                   K=0: "WILSON"
#%                                                   K=1: "Bulk-Prevention"
#%                                                   K=2: "Luscher-Weisz"   (default)
#%                                                   K=3: "IWASAKI"
#%                                                   K=4: "DBW2"

# Give the location of the top level distribution directory wrt. this location.
# Can be absolute or relative

HILA_DIR := ../..

# Number of colors
NCOL := 3

# Gauge action for HMC update: 
# "WILSON": 0, "BP": 1, "LW": 2, "IWASAKI": 3, "DBW2": 4
HMCS := 1

# Gauge action for gradient flow update: 
# "WILSON": 0, "BP": 1, "LW": 2, "IWASAKI": 3, "DBW2": 4
GFLOWS := 3

# Set default goal and arch
.DEFAULT_GOAL := suN_hmc

ifndef ARCH
ARCH := vanilla
endif

# Read in the main makefile contents, incl. platforms
include $(HILA_DIR)/libraries/main.mk

#HILAPP_OPTS += -insert-includes
#HILAPP_OPTS += -comment-pragmas
HILAPP_OPTS += -check-init

APP_OPTS += -DNDIM=4 -DNCOLOR=${NCOL} -DHMCACTION=${HMCS} -DGFLOWACTION=${GFLOWS}

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir

suN_hmc: build/su${NCOL}_hmc_hmcs${HMCS}_fs${GFLOWS}_${ARCH} ; @:

# Now the linking step for each target executable
build/su${NCOL}_hmc_hmcs${HMCS}_fs${GFLOWS}_${ARCH}: Makefile build/suN_hmc.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/suN_hmc.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS) 








