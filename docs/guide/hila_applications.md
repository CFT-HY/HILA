Building you first HILA application
==========

This section goes over how to get started with HILA and building your first HILA application

## HILA Application

Like most c++ applications, HILA applications require two thing, a makefile and application source code. Due to the functionality that HILA offers, the makefile and source code follow a well defined structure. Generally HILA applications are at their core c++ and the user is free to implement any methods and libraries they see fit. But to implement the functionality that the pre processor offers, a well defined skeleton is introduced.

### Makefile system

Each application requires a makefile to link the necessary HILA libraries and to allow specification of the target backend. An application makefile should define any target files and include the main makefile defined for the HILA libraries. The main makefile handles the HILA library linking and inclusion of the target backend.

Here is an example with comments:
~~~Make
#NECESSARY
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
# Allows the application folder to be defined anywhere in the machine
HILA_DIR := ../..

# Set default goal
.DEFAULT_GOAL := applications

#Set default architecture 
ifndef ARCH
ARCH := vanilla
endif

# Add an application specific header to the dependencies
APP_HEADERS := application.h

# Read in the main makefile contents to link and build HILA libraries
include $(TOP_DIR)/libraries/main.mk

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to build the dependencies in the build subdirectory
application: build/application ; @:

# Now the linking step for each target executable
build/application: Makefile build/application.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/application.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
~~~

TODO: **Should the above makefile illustrate which aspects are necessary and which are not**

The target backends are defined in the folder HILA/libraries/target_arch. There are two types of target backends. General ones defined for specific paralellization technologies:

| ARCH=   | Description                                                                                                            |
|---------|------------------------------------------------------------------------------------------------------------------------|
| `vanilla` | default CPU implementation                                                                                             |
| `AVX2   ` | AVX vectorization optimized program using [*vectorclass*](https://github.com/vectorclass)                              |
| `openmp ` | OpenMP parallelized program                                                                                            |
| `cuda   ` | Parallel [CUDA](https://developer.nvidia.com/cuda-toolkit) program                                                     |
| `hip    ` | Parallel [HIP](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3/page/Introduction_to_HIP_Programming_Guide.html) |

And ones which are defined for specific HPC platforms:

| ARCH       | Description                                               |
|------------|-----------------------------------------------------------|
| `lumi      ` | CPU-MPI implementation for LUMI supercomputer             |
| `lumi-hip  ` | GPU-MPI implementation for LUMI supercomputer using HIP   |
| `mahti     ` | CPU-MPI implementation for MAHTI supercomputer            |
| `mahti-cuda` | GPU-MPI implementation for MAHTI supercomputer using CUDA |

The latter definitions are due to the module systems and non-standard paths defined by supercomputing platforms.

### Simple hila application

A simple HILA application which computes a random gaussian field (f), it's derivative (g) and the average of the derivative field is given by:

~~~ C++
#include "hila.h"
static_assert(NDIM == 3, "NDIM must be 3");

int main(int argc, char * argv[]) {

    hila::initialize(argc,argv);

    // set up 32^3 lattice
    lattice.setup({32,32,32});

    // Random numbers are used here
    hila::seed_random(32345);

    Field<Complex<double>> f;
    Field<double> g = 0;

    // make f Gaussian random distributed
    onsites(ALL) f[X].gaussian();

    // calculate sum of 2nd derivatives of f in to g
    foralldir(d) {
        g[ALL] += abs(f[X+d] - 2*f[X] + f[X-d]);
    }

    // get average of g
    double average = 0;
    onsites(ALL) {
        average += g[X];
    }

    average = average/lattice.volume()
    hila::out0 << "Average of g is " << average << '\n';

    // make a clean exit
    hila::finishrun();    
}

~~~

One can compile this code at `HILA/applications/hila_example/` with `make simple` and run it with `./build/simple`.

TODO: **CONTINUE FROM HERE**