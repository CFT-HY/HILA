Creating your first HILA application
==========

This section goes over how to get started with HILA and creating your first HILA application

Like most c++ applications, HILA applications require two things, a makefile and application source code. Due to the functionality that HILA offers, the makefile and source code follow a well defined structure. Generally HILA applications are at their core c++ and the user is free to implement any methods and libraries they see fit. But to implement the functionality that the pre processor offers, a well defined default application structure is introduced:

~~~
applications/
├── hila_example
│   ├── build
│   │   ├── foo.cpt
│   │   ├── foo.o
|   |   .
|   |   ├── hila_simple_example.cpt
|   |   ├── hila_simple_example.o
|   |   .
│   │   ├── bar.cpt
│   │   └── bar.o
│   ├── Makefile
│   ├── parameters
│   └── src
│       ├── hila_example.cpp
│       └── hila_simple_example.cpp
.
.
.
~~~

In the structure HILA offers a directly, aptly named `applications`, where one can create their respective application directories. In here we have created a application directory `hila_example` which we will highlight in this section. Inside the `hila_example` directory we have the following necessary parts, visible in the above directory tree. 

The `build` directory is the location to which the `.o` object files are dumped into. The object files are compiled from the `.cpt` files — these will be discussed later in the documentation — created by the hila preprocessor. The resulting executable after compilation is also compiled into this directory.

The `Makefile` is self evidently the necessary makefile used by make to compile the hila application.

The `parameters` file is an optional file to define application parameters into. This is not necessary for the use of HILA applications, but is quite useful. This will be discussed later

Lastly the `src` directory is the directory where the user will define their HILA applications. In here we have two example HILA applications of which we will highlight `hila_simple_example.cpp`.

This file structure is necessary for the use of the makefile which handles the linking of HILA libraries used by the user.

## Table of contents

1. [Makefile system](#makefile-system)
2. [Simple HILA application](#simple-hila-application)
3. [Conclusion](#conclusion)

## Makefile system {#makefile-system}

Each application requires a makefile to link the necessary HILA libraries and to allow specification of the target backend. An application makefile should define any target files and include the main makefile defined for the HILA libraries. The main makefile handles the HILA library linking and inclusion of the target backend.

The following makefile handles the compilation of two seperate hila example applications `hila_example.cpp` and `hila_simple_example.cpp`:

~~~makefile
# Give the location of the top level distribution directory wrt. this location. Can be absolute or relative 
HILA_DIR := ../..

# A useful definition is to set the default target backend to be used for computing. In our example we set the default target backend to vanilla, which is the pure CPU MPI implementation. This allows one to skip the need of defining ARCH in the make process `make ARCH=vanilla -> make`.
ifndef ARCH
ARCH := vanilla
endif

# We then include the default makefile for hila applications which handles all the nitty gritty of defining paths for the target architecture and linking all the necessary libraries. This make file also handles use of the hila preprocessor:
include $(HILA_DIR)/libraries/main.mk

# One can also define options for the HILA preprocessor in this makefile by appending to the environment variable HILAPP_OPTS. In the example code with add the `-check-init` flag, but for now we will not explain what it's use is. We will discuss all the hila preprocessor flags later in the documentation.
HILAPP_OPTS += -check-init

# Additionally one can add HILA application options in the makefile. For example we set the system dimensions by appending to the `APP_OPTS` environment variable.
APP_OPTS += -DNDIM=3

# With multiple targets we want to use "make target", not "make build/target". This is needed to carry the dependencies to build-subdir
hila_example: build/hila_example ; @:
hila_simple_example: build/hila_simple_example ; @:

# Now the linking step for each target executable
build/hila_example: Makefile build/hila_example.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/hila_example.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

build/hila_simple_example: Makefile build/hila_simple_example.o $(HILA_OBJECTS) $(HEADERS)
	$(LD) -o $@ build/hila_simple_example.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
~~~

The only point of note is the definition for the respective object file locations with ` build/hila_simple_example.o`. For an applications this needs to be formatted in the same way as above, otherwise linking of c++ libraries and HILA objects will not be done correctly. The general format would be:

    build/{own application srouce name}: Makefile build/{own application srouce name}.o $(HILA_OBJECTS) $(HEADERS)
        $(LD) -o $@ build/{own application srouce name}.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)

### Target backends

The target backends are defined in the folder HILA/libraries/target_arch. There are two types of target backends. General ones defined for specific paralellization technologies:

| ARCH=   | Description                                                                                                            |
|---------|------------------------------------------------------------------------------------------------------------------------|
| `vanilla` | default CPU implementation with MPI                                                                                           |
| `AVX2   ` | AVX vectorization optimized program using [*vectorclass*](https://github.com/vectorclass)                              |
| `openmp ` | OpenMP parallelized program                                                                                            |
| `cuda   ` | Parallel [CUDA](https://developer.nvidia.com/cuda-toolkit) program                                                     |
| `hip    ` | Parallel [HIP](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3/page/Introduction_to_HIP_Programming_Guide.html) |

And ones which are defined for specific HPC platforms:

| ARCH=       | Description                                               |
|------------|-----------------------------------------------------------|
| `lumi      ` | CPU-MPI implementation for LUMI supercomputer             |
| `lumi-hip  ` | GPU-MPI implementation for LUMI supercomputer using HIP   |
| `mahti     ` | CPU-MPI implementation for MAHTI supercomputer            |
| `mahti-cuda` | GPU-MPI implementation for MAHTI supercomputer using CUDA |

The latter definitions are due to the module systems and non-standard paths defined by supercomputing platforms.

## Simple hila application {#simple-hila-application}

Now that we have discussed the appropriate makefile we can move on to a simple HILA application.

We offer a simple HILA application `hila_simple_example.cpp` which computes a random gaussian field (f), its laplacian (g) and the average of the laplacian field is given by:

~~~cpp
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
    onsites(ALL) f[X].gaussian_random();

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

Like all c++ applications, our program starts withing the main function. Before it, we need to include some necessary header files. At the beginning of the file we include the `hila.h` header file which contains all of the definitions for HILA libraries. This is necessary to gain access to HILA functionality. Additionally we use a static_assert to test our defined application option `-DNDIM=3`. This is useful redundancy so that we do not compile our application incorrectly.
~~~cpp
#include "hila.h"
static_assert(NDIM == 3, "NDIM must be 3");
~~~

After this process we call a few initialization and setup functions with the following lines of code:

~~~cpp
hila::initialize(argc,argv);

// set up 32^3 lattice
lattice.setup({32,32,32});

// Random numbers are used here
hila::seed_random(32345);
~~~

The first command `hila::initialize(argc,argv)` handles the initialization of mpi and reading in command line arguments and parameter files. This is a vital part of all HILA applications, but it is not too important for the user to understand what happens within it. 

Next we setup the lattice and it's size with the command `lattice.setup({32,32,32})`. The `lattice` object is defined globally within hila and contains all the information on how the lattice is split within MPI. As with initialization, this is also a vital part of any HILA application, but is designed in a way where the user need not worry about it. Note that due to NDIM option above passing for example `{32,32}` to `lattice.setup()` would result in a runtime error. TODO: CATCH THIS ERROR

Lastly for setup we initialize the random number generator with the command `hila::seed_random(32345)`. This will initialise the random number generator with the seed 32345.

Next in the application we define two Fields with:

~~~cpp
Field<Complex<double>> f;
Field<double> g = 0;
~~~

A Field in HILA is the numerical object which we operate on and iterate over. The size and MPI layout of the Fields are inherited from the lattice structure which was initialized before hand with the `lattice.setup()` command. Field is a c++ object which can be of many different data types and operations between them have been defined within HILA. This is implemented using standard c++ object oriented programming where we define the type within the brackets <T>. The available datatypes will be thoroughly documented later. For now we define one field of type `Complex<double>` and `double`. The latter Field g is initialized with the = constructor, where we set the Field to be uniformly 0. The f Field is initialized to null. 

We then introduce our first onsites loop which set's a complex gaussian random number for each point within the field:

~~~cpp
onsites(ALL) f[X].gaussian_random();
~~~

In essence this is the most important functionality that HILA offers. Onsites loops allow the user to very simply loop over the whole field without having to think about indexing, memory alignment, communication or any of the complications that writing c++ and MPI brings about. Essentially these loops are glorified for loops. With the HILA pre processor the above onsites loop expands to the following c++ code:

<details markdown="1">
<summary>onsites expansion</summary>

~~~cpp
// make f Gaussian random distributed
     //--  onsites(ALL) f[X].gaussian_random()
    {
      hila::check_that_rng_is_initialized();
      Field<Complex<double>> & _HILA_field_f = f;
      _HILA_field_f.check_alloc();
      const lattice_struct & loop_lattice = lattice;
      const int loop_begin = loop_lattice.loop_begin(Parity::all);
      const int loop_end   = loop_lattice.loop_end(Parity::all);
      for(int _HILA_index = loop_begin; _HILA_index < loop_end; ++_HILA_index) {
        Complex<double> _HILA_field_f_at_X;
        // Initial value of variable _HILA_field_f_at_X not needed
        _HILA_field_f_at_X.gaussian_random();
        _HILA_field_f.set_value_at(_HILA_field_f_at_X, _HILA_index);
      }
      hila::set_allreduce(true);
      _HILA_field_f.mark_changed(Parity::all);
    }
    //----------
~~~
AHHH SCARY PUT IT AWAY!!!
</details>

As we can see the expansion is complicated and scary, one can imagine how complicated it get's with different computing platforms. The X variable withing the onsites loop is a reserved variable within HILA applications. This variable is what defines the index of every point within the field. Appropriately the command `f[X].gaussian_random()` defines a gaussian random number for each point X within the field f. The ALL parameter within the onsites loop defines that we will iterate throughout the whole field. We will discuss variability of this parameter later in the documentation.

Next we compute:

 \f{align}{g(X) &= |\nabla^2 f(X)| \\
         &= \sum_{d \in \hat{e}} |f(X + d) - 2f(X) + f(X-d)|, \f}
 
where \f$\hat{e} = \{e_x,e_y,e_z\}\f$ is the set of unit vectors that allow us to iterate over all directions. In HILA to iterate over all directions we use the foralldir pragma. The resulting HILA code is:

```cpp
foralldir(d) {
    g[ALL] += abs(f[X+d] - 2*f[X] + f[X-d]);
}
```

We use a sum reduction assignment operator withing the foralldir pragma to indicate the sum in the laplacian equation. With assignment operators we can use the ALL variable directly to index the field f, which is equivalent to writing:

```cpp
foralldir(d) {
    onsites(ALL) g[X] += abs(f[X+d] - 2*f[X] + f[X-d]);
}
```

We then compute the average of this previously computed norm of the Laplacian of the field f using a similar sum reduction with the assignment operator:

```cpp
double average = 0;
onsites(ALL) {
    average += g[X];
}
```

We compute the average of each point with respect to the size of the system, which is give by `lattice.volume()`, since the lattice holds all the information of the systems structure. To output this value we use the default stream for text which limits the output only to the root node, so that we do not duplicate output from all mpi ranks. This default stream is held within the `hila::out0` command contained in the `hila` namespace. It is of type `std::ostream` meaning that it is essentially an alias to `std::cout` of the zeroth node.

```cpp
average = average/lattice.volume()
hila::out0 << "Average of g is " << average << '\n';
```

Lastly we wrap up the HILA application with the `hila::finishrun` command which cleans up mpi and performs a safe exit with a memory cleanup step. Additionally it prints out useful timing information coupled with a timestamp. Like `hila::initialize`, `lattice.setup` and `hila::seed_random`, this is a necessary method to call in any HILA application, especially when running with MPI.

## Conclusion {#conclusion}

This concludes the section on creating your first HILA application. We have gone through the basic structure of HILA applications, and how they are built and compiled. Additionally we have discussed basic functionality that HILA offers. With this foundational knowledge one can move on to reading the comprehensive guide on [HILA functionality](./hila_functionality.md)