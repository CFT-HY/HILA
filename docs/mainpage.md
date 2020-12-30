
# Description 

The Hila framework consists of 

1. the hilapp and 
2. a lattice simulation library.

The library is currently found in `programs/plumbing`.

Transformer contains a C++ preprocessing tool and framework for programming lattice field theory simulations: the main method for getting measurements from non-perturbative quantum field theories.  

Lattice field theory simulations involve up to 4 dimensional grids whose points are updated continuously according to some Monte Carlo update algorithm. The update at each grid point depends on the data stored at neighboring lattice points, and the ideal 
update order depends on the parity of the lattice points. Efficient parallel implementations using MPI and GPU's can be quite complicated to implement, modify and debug. Some kind of simplification was clearly needed. 

Transformer aims to make it easier for researchers to implement a broad class of these simulations by abstracting a lot of the technicalities involved, and by bringing the syntax from CUDA kernels and MPI calls to the essentials. The approach given 
here involves new datatypes and a preprocessing tool that converts c++ code with the new simplified syntax for loops and element accessors into working c++ code and gpu kernels. 

# HilaPP

## Generating this documentation

Build the documentation (with the git hash as the version number) using
~~~ bash
PROJECT_NUMBER=$(git rev-parse --short HEAD) doxygen
~~~

## Compiling the HilaPP and using it on c++ code

In short, the framework can be used in these steps: 

1. Write c++ code using the syntax and datatypes laid out below
2. Use the hilapp excecutable to convert this code into .cpt code 
3. Compile the new .cpt code into the final excecutable

\anchor imgDef
![imgDef]

[imgDef]: workflowV1.png "Workflow illustration"


To compile the hilapp, first create a build directory inside the main directory if it doesn't exist. 
Then, compile the hilapp by typing `make` in the main folder.
This will create an executable called `hilapp` in the build folder.

You can then use it to compile an extended C++ file into standard C++ using
~~~ bash
build/hilapp path/to/program.cpp
~~~
This will create a `cpt` file written in standard C++.

The `cpt` can be compiled with any c++ compiler, but must be linked against the headers and c++ files in the plumbing directory.

Check the example programs in the programs folder. You can use almost any standard C++ code, by there are a couple of new reserved names: the variable `X` and the function `onsites()`. In addition the framework defines a global `lattice` variable, which you should not overwrite.

In order to use the additional features for field type variables, you should inlude `plumbing/field.h` in you program. You can also include one or more of the files in the `datatypes` folder, which contains predefined datatypes that can be used to construct a field.

After this you need to define an output stream. To use `stdout`, add the line
~~~ C++
std::ostream &hila::output = std::cout;
~~~
Next, we need to initialize a lattice. The following constructs a lattice and sets it as default
~~~ C++
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;
~~~


## Using the Makefile system

Each of the example applications has a makefile for compiling the application with a
given target backend. To compile it, run
~~~ bash
make TARGET=target program_name
~~~
The lower case target should be replaced by one of `vanilla`, `AVX` or `CUDA`. This
will create a `build` directory and compile the application there.

An application makefile should define any target files and include the main makefile.
Here is an example with comments:
~~~
# Give the location of the top level distribution directory wrt. this.
# Can be absolute or relative
TOP_DIR := ../..

# Add an application specific header to the dependencies
APP_HEADERS := application.h

# Read in the main makefile contents, incl. platforms
include $(TOP_DIR)/libraries/main.mk

# With multiple targets we want to use "make target", not "make build/target".
# This is needed to carry the dependencies to build-subdir
application: build/application ; @:

# Now the linking step for each target executable
build/application: Makefile build/application.o $(HILA_OBJECTS) $(HEADERS) 
	$(LD) -o $@ build/application.o $(HILA_OBJECTS) $(LDFLAGS) $(LDLIBS)
~~~



### Compiling on Puhti

There is a separate makefile for compiling hilapp on Puhti.
To use it, run
~~~
module load gcc
make -f Makefile_puhti
~~~

This will link against the llvm installation in the hila development project folder.


## Syntax - What works

### Single line statements

You can operate on fields using statements like
~~~ C++
my_field[ALL] = my_other_field[X] + my_third_field[X];
~~~
On the left-hand side of the statement you should specify
either `[ALL]` lattice sites, `[EVEN]` sites or `[ODD]` sites.
The statement will apply only to this collection of sites.
On the right hand side, use `[X]` to refer to this collection
of sites.

You can refer to neighbouring sites by adding a direction (`e_x`, `-e_x`, `e_y`, `-e_y`, `e_z`, `-e_z`, `e_t`, `-e_t`, ...):
~~~ C++
my_field[EVEN] = my_field[X+e_y];
~~~

You can also operate on fields directly,
~~~ C++
my_field = my_field + 1;
~~~
This will operate on all sites and is equivalent to 
~~~ C++
my_field[ALL] = my_field[X] + 1;
~~~


### General loops 
Loops over all sites or a parity:
~~~ C++
forsites(ALL){}
forsites(EVEN){}
forsites(ODD){}
~~~
Inside the loop, refer to the sites using X:
~~~ C++
forsites(ALL){
    my_field[X] = 1;
}
~~~

As before, you can refer to neighbouring sites by adding a direction:
~~~ C++
forsites(EVEN){
    my_field[X] = my_field[X+e_y];
}
~~~



## What doesn't work (as expected)

Functions that implicitly depend on the site and return a number. For example
~~~ C++
forsites(EVEN){
    matrix_field[X].gaussian();
}
~~~
runs incorrectly with AVX. It does not actually run once for each site, but only
once for each vector.


# Extensions

## HMC

### Gauge field

The [gauge field](@ref gauge_field_base) class is mainly a convenient wrapper containing a matrix
field for each direction. It allows us to refer to the gauge field as ```gauge_field<SUN> U```
rather than ```SUN U[NDIM]```, which is inconvenient to pass as a reference.

The [gauge fields](@ref gauge_field_base) also contains a momentum field. Since fields are only
allocated if necessary, this is not a large amount of data.
[Fundamental gauge fields](@ref gauge_field) can also store a copy of the gauge field for
HMC.

### Actions

[Actions](@ref action_base) represent terms in the full action of a quantum field theory and
are used to set up the HMC simulation. Each [action](@ref action_base) implements at least

 1. `double action()`: returns the current value of the action
 2. `void force_step()`: calculates the derivative and updates canonical momenta
 3. `void draw_gaussian_fields()`: draws random values for any gaussian fields
 4. `void backup_fields()`: make a backup of gauge fields at the beginning of HMC
 5. `void restore_backup()`: restore the original field from the backup if the update
    is rejected

For example, the [gauge action](@ref gauge_action) represents the Wilson plaquette action
\f[
    S_{gauge} = \sum_x \beta\left [ 1 - Re \sum_{\mu>\nu} U_{x,\mu} U_{x+\mu,\nu} U^\dagger_{x+\nu, \mu} U^\dagger_{x,\nu} \right ]
\f]

The [fermion action](@ref fermion_action) represents the pseudo fermion action
\f[
    S_{fermion} = e^{-\sum_{x,y} \chi_x^\dagger \left(\frac{1}{D^\dagger D}\right)_{x,y} \chi_y}.
\f]
The Dirac operator can be any of the implemented [Wilson Dirac](@ref Dirac_Wilson) operator,
the [even-odd preconditioned Wilson Dirac](@ref Dirac_Wilson_evenodd) operator,
the [staggered Dirac](@ref dirac_staggered) operator or
the [even-odd preconditioned staggered Dirac](@ref dirac_staggered_evenodd) operator.
See operators below for more detail about how these and the matrix inversion are used.

At small mass it is often more efficient to split the fermion determinant
\f[
    Z_{fermion} = \int d\chi e^{-S_{fermion}} = det\left( D^\dagger D \right)
\f]
to
\f[
    det\left( D^\dagger D \right) = det\left( (D + m)^\dagger (D + m) \right)
    det\left( \frac{ D^\dagger D }{(D + m)^\dagger (D + m)} \right)
\f]
To use this, you need two actions, [Hasenbusch action 1](@ref Hasenbusch_action_1)
and [Hasenbusch action 2](@ref Hasenbusch_action_2).


### Integrators

An [integrator](@ref integrator_base) updates the gauge fields and their canonical
momenta keeping the action approximately constant. Two integrators are defined,
the [leapfrog](@ref leapfrog_integrator) and the [O2](@ref O2_integrator) (aka Omelyan)
integrators.

[Integrators](@ref integrator_base) are constructed from an action term and a lower level
integrator (or the momentum action on the lowest level). An integrator step updates
the gauge field keeping the action approximately constant.

[Integrators](@ref integrator_base) form a hierarchy, where lowest levels are run more often
in a trajectory.
The momentum action is also an integrator and forms the lowest level.
Generally the force of the gauge action is fast to calculate and should be added second.
The fermion action is the most expensive due to the inversion of the Dirac matrix and 
should be added on a high level.

The [leapfrog](@ref leapfrog_integrator) integrator requires a single evaluation of the
derivative of the action term and conserves the action to second order in the step size.
The [O2](@ref O2_integrator) integrator conserves the action to the third order in the,
but requires two evaluations.


### Full HMC

The full process of setting up HMC is
~~~ C++
// First define a gauge field
gauge_field<SU<N, double>> gauge;

// Let's just start from unity
gauge.set_unity();

// Set up the action of the gauge and momentum actions
gauge_action ga(gauge, beta);
gauge_momentum_action ma(gauge);

// Set up the first level in the intergator hierarchy
O2_integrator integrator_level_1(ga, ma);

// Define the Dirac operator
dirac_staggered_evenodd D(mass, gauge);

// and the fermion action
fermion_action fa(D, gauge);

// and finally the second level of the integrator
O2_integrator integrator_level_2(fsum, integrator_level_1);

// Now we can run an HMC trajectory
update_hmc(integrator_level_2, hmc_steps, traj_length);
~~~


## Operators

Operators are classes that define an `apply(field<type> input, field<type> output)` method.
The method takes the a field and runs a transformation on it, returning the result in
the output field.

The [Wilson Dirac](@ref Dirac_Wilson) and [staggered Dirac](@ref dirac_staggered) operators
are defined in libraries/dirac. They implement the two most common lattice Dirac operators.
These files also have the even-odd preconditioned versions of these operators.

The Dirac operators also have a `dagger(field<type> input, field<type> output)` method, which
implements the conjugate of the operator.

The [conjugate gradient](@ref CG) operator calculates the inverse of the square of an operator
applied to a vector
\f[
    out = CG(D) in = \frac{1}{D^\dagger D} in.
\f]
In is defined in libraries/dirac/conjugate_gradient.h

Note that the [Hasenbusch preconditioned operator](@ref Hasenbusch_operator) in
libraries/dirac/conjugate_gradient.h is a utility class used in the Hasenbusch action.

## Backends

Backends are primarily implemented in three places.
First, in HilaPP, loop generation and loop function handling code is in the files
`hilapp/src/codegen_*.cpp`.
The code generation functions are called in
[backend_handle_loop_function](@ref MyASTVisitor::backend_handle_loop_function)
and [backend_generate_code](@ref MyASTVisitor::backend_generate_code).

In order to define a new backend, you should edit the two functions above, implement the code
generation function and add any new files to `hilapp/Makefile`.

Second, in the library in the folders `libraries/plumbing/backend_*`. These implement
field storage in (usually in `field_storage_backend.h`), any other necessary top level
definitions in `defs.h` and possible an extension of the lattice class in `lattice.h`.
These are included in `libraries/plumbing/field_storage.h`, `libraries/plumbing/defs.h`
and `libraries/plumbing/lattice.h` respectively.

A new backend should implement at least the [field storage](@ref MyASTVisitor::field_storage)
class. The new file needs to be included in `libraries/plumbing/field_storage.h`.

Finally, `libraries/platforms` has a collection of makefiles, chosen by the `PLATFORM`
flag in the standard Makefile. These include combinations of a specific system and 
a backend. New backend requires a new makefile that defines the necessary flags
to produce and compile the correct code.


# Testing

In the `programs/test_cases` folder you can find a collection of simple test programs. To test whether the translations work on the cpu, type:

~~~ bash
./test.sh 
~~~

This tests the transform, compilation and run process for the test_*.cpp files for dimensions from 1 to 4, and outputs the exit status of each step. 
If you're on a machine with GPU's, you can test the GPU transformations with:

~~~ bash
./test_GPU.sh
~~~



# Goals

 1. Write tests for existing and new features
     * Test also things that fail. The test basically defines how things should work.
 2. Extend to support:
     * If (or where) statements in a loop
 3. Implement HIP backend (should be a simple replacement of CUDA syntax and functions)
 4. Multiple lattices in one program
 5. Expand and update documentation
 6. Implement an OpenACC backend

