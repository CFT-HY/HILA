Snippets of important information.
=========

Old guide. Saving for if there is some useful guidance here or important informtion to add to the new guide.

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

You can refer to neighbouring sites by adding a Direction (`e_x`, `-e_x`, `e_y`, `-e_y`, `e_z`, `-e_z`, `e_t`, `-e_t`, ...):
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

As before, you can refer to neighbouring sites by adding a Direction:
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
field for each Direction. It allows us to refer to the gauge field as ```gauge_field<SUN> U```
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

Operators are classes that define an `apply(Field<type> input, Field<type> output)` method.
The method takes the a field and runs a transformation on it, returning the result in
the output field.

The [Wilson Dirac](@ref Dirac_Wilson) and [staggered Dirac](@ref dirac_staggered) operators
are defined in libraries/dirac. They implement the two most common lattice Dirac operators.
These files also have the even-odd preconditioned versions of these operators.

The Dirac operators also have a `dagger(Field<type> input, Field<type> output)` method, which
implements the conjugate of the operator.

The [conjugate gradient](@ref CG) operator calculates the inverse of the square of an operator
applied to a vector
\f[
    out = CG(D) in = \frac{1}{D^\dagger D} in.
\f]
In is defined in libraries/dirac/conjugate_gradient.h

Note that the [Hasenbusch preconditioned operator](@ref Hasenbusch_operator) in
libraries/dirac/conjugate_gradient.h is a utility class used in the Hasenbusch action.


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

# OLD HILA README

# Overview

## A simple hila application

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

    // calculate sum of 2nd derivatives of f to g
    foralldir(d) {
        g[ALL] += abs(f[X+d] - 2*f[X] + f[X-d]);
    }

    // get average of g
    double ave = 0;
    onsites(ALL) {
        ave += g[X];
    }

    hila::out0 << "Average of g is " << ave/lattice.volume() << '\n';

    // make a clean exit
    hila::finishrun();    
}

~~~
You can compile this at `hila/applications/hila_example/` with `make simple` and run it with `build/simple`

## Datatypes

- NDIM: number of dimensions, values 2,3,4  (TODO: NDIM=1?).  Typically set in application Makefile

- Standard types: `int`, `int64_t`, `float`, `double` (`long double`?)

- Hila provided basic types: `Complex\<S\>`, `Vector\<n,T\>`, `Matrix\<n,m,T\>`, `SquareMatrix\<n,T\>`, `Array\<n,m,T\>`

  Here S is any standard type, and T includes S and Complex\<S\>.  C++ or C standard complex types should not be used (not
  AVX vectorizable).   See docs for functions/methods (TODO Doxygen docs)

- Special types: 
    - `Parity`: enum with values EVEN, ODD, ALL; refers to parity of the site.
       Parity of site (x,y,z,t) is even if `(x+y+z+t)` is even, odd otherwise.

    - `Direction`: conceptually unit vector with values `±e_x, ±e_y, ±e_z, ±e_t`  (if NDIM==4).
       Implemented as an enum class.  Can be used to index arrays of size NDIM.

    - `CoordinateVector`: derived from Vector<NDIM,int>. 
            
       Direction variable acts as an unit vector in vector algebra:
       (assume below NDIM=4)

~~~ C++            
     CoordinateVector v;
     Direction d = e_x;
     v = d;             // v = [1,0,0,0]
     v += e_y - 3*d;    // v = [-2,1,0,0]
     v = {0,1,-1,0};    // v = [0,1,-1,0] equivalent v = e_y - e_z;
     hila::out0 << v.dot({1,2,3,4});  // dot product of 2 vectors, prints -1
     int j = d;         // ok
     d = j;             // ERROR: cannot assign int to Direction
     ++d;               // e_x -> e_y
     is_up_dir(d);      // true if d is along positive x,y,z,t -dir.
     
~~~            

## Field access and traversal

The principal traversal of the lattice is with *site loops* `onsites(Parity)`, and a
special location identifier `X` (effectively a new keyword).  

~~~ C++
   using mytype = Matrix<3,3,Complex<double>>;   // use type alias
   Field<mytype> f,g,h;
   . . .

   onsites(ALL) f[X] = 2 + g[X];          // 2 acts as 2*I for square matrices
   f[ALL] = 2 + g[X];                     // equivalent shorter form for simple 1-line assignments
   f = 2 + g;                             // this is also equivalent!

   parity p = EVEN;
   Direction d = e_x;

   onsites(p) {
       auto t = g[X + d];                 // X +- Direction fetches from neighbour site
       f[X] += t + t*t;                   // can define variables in the loop   

       h[X] = g[X + e_x - 2*e_y];         // non-nearest neighbour fetch (TODO:optimize!)

       if (X.coordinate(e_t) == 0) {      // Do this on 1st timeslice only
            h[X] *= 0.5;
       }
   }
   
~~~

`X` can be used only inside site loops.  
Access operation `f[X]` can be applied only to field variables, and has the type of the
field element (in the case above `mytype`).

`X` has methods:

- `CoordinateVector X.coordinates()`: CoordinateVector of the current site

- `int X.coordinate(Direction)`: coordinate to direction

- `Parity X.parity()`: parity of current site


The assignment `f[ALL] = 2 + g[X];` can also be done with `f = 2 + g`.
The main difference is in sequencing: the first form goes through the lattice sites in one *site loop*,
whereas the second stores the result of 2 + g to a temporary field variable which is copied to f (in this case
std::moved).  The site loop form is faster since it minimizes temporaries and memory accesses.  

Because `f[X]` is of type field element, the methods defined for the element type can be used.  
`f[X].dagger()` is ok, `f.dagger()` is not.

`f[X]` also serves as a visual identifier for a field variable access.

Reduction:

~~~ C++
    mytype d = 0;
    onsites(ALL) d += f[X] - g[X+e_x];

    hila::out0 << "The reduction is << d << std::endl;
~~~

Other features:

~~~ C++
    double a = 3, b = 5;
    Field<double> f, g=0;

    onsites(ALL) {
        f[X] = (a + b);            // ok, loop extern variables a,b do not change within the loop
        b = f[X];                  // ERROR: cannot change a loop extern non-field variable (except reductions)
        double c = sin(f[X]);      // ok, variable c defined within the loop
        f[X] = c + g;              // ERROR: using field variable g without [X]
    }

    CoordinateVector v = {0,1,1,0};

    f = g.shift(v);                // these two
    f[ALL] = g[X + v];             // are equivalent

    f[EVEN] = g[X + v];            // Cannot be done with g.shift() alone
~~~

Access field at a single point: `f[CoordinateVector]`.  This can be used only outside site loops.

~~~ C++
  CoordinateVector v = {2,3,4,5};
  auto value = f[v];              // "value" is broadcast to all nodes!
  f[v] = 1;                       // In assignment, values are not broadcast: the node which
                                  // owns site v must have correct rhs.
~~~


## Input library

Class hila::input can be used to read parameters and other data for simulation programs.
It matches key-value pairs from input files.  As an example, if the file `parameters.dat` contains

```
    # this is a comment
    # Run parameters for run XYZ

    lattice size  64, 64, 64, 128
    beta          5.4
    clover        perturbative

    loops         25000
    seed          3474212

    coefficients    0.5, 0.7, 0.85, 1.3, 1.6, 2
    labels        setA, setB, setC
```

it can be read (mostly) using the method `input::get(std::string key)`:

~~~ C++
#include "hila.h"

int main(int argc, char * argv[]) {

    hila::initialize(argc,argv);

    // open file after hila::initialize
    // here to variable p
    hila::input p("parameters.dat");

    CoordinateVector lsize = p.get("lattice size");
    double beta            = p.get("beta");

    // Calling get_item() as below means that allowed values for 
    // "clover" are:  "tree", "perturbative", or a float/double value.
    // Return value is 0, 1, 2 respectively.
    int i = p.get_item("clover",{"tree","perturbative","%f"});
    double clover;
    if (i == 0) 
        clover = 1;
    else if (i == 1) 
        clover = <perturbative expression>;
    else 
        clover = p.get();  // the number is read here without key argument

    int loops       = p.get("loops");
    long rng_seed   = p.get("seed");

    // reading a std::vector<> reads in comma-separated values
    // this reads in a vector of 6 doubles
    std::vector<double> run_coefficients = p.get("coefficients");

    // and this a vector of 3 strings
    std::vector<std::string> labels      = p.get("labels");

    // Close the file. File is also closed when p gets out of scope
    p.close();   

    // lattice setup is convenient to do after parameters have been read
    lattice.setup(lsize);

~~~

- The method `input::get()` above deduces the type to be 
  read in from the expected return value.  The order is fixed, the items (lines)
  cannot be swapped (TODO: should this be allowed?).
  If an error occurs (wrong keys or values), program exits with an error message.

- Because the order is fixed, the keys don't really carry information for the program.  However, they
  help to ensure that the values are as intended.

- The method `input::get()` broadcasts the values to all nodes.  They have to be called by
  all nodes simultaneously.

- Method `input::get_value()` has more options for synchronization and error returns.  See 
  documentation in `input.h`

## Check input and layout

The input files and the lattice layout can be checked with the 
commands (after the application program has been built)
~~~ bash
   <hila-program-name> check
   <hila-program-name> check=<number-of-nodes>        # without spaces
~~~
This runs the program without initializing MPI, Cuda or other hardware features and
exits at `lattice.setup()` before any large memory allocations are made.  If the 
number-of-nodes argument is given, program reports how the node layout is done.

Example: if you built the `hila_example` program above, in directory `hila/applications/hila_example`
the command `build/hila_example check=32` checks the input file and the layout to 32 nodes.





# Stale Instructions

## Generating documentation

Build the documentation (with the git hash as the version number) using
~~~ bash
PROJECT_NUMBER=$(git rev-parse --short HEAD) doxygen
~~~

## Compiling the preprocessing tool and using it on c++ code

In short, the framework can be used in these steps: 

1. Write c++ code using the syntax and datatypes laid out below
2. Use the hilapp excecutable to convert this code into .cpt code 
3. Compile the new .cpt code into the final excecutable

![Workflow illustration](/docs/workflowV1.png)
 

You can then use it to compile an extended C++ file into standard C++ using
~~~ bash
bin/hilapp path/to/program.cpp
~~~
This will create a `cpt` file written in standard C++.

The `cpt` can be compiled with any c++ compiler, but must be linked against the headers and c++ files in the plumbing directory.

Check the example programs in the programs folder. You can use almost any standard C++ code, by there are a couple of new reserved names: the variable `X` and the function `onsites()`. In addition the framework defines a global `lattice` variable, which you should not overwrite.

In order to use the additional features for field type variables, you should inlude `plumbing/field.h` in you program. You can also include one or more of the files in the `datatypes` folder, which contains predefined datatypes that can be used to construct a field.

