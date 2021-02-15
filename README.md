
# Description 

Hila ("lattice" in Finnish) is a C++ lattice field theory programming framework, aimed at HPC simulations.  

Purpose: make writing applications straightforward and intuitive, while producing optimized executables for 
different (super)computing platforms (parallelisation with MPI, GPU computing with Cuda or HIP, AVX vectorization, 
etc.).  Details of the parallelisation and computing architecture are hidden from the application layer.
Write once -- run anywhere.

Hila is based on hila preprocessor "hilapp", which is a C++ source-to-source transformer using the 
[libtooling](https://clang.llvm.org/docs/LibTooling.html) toolbox of the
[Clang](https://clang.llvm.org/) compiler.
It converts application C++ to platform-specific C++ code,
which is passed to appropriate compilers for the platforms.


## Quick start

- Clone hila repository (TODO: new repo address?)

~~~ bash
       git clone git@bitbucket.org:Kari_Rummukainen/hila.git
~~~

- For building *hilapp*, you need [clang](https://clang.llvm.org/) development tools (actually, only include
  files).
  These can be found in most Linux distribution repos, e.g. in Ubuntu 20.04:

~~~ bash
       apt install clang-11 llvm-11 clang-tools-11 libclang-common-11-dev libclang-cpp11-dev libclang-11-dev clang-format-11
~~~

    Change version number as needed; at least 8 required.  (TODO: what is needed for Macs?)

- Compile *hilapp*:

~~~ bash
        cd hila/hilapp
        make [-j4]
        make install
~~~

   This builds *hilapp* in hila/hilapp/build, and `make install` moves it to hila/hilapp/bin, which is the
   default location for the program.  Build takes 1-2 min.
   By default, hilapp Makefile uses clang++ installed in stage 1. You can also use g++ with `make CXX=g++`. 

- *NOTE: clang dev libraries are not installed in most supercomputer systems.  However, if the system has x86_64 
  processors (by far most common), you can use `make static` -command to build statically linked hilapp. 
  Copy `hila/hilapp/build/hilapp` to directory `hila/hilapp/bin` on the target machine.*
  

- Test `bin/hilapp -help`

- Build an application:
~~~ bash
       cd ../applications/hila_example
       make
       build/hila_example  or  mpirun -np 4 build/hila_example
~~~

- Computing platform is chosen by `make ARCH=<platform>`:
    - `make [ ARCH=vanilla ]` (often default) builds a standard MPI-parallelized program.
    - `make ARCH=AVX2` builds AVX-optimized program using [*vectorclass*](https://github.com/vectorclass) 
       library.
    - `make ARCH=cuda` builds parallel Cuda-program.  Requires nvcc compiler.

   Typically these need to be customized for supercomputing platforms.  See directory 
   hila/libraries/target_arch

- Linking *hilapp* statically: if the target machine does not have clang dev libraries, hilapp can be linked
  statically on a workstation/laptop.  Use comm
    - Move to `hila/hilapp` -directory on the workstation
    

# Overview

## A simple hila application

~~~ C++
#include "hila.h"
static_assert(NDIM == 3, "NDIM must be 3");

int main(int argc, char * argv[]) {

    hila::initialize(argc,argv);

    // set up 32^3 lattice
    lattice->setup({32,32,32});

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

    output0 << "Average of g is " << ave/lattice->volume() << '\n';

    // make a clean exit
    hila::finishrun();    
}

~~~

## Datatypes

- NDIM: number of dimensions, values 2,3,4  (TODO: NDIM=1?).  Typically set in application Makefile

- Standard types: `int`, `int64_t`, `float`, `double` (`long double`?)

- Hila provided basic types: `Complex<S>`, `Vector<n,T>`, `Matrix<n,m,T>`, `SquareMatrix<n,T>`, `Array<n,m,T>`, `Array1d<n,T>`

  Here S is any standard type, and T includes S and Complex<S>.  C++ or C standard complex types should not be used (do not
  AVX vectorize).  These types have a large number of useful methods, see (TODO Doxygen docs)

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
            output0 << v.dot({1,2,3,4});  // dot product of 2 vectors, prints -1
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
Access operation `f[X]` can be applied only to Field variables, and has the type of the
Field element (in the case above `mytype`).

`X` has methods:

- `CoordinateVector X.coordinates()`: CoordinateVector of the current site

- `int X.coordinate(Direction)`: coordinate to direction

- `Parity X.parity()`: parity of current site


The assignment `f[ALL] = 2 + g[X];` can also be done with `f = 2 + g`.
The main difference is in sequencing: the first form goes through the lattice sites in one *site loop*,
whereas the second stores the result of 2 + g to a temporary Field variable which is copied to f (in this case
std::moved).  The site loop form is faster since it minimizes temporaries and memory accesses.  

Because `f[X]` is of type Field element, the methods defined for the element type can be used.  
`f[X].dagger()` is ok, `f.dagger()` is not.

`f[X]` also serves as a visual identifier for a Field variable access.

Reduction:

~~~ C++
    mytype d = 0;
    onsites(ALL) d += f[X] - g[X+e_x];

    output0 << "The reduction is << d << std::endl;
~~~

Other features:

~~~ C++
    double a = 3, b = 5;
    Field<double> f, g=0;

    onsites(ALL) {
        f[X] = (a + b);            // ok, loop extern variables a,b do not change within the loop
        b = f[X];                  // ERROR: cannot change a loop extern non-Field variable (except reductions)
        double c = sin(f[X]);      // ok, variable c defined within the loop
        f[X] = c + g;              // ERROR: using Field variable g without [X]
    }

    CoordinateVector v = {0,1,1,0};

    f = g.shift(v);                // these two
    f[ALL] = g[X + v];             // are equivalent

    f[EVEN] = g[X + v];            // Cannot be done with g.shift() alone
~~~




# Instructions

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




## Using the Makefile system



## Syntax - What works

### Single line statements

# Goals

 1. Write tests for existing and new features
     * Test also things that fail. The test basically defines how things should work.
 2. Extend both to support:
     * Fourier transform of field variable
     * If (or where) statements in a loop
     * Reduction on dimension (so f(t) = sum_x g(x,t))
     * Array-of-Struct-of-Arrays layout
     * MPI
 3. Implement OpenACC once the compiler version is updated
 4. Get rid of NDIM
 5. Extend field to allow lattice as input
 6. Document the library
