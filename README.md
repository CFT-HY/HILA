
# Description 

Hila (lattice in Finnish) is a C++ lattice field theory programming framework, aimed at HPC simulations.  

Purpose: make writing applications straightforward and intuitive, while producing optimized executables for 
different (super)computing platforms (parallelization with MPI, GPU computing with Cuda or HIP, AVX vectorization, 
etc.).  Details of the parallelization and computing architecture are hidden from the user's view, and 
all applications automatically run on present or future platform.
Write once -- run anywhere.

Hila is based on hila preprocessor "hilapp", which is a C++ source-to-source transformer using the 
[libtooling](https://clang.llvm.org/docs/LibTooling.html) toolbox of the
[Clang](https://clang.llvm.org/) compiler.
It converts application C++ to platform-specific C++ code,
which is passed to appropriate compilers for the platforms.

Behind the scenes hila takes care of MPI layout and communications.  It lays out the 
lattice fields differently for different computing platforms: 'array of structures' (standard),
'array of structures of vectors' (AVX-type), or 'structure of arrays' (GPU-type).

## Dependencies

### Hilapp

| Dependencies | Minimum Version   | Required  |
|--------------|-------------------|-----------|
| Clang        | 8 -               | Yes       |

### HILA applications

| Dependencies | Minimum Version   | Required  |
|--------------|-------------------|-----------|
| Clang / GCC  | 8 -    /  x       | Yes       |
| FFTW3        | x                 | Yes       |
| MPI          | x                 | Yes       |
| OpenMP       | x                 | No        |
| CUDA         | x                 | No        |
| HIP          | x                 | No        |

**Installing dependencies**

Installing all dependencies on ubuntu:
```bash
apt install build-essential \
            libopenmpi-dev \
            libfftw3-dev \
            libomp-dev
```

**CUDA:**

See NVIDIA drivers and CUDA documentation: https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html

**HIP:**

See ROCm and HIP documentation: https://docs.amd.com/, https://rocmdocs.amd.com/en/latest/Installation_Guide/HIP-Installation.html
## Quick start

Begin by cloning HILA repository:

``` bash
git clone https://haaaaron@bitbucket.org/Kari_Rummukainen/hila.git
```

The HILA working method is split into two parts. The first part is getting access to the hilapp preprocessor. And the second part is building simulations for targeted architectures and technologies using the hilapp preprocessor.

## 1. HILA preprocessor
The preprocessor can be accessed either by compiled from source using the clang libtooling toolbox or by using the offered [singularity container](https://bitbucket.org/haaaaron/hila-singularity/src/main/). (Future support for packaged .rpm and .deb? Maybe docker aswell for MACOSX and windows support RPM is supported by macosx)

#### **Compiling from source**
For building *hilapp*, you need [clang](https://clang.llvm.org/) development tools (actually, only include files). These can be found in most Linux distribution repos, e.g. in Ubuntu 20.04:

~~~ bash
export LLVM_VERSION=12
apt install clang-$LLVM_VERSION \
            llvm-$LLVM_VERSION \
            clang-tools-$LLVM_VERSION \
            libclang-common-$LLVM_VERSION-dev \
            libclang-cpp$LLVM_VERSION-dev \
            libclang-$LLVM_VERSION-dev \
            clang-format-$LLVM_VERSION
~~~

Compile *hilapp*:

~~~ bash
cd hila/hilapp
make [-j4]
make install
~~~

This builds *hilapp* in hila/hilapp/build, and `make install` moves it to hila/hilapp/bin, which is the default location for the program.  Build takes 1-2 min. 

("By default, hilapp Makefile uses clang++ installed in stage 1. You can also use g++ with `make CXX=g++`." Is this detail too complicated? Should just stick to clang in this part.) 

- *NOTE: clang dev libraries are not installed in most supercomputer systems.  However, if the system has x86_64 
  processors (by far most common), you can use `make static` -command to build statically linked hilapp. 
  Copy `hila/hilapp/build/hilapp` to directory `hila/hilapp/bin` on the target machine. Simpler approach for HPC platforms is use of singularity containers*
  
Test that hilapp works

    ./bin/hilapp --help

#### **Singularity container**

The singularity container offers a more packaged approach where one doesn't need to worry about clang libtoolbox support. Hence for HPC platforms where the access of such compiler libraries can be tedious one can simply opt to use the container. 

Since the HILA preprocessor doesn't require optimized performance any overhead is irrelevant. In theory there should be none.

See separate [repository](https://bitbucket.org/haaaaron/hila-singularity/src/main/) for use of singularity container.

(Add downloadable .sif file to git page?)

## 2. Building HILA applications

The second part is building HILA applications. Here we will go over a test example. All applications should lie in the applications folder.

- *NOTE: that at this point one will need to install the FFTW3 and OpenMPI development libraries, see dependencies section* 

Build an application:
``` bash
cd hila/applications/hila_example
make [-j4]
./build/hila_example
```
By default all HILA applications are built using MPI so one can run directly:

    mpirun -n 4 ./build/hila_example

 
Computing platform is chosen by 

    make ARCH=<platform>

- `[ ARCH=vanilla ]` (often default) builds a standard MPI-parallelized program
- `ARCH=AVX2` builds AVX-optimized program using [*vectorclass*](https://github.com/vectorclass)
- `ARCH=openmp` builds OpenMP parallelized program
- `ARCH=cuda` builds parallel CUDA-program
- `ARCH=hip` builds parallel HIP-program

Typically these need to be customized for supercomputing platforms due to stack limitations of said platforms. ~~See directory hila/libraries/target_arch~~ -> TODO: should have a list of all system specific target architectures

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

    output0 << "The reduction is << d << std::endl;
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
    lattice->setup(lsize);

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
exits at `lattice->setup()` before any large memory allocations are made.  If the 
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


