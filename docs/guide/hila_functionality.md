HILA Functionality
====================

This section is comprehensive description of the functionality it offers. For technical documentation each class, method, function etc. has been (work in progress) documented with standard docstring documentation which has been generated with doxygen.

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

~~~cpp       
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

~~~cpp
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

~~~cpp
    mytype d = 0;
    onsites(ALL) d += f[X] - g[X+e_x];

    hila::out0 << "The reduction is << d << std::endl;
~~~

Other features:

~~~cpp
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

~~~cpp
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

~~~cpp
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
~~~bash
   <hila-program-name> check
   <hila-program-name> check=<number-of-nodes>        # without spaces
~~~
This runs the program without initializing MPI, Cuda or other hardware features and
exits at `lattice.setup()` before any large memory allocations are made.  If the 
number-of-nodes argument is given, program reports how the node layout is done.

Example: if you built the `hila_example` program above, in directory `hila/applications/hila_example`
the command `build/hila_example check=32` checks the input file and the layout to 32 nodes.

## HILA pre-processor tool

In short, the framework can be used in these steps: 

1. Write c++ code using the syntax and datatypes laid out below
2. Use the hilapp excecutable to convert this code into .cpt code 
3. Compile the new .cpt code into the final excecutable

![Workflow illustration](../images/workflowV1.png)
 

You can then use it to compile an extended C++ file into standard C++ using
~~~bash
bin/hilapp path/to/program.cpp
~~~
This will create a `cpt` file written in standard C++.

The `cpt` can be compiled with any c++ compiler, but must be linked against the headers and c++ files in the plumbing directory.

Check the example programs in the programs folder. You can use almost any standard C++ code, by there are a couple of new reserved names: the variable `X` and the function `onsites()`. In addition the framework defines a global `lattice` variable, which you should not overwrite.

In order to use the additional features for field type variables, you should inlude `plumbing/field.h` in you program. You can also include one or more of the files in the `datatypes` folder, which contains predefined datatypes that can be used to construct a field.