# Compiling

Compile the transformer by typing `make` in the main folder.
This will create an executable called `transformer` in the build folder.
You can compile an extended c++ file into standard c++ using
~~~
build/transformer path/to/program.cpp
~~~
This will create a cpt file written in standard C++.

This can be compiled with any c++ compiler, but must be linked against the headers and c++ files in the plumbing directory.

Check the example programs in the programs folder. You can use almost any standard C++ code, by there are a couple of new reserved names: the variable `X` and the function `onsites()`. In addition the framework defines a global `lattice` variable, which you should not overwrite.

In order to use the additional features for field type variables, you should inlude `plumbing/field.h` in you program. You can also include one or more of the files in the `datatypes` folder, which contains predefined datatypes that can be used to construct a field.

After this you need to define an output stream. To use `stdout`, add the line
~~~
std::ostream &hila::output = std::cout;
~~~
Next, we need to initialize a lattice. The following constructs a lattice and sets it as default
~~~
lattice_struct my_lattice;
lattice_struct * lattice = & my_lattice;
~~~



# What works:

## Single line statements

You can operate on fields using statements like
~~~
my_field[ALL] = my_other_field[X] + my_third_field[X];
~~~
On the left-hand side of the statement you should specify
either `[ALL]` lattice sites, `[EVEN]` sites or `[ODD]` sites.
The statement will apply only to this collection of sites.
On the right hand side, use `[X]` to refer to this collection
of sites.

You can refer to neighbouring sites by adding a direction (`XUP`, `XDOWN`, `YUP`, `YDOWN`, `ZUP`, `ZDOWN`, `TUP`, `TDOWN`, ...):
~~~
my_field[EVEN] = my_field[X+YUP];
~~~


## General loops 
Loops over all sites or a parity:
~~~
forsites(ALL){}
forsites(EVEN){}
forsites(ODD){}
~~~
Inside the loop, refer to the sites using X:
~~~
forsites(ALL){
    my_field[X] = 1;
}
~~~

As before, you can refer to neighbouring sites by adding a direction:
~~~
forsites(EVEN){
    my_field[X] = my_field[X+YUP];
}
~~~


