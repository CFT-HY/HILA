# Compiling

Compile the transformer by typing ´make´ in the main folder.
This will create an executable called ´transformer´ in the build folder.
You can compile an extended c++ file into standard c++ using
~~~
build/transformer path/to/program.cpp
~~~
This will create a cpt file written in standard C++.

This can be compiled with any c++ compiler, but must be linked against the headers and c++ files in the plumbing directory.
