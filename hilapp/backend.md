Backends
========

Backends are primarily implemented in three places.
First, in HilaPP, loop generation and loop function handling code is in the files
`hilapp/src/codegen_*.cpp`.
The code generation functions are called in
[backend_handle_loop_function](@ref TopLevelVisitor::backend_handle_loop_function)
and [backend_generate_code](@ref TopLevelVisitor::backend_generate_code).

In order to define a new backend, you should edit the two functions above, implement the code
generation function and add any new files to `hilapp/Makefile`.

Second, in the library in the folders `libraries/plumbing/backend_*`. These implement
field storage in (usually in `field_storage_backend.h`), any other necessary top level
definitions in `defs.h` and possible an extension of the lattice class in `lattice.h`.
These are included in `libraries/plumbing/field_storage.h`, `libraries/plumbing/defs.h`
and `libraries/plumbing/lattice.h` respectively.

A new backend should implement at least the [field storage](@ref TopLevelVisitor::field_storage)
class. The new file needs to be included in `libraries/plumbing/field_storage.h`.

Finally, `libraries/platforms` has a collection of makefiles, chosen by the `ARCH`
flag in the standard Makefile. These include combinations of a specific system and 
a backend. New backend requires a new makefile that defines the necessary flags
to produce and compile the correct code.