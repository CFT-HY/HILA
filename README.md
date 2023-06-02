
HILA
========= 
![axions](./images/AxionStringNetwork.png)

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

## User guide

Installation and user guide can be viewed on the web documentation hosted on: TODO: add link to github pages or hosted guide somewhere

To generate the user guide and technical documentation locally one can run:

    doxygen /docs/config 

To open the documentation locally with any browser:

    firefox /docs/html/index.html

## License 

The GPL-2.0 license information can be reviewed at the [license](./LICENSE) page