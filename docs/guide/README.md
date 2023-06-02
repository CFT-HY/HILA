Description and Installation  {#mainpage}
========= 

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

## Table of Contents
1. [Dependencies](#dependencies)
2. [Installation](#installation)
3. [User Guide](#user-guide)

## Dependencies {#dependencies}

### Hilapp

| Dependencies | Minimum Version   | Required  |
|--------------|-------------------|-----------|
| Clang        | 8 -               | Yes       |

#### Installing dependencies for HILA preprocessor:

If one opts to use a docker or singularity container, skip  directly to the [installation](#installation) section.

For building *hilapp*, you need [clang](https://clang.llvm.org/) development tools (actually, only include files). These can be found in most Linux distribution repos, e.g. in Ubuntu 22.04:

~~~bash
export LLVM_VERSION=15
sudo apt-get -y install clang-$LLVM_VERSION \
                   libclang-$LLVM_VERSION-dev
~~~

### HILA applications

| Dependencies | Minimum Version   | Required  |
|--------------|-------------------|-----------|
| Clang / GCC  | 8 -    /  x       | Yes       |
| FFTW3        | x                 | Yes       |
| MPI          | x                 | Yes       |
| OpenMP       | x                 | No        |
| CUDA         | x                 | No        |
| HIP          | x                 | No        |

#### Installing dependencies for HILA applications:

If one opts to use docker, skip directly to the [installation](#installation) section.

Installing all non GPU dependencies on ubuntu:
~~~bash
sudo apt install build-essential \
            libopenmpi-dev \
            libfftw3-dev \
            libomp-dev
~~~

**CUDA:**

See NVIDIA drivers and CUDA documentation: https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/index.html

**HIP:**

See ROCm and HIP documentation: https://docs.amd.com/, https://rocmdocs.amd.com/en/latest/Installation_Guide/HIP-Installation.html

## Installation{#installation}

Begin by cloning HILA repository:

~~~bash
git clone https://github.com/CFT-HY/HILA
~~~

The installation process is split into two parts. Building the HILA preprocessor and compiling HILA applications. Both can be installed from source, and both steps have their respective containerization options available. The variety in options is to address differing issues which arise in platform dependencies.

When it comes to installing HILA applications there are many avenues one can take depending on their platform. The available platforms and offered methods are listed below, which link to the necessary section in the installation guide.

#### LINUX

HILA has originally been developed on linux, hence all of the available options can be used. The HILA preprocessor can be built from [source](#hila-preprocessor) or with the use of a [singualrity](#singularity) container. Additionally one can opt to use the [docker](#docker) container which installs the hila preprocessor directly.

#### MAC

On mac the installation of the HILA preprocessor dependencies and HILA application dependencies can be tedious, and in some cases impossible. Availability of clang libtoolbox is open ended. For this reason the best option is to use the available [docker](#docker) container.

#### WINDOWS

On windows the installation of the HILA preprocessor dependencies and HILA application dependencies are untested. For this reason the best option is to use the available [docker](#docker) container. 

#### HPC

On supercomputing platforms the HILA application dependencies are most likely available. The only issue is the availability of the clang libtoolbox which is used in building the HILA preprocessor. Due to the availability of singularity on supercomputing platfroms the best solution is to opt to use the [singularity](#singularity) container. 

**After installing the hila preprocessor with one of the above options one can move on to the [building HILA applications](#building-hila-applications) section.**

### Containers

HILA comes with both a singularity and docker container for differing purposes. The aim is to make use easy on any platform be it linux, mac, windows or a supercomputer.

#### Docker {#docker}

The docker container is meant to develop and produce HILA applications, libraries and hilapp with ease. One can produce HILA applications on their local machine and run them in a container without having to worry about dependencies. Note that there is overhead when running MPI communication in docker, thus one will not get optimal simulation performance when running highly paralelled code in a container. This is a non issue with small scale simulations or testing.

For instructions on using the docker container have a look at the [README.md](../../docker/README.md) in the docker folder

#### Singularity {#singularity}

The singularity container offers a more packaged approach where one doesn't need to worry about clang libtoolbox support for compiling the HILA pre processor. Hence for HPC platforms where the access of such compiler libraries can be tedious one can simply opt to use the container version of hilapp. This approach is mainly meant to be used for pre processing applications on an HPC platform.

For instructions on installing singularity and building containers have a look at the [README.md](../../singularity/README.md) in the singularity folder

### HILA preprocessor {#hila-preprocessor}

Before building the preprocessor one must first install the dependencies. See the [dependencies](#dependencies)

Compile *hilapp*:

~~~bash
cd hila/hilapp
make [-j4]
make install
~~~

This builds *hilapp* in hila/hilapp/build, and `make install` moves it to hila/hilapp/bin, which is the default location for the program.  Build takes 1-2 min. 

("By default, hilapp Makefile uses clang++ installed in stage 1. You can also use g++ with `make CXX=g++`." Is this detail too complicated? Should just stick to clang in this part.) 

- *NOTE: clang dev libraries are not installed in most supercomputer systems.  However, if the system has x86_64 
  processors (by far most common), you can use `make static` -command to build statically linked hilapp. 
  Copy `hila/hilapp/build/hilapp` to directory `hila/hilapp/bin` on the target machine. Simpler approach for HPC platforms is use of singularity containers*
  
Test that hilapp works:

    ./bin/hilapp --help

<details>
<summary> Expected output </summary>

~~~bash
$ ./bin/hilapp --help
USAGE: hilapp [options] <source files>

OPTIONS:

Generic Options:

  --help                    - Display available options (--help-hidden for more)
  --help-list               - Display list of available options (--help-list-hidden for more)
  --version                 - Display the version of this program

hilapp:

  --AVXinfo=<int>           - AVX vectorization information level 0-2. 0 quiet, 1 not vectorizable loops, 2 all loops
  -D <macro[=value]>        - Define name/macro for preprocessor
  -I <directory>            - Directory for include file search
  --allow-func-globals      - Allow using global or extern variables in functions called from site loops.
                              This will not work in kernelized code (for example GPUs)
  --check-init              - Insert checks that Field variables are appropriately initialized before use
  --comment-pragmas         - Comment out '#pragma hila' -pragmas in output
  --dump-ast                - Dump AST tree
  --function-spec-no-inline - Do not mark generated function specializations "inline"
  --gpu-slow-reduce         - Use slow (but memory economical) reduction on gpus
  --ident-functions         - Comment function call types in output
  --insert-includes         - Insert all project #include files in .cpt -files (portable)
  --method-spec-no-inline   - Do not mark generated method specializations "inline"
  --no-include              - Do not insert any '#include'-files (for debug, may not compile)
  --no-interleave           - Do not interleave communications with computation
  --no-output               - No output file, for syntax check
  -o <filename>             - Output file (default: <file>.cpt, write to stdout: -o - 
  --syntax-only             - Same as no-output
  --target:AVX              - Generate AVX vectorized loops
  --target:AVX512           - Generate AVX512 vectorized loops
  --target:CUDA             - Generate CUDA kernels
  --target:HIP              - Generate HIP kernels
  --target:openacc          - Offload to GPU using openACC
  --target:openmp           - Hybrid OpenMP - MPI
  --target:vanilla          - Generate loops in place
  --target:vectorize=<int>  - Generate vectorized loops with given vector size 
                              For example -target:vectorize=32 is equivalent to -target:AVX
  --verbosity=<int>         - Verbosity level 0-2.  Default 0 (quiet)
~~~

</details>


### Building HILA applications {#building-hila-applications}

The second part is building HILA applications. Here we will go over an example with a health check test application. All applications should lie in the applications folder.

- *NOTE: that at this point one will need to install the FFTW3 and OpenMPI development libraries, see dependencies section* 

Build an application:
~~~bash
cd hila/applications/hila_healthcheck
make [-j4]
./build/hila_healthcheck
~~~

<details>
<summary> Expected output </summary>

~~~bash
$ ./build/hila_healthcheck 
----- HILA ⩩ lattice framework ---------------------------
Running program ./build/hila_healthcheck
with command line arguments ''
Code version: git SHA d0222bca
Compiled Jun  1 2023 at 11:13:10
with options: EVEN_SITES_FIRST SPECIAL_BOUNDARY_CONDITIONS
Starting -- date Thu Jun  1 11:13:28 2023  run time 8.328e-05s
No runtime limit given
GNU c-library performance: not returning allocated memory
----- Reading file parameters ------------------------------
lattice size         256,256,256
random seed          0
------------------------------------------------------------
------------------------------------------------------------
LAYOUT: lattice size  256 x 256 x 256  =  16777216 sites
Dividing to 1 nodes

Sites on node: 256 x 256 x 256  =  16777216
Processor layout: 1 x 1 x 1  =  1 nodes
Node remapping: NODE_LAYOUT_BLOCK with blocksize 4
Node block size 1 1 1  block division 1 1 1
------------------------------------------------------------
Communication tests done -- date Thu Jun  1 11:13:31 2023  run time 3.11s
------------------------------------------------------------
Random seed from time: 3871436182438
Using node random numbers, seed for node 0: 3871436182438
--- Complex reduction value ( -2.7647453e-17, 5.5294928e-17 ) passed
--- Vector reduction, sum ( -7.1331829e-15, -1.4328816e-15 ) passed
--- Setting and reading a value at [ 37 211 27 ] passed
--- Setting and reading a value at [ 251 220 47 ] passed
--- Setting and reading a value at [ 250 249 134 ] passed
--- Maxloc is [ 112 117 164 ] passed
--- Max value 2 passed
--- Minloc is [ 192 135 27 ] passed
--- Min value -1 passed
--- Field set_elements and get_elements with 51 coordinates passed
--- SiteSelect size 51 passed
--- SiteValueSelect size 51 passed
--- SiteSelect content passed
--- SiteValueSelect content passed
--- SiteIndex passed
--- 2-dimensional slice size 65536 passed
--- slice content passed
--- 1-dimensional slice size 256 passed
--- slice content passed
--- FFT constant field passed
--- FFT inverse transform passed
--- FFT of wave vector [ 132 159 243 ] passed
--- FFT of wave vector [ 167 161 208 ] passed
--- FFT of wave vector [ 152 87 255 ] passed
--- FFT of wave vector [ 156 86 229 ] passed
--- FFT of wave vector [ 78 246 141 ] passed
--- FFT real to complex passed
--- FFT complex to real passed
--- Norm of field = 44434.862 and FFT = 44434.862 passed
--- Norm of binned FFT = 44434.862 passed
--- Binning test at vector [ 100 220 7 ] passed
--- Spectral density test with above vector  passed
--- Binning test at vector [ 193 10 49 ] passed
--- Spectral density test with above vector  passed
--- Binning test at vector [ 235 241 96 ] passed
--- Spectral density test with above vector  passed
TIMER REPORT:             total(sec)          calls     time/call  fraction
---------------------------------------------------------------------------
MPI broadcast       :          0.000             40      0.263 μs   0.0000
MPI reduction       :          0.000             34      2.003 μs   0.0000
FFT total time      :         44.544             14      3.182 s    0.6449
 copy pencils       :          3.261             15      0.217 s    0.0472
 MPI for pencils    :          0.000             90      1.298 μs   0.0000
 FFT plan           :          0.003             42     73.150 μs   0.0000
 copy fft buffers   :          2.412        5505024      0.438 μs   0.0349
 FFT execute        :          2.356        2752512      0.856 μs   0.0341
 pencil reshuffle   :         12.967             30      0.432 s    0.1878
 save pencils       :         26.043             15      1.736 s    0.3771
bin field time      :          9.014              7      1.288 s    0.1305
---------------------------------------------------------------------------
 No communications done from node 0
Finishing -- date Thu Jun  1 11:14:37 2023  run time 69.07s
------------------------------------------------------------
~~~

**NOTE: Naturally the run time depends on your system**

</details>

By default all HILA applications are built using MPI so one can run:

    mpirun -n 4 ./build/hila_healthcheck

<details>
<summary> Expected output </summary>

~~~bash
$ mpirun -n 4 ./build/hila_healthcheck
----- HILA ⩩ lattice framework ---------------------------
Running program ./build/hila_healthcheck
with command line arguments ''
Code version: git SHA d0222bca
Compiled Jun  1 2023 at 11:13:10
with options: EVEN_SITES_FIRST SPECIAL_BOUNDARY_CONDITIONS
Starting -- date Thu Jun  1 11:18:22 2023  run time 0.0001745s
No runtime limit given
GNU c-library performance: not returning allocated memory
----- Reading file parameters ------------------------------
lattice size         256,256,256
random seed          0
------------------------------------------------------------
------------------------------------------------------------
LAYOUT: lattice size  256 x 256 x 256  =  16777216 sites
Dividing to 4 nodes

Sites on node: 256 x 128 x 128  =  4194304
Processor layout: 1 x 2 x 2  =  4 nodes
Node remapping: NODE_LAYOUT_BLOCK with blocksize 4
Node block size 1 2 2  block division 1 1 1
------------------------------------------------------------
Communication tests done -- date Thu Jun  1 11:18:23 2023  run time 1.046s
------------------------------------------------------------
Random seed from time: 4184648360436
Using node random numbers, seed for node 0: 4184648360436
--- Complex reduction value ( -2.7539926e-17, 5.5079939e-17 ) passed
--- Vector reduction, sum ( 1.4328816e-15, -7.4627804e-15 ) passed
--- Setting and reading a value at [ 139 215 41 ] passed
--- Setting and reading a value at [ 231 44 102 ] passed
--- Setting and reading a value at [ 238 201 150 ] passed
--- Maxloc is [ 80 69 74 ] passed
--- Max value 2 passed
--- Minloc is [ 219 105 178 ] passed
--- Min value -1 passed
--- Field set_elements and get_elements with 51 coordinates passed
--- SiteSelect size 51 passed
--- SiteValueSelect size 51 passed
--- SiteSelect content passed
--- SiteValueSelect content passed
--- SiteIndex passed
--- 2-dimensional slice size 65536 passed
--- slice content passed
--- 1-dimensional slice size 256 passed
--- slice content passed
--- FFT constant field passed
--- FFT inverse transform passed
--- FFT of wave vector [ 239 139 86 ] passed
--- FFT of wave vector [ 218 12 247 ] passed
--- FFT of wave vector [ 94 206 99 ] passed
--- FFT of wave vector [ 34 78 96 ] passed
--- FFT of wave vector [ 221 224 199 ] passed
--- FFT real to complex passed
--- FFT complex to real passed
--- Norm of field = 44418.915 and FFT = 44418.915 passed
--- Norm of binned FFT = 44418.915 passed
--- Binning test at vector [ 106 69 123 ] passed
--- Spectral density test with above vector  passed
--- Binning test at vector [ 240 142 174 ] passed
--- Spectral density test with above vector  passed
--- Binning test at vector [ 226 28 118 ] passed
--- Spectral density test with above vector  passed
TIMER REPORT:             total(sec)          calls     time/call  fraction
---------------------------------------------------------------------------
MPI broadcast       :          0.002             40     49.358 μs   0.0001
MPI reduction       :          0.289             34      8.508 ms   0.0120
MPI post receive    :          0.000              4      1.782 μs   0.0000
MPI start send      :          0.000              4      3.923 μs   0.0000
MPI wait receive    :          0.001              4      0.277 ms   0.0000
MPI wait send       :          0.002              4      0.404 ms   0.0001
MPI send field      :          0.001             15     67.812 μs   0.0000
FFT total time      :         14.922             14      1.066 s    0.6182
 copy pencils       :          1.941             15      0.129 s    0.0804
 MPI for pencils    :          1.644             90     18.263 ms   0.0681
 FFT plan           :          0.006             42      0.140 ms   0.0002
 copy fft buffers   :          1.164        1376256      0.846 μs   0.0482
 FFT execute        :          0.933         688128      1.355 μs   0.0386
 pencil reshuffle   :          7.246             30      0.242 s    0.3002
 save pencils       :          2.994             15      0.200 s    0.1240
bin field time      :          2.792              7      0.399 s    0.1157
---------------------------------------------------------------------------
 COMMS from node 0: 4 done, 0(0%) optimized away
Finishing -- date Thu Jun  1 11:18:46 2023  run time 24.14s
------------------------------------------------------------
~~~
**NOTE: Naturally the run time depends on your system**

</details>

 
Computing platform is chosen by 

    make ARCH=<platform>

**List of computing platforms:**

| ARCH=   | Description                                                                                                            |
|---------|------------------------------------------------------------------------------------------------------------------------|
| `vanilla` | default CPU implementation                                                                                             |
| `AVX2` | AVX vectorization optimized program using [*vectorclass*](https://github.com/vectorclass)                              |
| `openmp ` | OpenMP parallelized program                                                                                            |
| `cuda` | Parallel [CUDA](https://developer.nvidia.com/cuda-toolkit) program                                                     |
| `hip` | Parallel [HIP](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3/page/Introduction_to_HIP_Programming_Guide.html) |

For cuda compilation one needs to define their CUDA version and architercure either as environment variables or during the make process:

~~~bash
export CUDA_VERSION=11.6
export CUDA_ARCH=61
make ARCH=cuda
or
make ARCH=cuda CUDA_VERSION=11.6 CUDA_ARCH=61
~~~
*NOTE: Default cuda version is 11.6 and compute architecture is sm_61*

**HPC platforms**:

| ARCH       | Description                                               |
|------------|-----------------------------------------------------------|
| `lumi` | CPU-MPI implementation for LUMI supercomputer             |
| `lumi-hip` | GPU-MPI implementation for LUMI supercomputer using HIP   |
| `mahti` | CPU-MPI implementation for MAHTI supercomputer            |
| `mahti-cuda` | GPU-MPI implementation for MAHTI supercomputer using CUDA |

## User guide {#user-guide}

Now that HILA has been built successfully, the next step is to build your first HILA application: [hila application guide](./hila_applications.md)

After building your first HILA application one can move on to the comprehensive guide, which describes everything that HILA has to offer: [comprehensive guide](./hila_functionality.md)

Both of these resources can be viewed on the web guide hosted on: TODO: add link to github pages or hosted guide somewhere

To generate the user guide and technical documentation locally one can run:

    doxygen /docs/config 

To open the documentation locally with any browser:

    firefox /docs/html/index.html