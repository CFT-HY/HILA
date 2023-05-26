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

## Dependencies

### Hilapp

| Dependencies | Minimum Version   | Required  |
|--------------|-------------------|-----------|
| Clang        | 8 -               | Yes       |

#### Installing dependencies for HILA preprocessor:

If one opts to use a docker or singularity container, skip  directly to the [installation](#installation) section.

For building *hilapp*, you need [clang](https://clang.llvm.org/) development tools (actually, only include files). These can be found in most Linux distribution repos, e.g. in Ubuntu 22.04:

~~~ bash
export LLVM_VERSION=15
apt-get -y install clang-$LLVM_VERSION \
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

## Installation

Begin by cloning HILA repository:

``` bash
git clone https://github.com/CFT-HY/HILA
```

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

After installing the hila preprocessor with one of above options one can move on to the [building HILA applications](#building-hila-applications) section.

### Containers

HILA comes with both a singularity and docker container for differing purposes. The aim is to make use easy on any platform be it linux, mac, windows or a supercomputer.

#### Docker

The docker container is meant to develop and produce HILA applications, libraries and hilapp with ease. One can produce HILA applications on their local machine and run them in a container without having to worry about dependencies. Note that there is overhead when running MPI communication in docker, thus one will not get optimal simulation performance when running highly paralelled code in a container. This is a non issue with small scale simulations or testing.

For instructions on using the docker container have a look at the [README.md](docker/README.md) in the docker folder

#### Singularity

The singularity container offers a more packaged approach where one doesn't need to worry about clang libtoolbox support for compiling the HILA pre processor. Hence for HPC platforms where the access of such compiler libraries can be tedious one can simply opt to use the container version of hilapp. This approach is mainly meant to be used for pre processing applications on an HPC platform.

For instructions on installing singularity and building containers have a look at the [README.md](singularity/README.md) in the singularity folder

### HILA preprocessor

Before building the preprocessor one must first install the dependencies. See the [dependencies](#dependencies)

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


### Building HILA applications

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

**List of computing platforms:**

| ARCH=   | Description                                                                                                            |
|---------|------------------------------------------------------------------------------------------------------------------------|
| `vanilla` | default CPU implementation                                                                                             |
| `AVX2   ` | AVX vectorization optimized program using [*vectorclass*](https://github.com/vectorclass)                              |
| `openmp ` | OpenMP parallelized program                                                                                            |
| `cuda   ` | Parallel [CUDA](https://developer.nvidia.com/cuda-toolkit) program                                                     |
| `hip    ` | Parallel [HIP](https://docs.amd.com/bundle/HIP-Programming-Guide-v5.3/page/Introduction_to_HIP_Programming_Guide.html) |

For cuda compilation one needs to define their CUDA version and architercure either as environment variables or during the make process:

```bash
export CUDA_VERSION=11.6
export CUDA_ARCH=61
make ARCH=cuda
or
make ARCH=cuda CUDA_VERSION=11.6 CUDA_ARCH=61
```
*NOTE: Default cuda version is 11.6 and compute architecture is sm_61*

**HPC platforms**:

| ARCH       | Description                                               |
|------------|-----------------------------------------------------------|
| `lumi      ` | CPU-MPI implementation for LUMI supercomputer             |
| `lumi-hip  ` | GPU-MPI implementation for LUMI supercomputer using HIP   |
| `mahti     ` | CPU-MPI implementation for MAHTI supercomputer            |
| `mahti-cuda` | GPU-MPI implementation for MAHTI supercomputer using CUDA |

## User guide

Now that HILA has been built successfully, the next step is to build your first HILA application: [hila application guide](./docs/guide/hila_applications.md)

After building your first HILA application one can move on to the comprehensive guide, which describes everything that HILA has to offer: [comprehensive guide](./docs/guide/hila_functionality.md)

Both of these resources can be viewed on the web guide hosted on: TODO: add link to github pages or hosted guide somewhere

To generate the user guide and technical documentation locally one can run:

    doxygen /docs/config 

To open the documentation locally with any browser:

    firefox /docs/html/index.html