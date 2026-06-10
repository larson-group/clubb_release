# Building

## CMake Build System Quickest Start

### Dependencies

| Dependency | Version |
|------------|---------|
| CMake      | >= 3.18 |
| Unix Style Make | System version is fine |

optionally, [Ninja](https://ninja-build.org/) is another build system CMake is capable of generating 
which is capable of highly parallelized builds often beating make's j flag in terms of build times.
For larger builds, Ninja >= 1.10 is recommended but not required.

### Quick and Easy CLUBB Compilation (Quickest Start)

Other sections detail more information about the underlying CMake system and elaborate on 
some of the details/more complex usage, but if you're just interested in getting CLUBB 
compiled and running as quickly and easily as possible, read on as this section is for you.

To compile CLUBB, simply navigate to CLUBB's source tree root, set up your 
environment, and invoke the compile script. For example, compiling with intel could be done via
```
module load intel netcdf-fortran
./compile.py
```

Or, e.g., to compile with nvfortran and for the GPU using OpenACC
```
module load nvhpc netcdf-fortran
./compile.py -gpu openacc
```

However, if you don't have `module` set up, you can just do the following
```
FC=<compiler> ./compile.py
```
The following 4 compilers are supported, each of which has 2 or more aliases: 
 - gfortran, gcc, gnu
 - ifx, ifort, intel, intel-oneapi, intel-classic
 - nvfortran, nvhpc
 - crayftn, cce, ftn 

If no compiler is specified, `compile.py` will default to gfortran! 

`compile.py` will create a subfolder named `build/<compiler>_*` and place all compiled libraries and 
binaries in subdirectories corresponding to the corresponding source code directories. o

A quick set of tests can be run after compilation by including the -run_tests flag 
```
./compile.py -run_tests
```

Other various settings are possible as well, run with -h/--help to see a list 
```
./compile.py -h
```

If you want a more nuanced CLUBB build or want more control over various facets of compilation or 
are interested in cross compiling, see the other sections which dive further into detail regarding 
the CMake system.

### Specifying the compiler

The compile script (compile.py) attempts to detect a compiler from the environment variable "FC", 
and tries to set the toolchain appropriately. If neither `FC` nor `LMOD_FAMILY_COMPILER` is set, it 
will also try to find `gfortran` on `PATH` and use the corresponding default toolchain.

If no compiler can be found, a specific toolchain will need to be specified. These can be found 
in <CLUBB_ROOT>/cmake/toolchains.

This system works well when the environment is setup using lmod (module system), but more work may 
be needed to make compilation function correctly if the environment is not set properly.

### Important CMake options

These options define compiler settings that are rarely changed. To modify them, the easiest way is 
to use compiler flags: 

| Compiler Flags          | Description                                   | Default Value         |
|-------------------------|-----------------------------------------------|-----------------------|
|  -disable_netcdf        |  Disable NetCDF output support                |  Disabled (NetCDF on) |
|  -openmp                |  Enable OpenMP threading                      |  Disabled             |
|  -tuning                |  Enable TUNING mode with extra runtime checks |  Disabled             |
|  -disable_silhs         |  Disable SILHS                                |  Disabled (SILHS on)  |
|  -run_tests             |  Run ctests after compilation                 |  Disabled             |


However, you can also edit the top level CMakeLists.txt (clubb/CMakeLists.txt):

| CMake Option            | Description                     |    Default Value  |
|-------------------------|---------------------------------|-------------------|
|   USE_NetCDF            | Use netcdf for model output     |       ON          |
|   ENABLE_OMP            | Enable CLUBB to run with OpenMP |       OFF         |
|   TUNING                | Compile CLUBB for tuning        |       OFF         |
|   SILHS                 | Compile for SILHS capability    |       ON          |
|   ENABLE_TESTS          | Enable running tests            |       ON          |


## CMake Build System

CMake is a build system generator that creates build systems dynamically given the system of 
execution and some configuration parameters. This functionality is leveraged here to facilitate 
building the CLUBB executables and libraries and is intended to replace the previous 
system in <CLUBB_ROOT>/compile/.

For a primer on CMake, see the CMake course: <https://cmake.org/cmake/help/book/mastering-cmake/> 
this will provide boilerplate coverage on basic concepts that will be required to maintain, 
modify, and expand the current CMake system
