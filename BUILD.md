# Building

## CMake Build System Quickest Start

### Dependencies

| Dependency | Version |
|------------|---------|
| CMake      | >= 3.18 |
| Unix Style Make | System version is fine |

optionally, [Ninja](https://ninja-build.org/) is another build system CMake is capable of generating which is capable of highly parallelized builds
often beating make's j flag in terms of build times.
For larger builds, Ninja >= 1.10 is recommended but not required.

### Quick and Easy CLUBB Compilation (Quickest Start)

Other sections detail more information about the underlying CMake
system and elborate on some of the details/more complex useage, but if you're just interested in getting CLUBB compiled and running as quickly and easily as possible, read on as this section is for you.

To just drive a compilation of CLUBB, simply navigate to the source tree root, setup up environment 
and invoke the compile script. For example, compiling with intel could be done via
```
module load intel netcdf-fortran
./compile.py
```

Or to compile with nvfortran and for the GPU using OpenACC
```
module load nvhpc netcdf-fortran
./compile.py -gpu openacc
```

A quick set of tests can be run after compilation by including the -run_tests flage 
```
./compile.py -run_tests
```

This script will detect the compiler, create a subfolder named `build/COMPILER` and drive the build from there, with all compiled libraries/binaries present under directories corresponding to the location of the source code used to build that library/binary.

If you would like all build artifacts in a centralized location, simply specify an installaton location, i.e.

```
./compile.py -install <path-to-install-location>
```

where the path to the installation can be realative to the source tree root or absolute pointing to arbitrary points on the filesystem. If this option is specified, all binaries/libraries will be installed to this location. By default, the script will install all build artifacts into a directory `install/COMPILER` in the root of the CLUBB source tree.


Other various settings are possible as well, run with -h/--help to see a list 
```
./compile.py -h
```

If you want a more nuanced CLUBB build or want more control over various facets of compilation or are interested in cross compiling, see the other sections which dive further into detail regarding the CMake system.

### Specifying the compiler

The compile script (compile.py) attempts to detect a compiler from the environment variable "FC", and tries to set the 
toolchain appropriately. 
If no compiler can be found, a specific toolchain will need to be specified. These can be found in <CLUBB_ROOT>/cmake/toolchains.
This system works well when the environment is setup using lmod (module system), but more work may be needed to make compilation
function correctly if the environment is not set properly.

### Important CMake options

These options define compiler settings that are rarely changed. To modify them, edit 
the top level CMakeLists.txt (clubb/CMakeLists.txt) 

| CMake Option            | Description                     |    Default Value  |
|-------------------------|---------------------------------|-------------------|
|   USE_NetCDF            | Use netcdf for model output     |       ON          |
|   ENABLE_OMP            | Enable CLUBB to run with OpenMP |       OFF         |
|   TUNING                | Compile CLUBB for tuning        |       OFF         |
|   SILHS                 | Compile for SILHS capability    |       ON          |
|   ENABLE_TESTS          | Enable running tests            |       ON          |


## CMake Build System

CMake is a build system generator that creates build systems dynamically given the system of execution and some configuration
parameters. This functionality is leveraged here to faciliate building the CLUBB executables and libraies and is intended to
replace the current system in <CLUBB_ROOT>/compile/.

For a primer on general CMake, see the CMake course: <https://cmake.org/cmake/help/book/mastering-cmake/>
this will provide boilerplate coverage on basic concepts that will be required to maintain, modify, and expand the current CMake system