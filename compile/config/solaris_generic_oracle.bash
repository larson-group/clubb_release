# $Id$
# Makefile definitions for Sun Studio's Fortran compiler (SPARC or x64)

# Fortran 95 compiler and linker
FC=f95 # Sun Fortran
LD=f95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
# In Sun f95: remove the reference to ftrap, or use -ftrap=common
DEBUG="-g -C -fns=no -ftrap=%none" #-stackvar -xcheck=init_local -ftrap=common

WARNINGS="-w3 -ansi"
# == Machine specific flags ==
# The 2nd line is specific to ketley.math.uwm.edu.
# Note that when linking to sunperf (for LAPACK) you must use -dalign
ARCH="-xarch=native64 -xcache=native -xchip=native -dalign"
#ARCH="-m64 -xarch=sse3a -xcache=native -xchip=opteron -dalign"

# == Optimization ==
# Using -xtypemap it is possible to compile for quad precision floating 
# point on SPARC, but on AMD64/EMT64 this is still not implemented.
# These are all pretty conservative options, check the your compiler manual 
# for information on using more aggressive techniques (inlining, etc.)
# -KPIC would need to be here to to build shared libraries.
# Note that -g is needed for profiling.
# These options are valid on a SPARC and AMD64 machine
OPTIMIZE="-g -xO3 -xvector=lib -ftrap=%none"
# This works on x86/x64 only
#OPTIMIZE="-g -xO3 -xvector=simd -ftrap=%none"
# These are more aggresive options that enable OpenMP
#OPTIMIZE="-g -xopenmp -xparallel -stackvar -xipo=2 -xvector=simd -xO4 -ftrap=%none"

# == NetCDF Location  ==
NETCDF="/usr/local/netcdf-sun64"

# == LAPACK libraries ==
LAPACK="-xlic_lib=sunperf"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# -Bstatic; for static linking of libraries
# -xlibmopt; link the optimized version of libm (this can affect results)
LDFLAGS="-L$NETCDF/lib -lnetcdf -xlibmopt $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$OPTIMIZE $ARCH"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# You will need to `make clean' if you change these
# Use -M<include path> to set a module file directory
CPPDEFS="-DNETCDF -Dnooverlap -Dradoffline"
CPPFLAGS="-M$NETCDF/include"

# == Static library processing ==
AR="ar"
ARFLAGS="cru"
RANLIB="ranlib"

# == Shared library processing ==
SHARED="$FC"
SHAREDFLAGS="-G"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="gmake"
