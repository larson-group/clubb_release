# $Id$
# Makefile definitions for Oracle Solaris Studio Fortran compiler (SPARC or x64)
# Tested with version 12.3

# Fortran 95 compiler and linker
FC=sunf95 # Sun/Oracle Fortran
LD=sunf95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using clubb_tuner.
# In Sun f95: remove the reference to ftrap, or use -ftrap=common
DEBUG="-g -C -fns=no -ftrap=%none" #-stackvar -xcheck=init_local -ftrap=common

WARNINGS="-w3 -ansi"
# == Machine specific flags ==
# Note that when linking to sunperf (for LAPACK) you must use -dalign
ARCH="-m64 -xarch=native -dalign"

# == Used to promote all real's to double precision ==
# Note the on SPARC the compiler will allow as large as quad precision if desired, 
# but the included LAPACK library only allows up to double precision.
DOUBLE_PRECISION="-xtypemap=real:64,double:64,integer:32"

# == Optimization ==
# These are all pretty conservative options, check the your compiler manual 
# for information on using more aggressive techniques (inlining, etc.)
# -KPIC would need to be here to to build shared libraries.
# Note that -g is needed for both profiling and debugging.
# These options are valid on a SPARC and x64 machine
OPTIMIZE="-g -xO3 -xvector=lib -ftrap=%none"

# == NetCDF Location  ==
# Note: Solaris studio 12.3 does not appear to compile the latest netCDF.
# The 3.6.x version will compile and pass "make test".
NETCDF="/usr/local/netcdf-3.6.3"

# == LAPACK libraries ==
LAPACK="-xlic_lib=sunperf"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$OPTIMIZE $ARCH"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# -Bstatic; for static linking of libraries
# -xlibmopt; link the optimized version of libm (this can affect results)
LDFLAGS="$FFLAGS -L$NETCDF/lib -lnetcdf -xlibmopt $LAPACK"


# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# You will need to `make clean' if you change these
# Use -M<include path> to set a module file directory
CPPDEFS="-DNETCDF -DNO_LAPACK_ISNAN -DCLUBB_REAL_TYPE=8"
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
