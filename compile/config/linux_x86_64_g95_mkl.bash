# $Id: linux_ia32_g95_optimize.bash 4627 2010-06-30 17:23:50Z meyern@uwm.edu $
# Configuration file for a Linux machine using g95
# An MS Windows machine with MSYS or Cygwin should be similar
# This configuration file links to the Intel MKL.

# Fortran 95 compiler and linker
FC=g95
LD=g95

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# It is sometimes helpful to turn on floating-point trapping for the 
#  standalone program, but this will not work when using the tuner.
# In g95: `setenv G95_FPU_INVALID TRUE' or `export G95_FPU_INVALID=TRUE'
# These are the options for debugging symbols, bounds checking, IEEE-754,
# and a special g95 flag for debugging, respectively. 
DEBUG="-g -fbounds-check -mieee-fp -ftrace=full"

# == Warnings ==
# This is the set of preferred warnings to use when compiling CLUBB with g95.
# Warnings ignored are:
# 	142 - Nonblock DO statement is obsolescent 
# 	165 - Implicit interface is called
# 	167 - PRIVATE module procedure is never invoked
WARNINGS="-Wall -Wextra -Wno=142,165,167 -pedantic"

# == Machine specific flags ==
# Note: some of these are 64 bit architectures, so make sure NetCDF is
# compiled accordingly.
#ARCH="-march=pentium4 -msse2 -mfpmath=sse" # Old P4s
#ARCH="-march=nocona -msse3 -mfpmath=sse" # New P4s
#ARCH="-march=nocona -msse3 -mfpmath=sse -r8"# New P4s, double precision
#ARCH="-march=k8 -msse3 -mfpmath=sse -r8" # New Opterons, double precision
ARCH="-march=k8 -msse3 -mfpmath=sse" # New Opterons

# == Optimization ==
# These are all pretty conservative options. Check your compiler manual
# for information on using more aggressive techniques (inlining, etc.)
OPTIMIZE="-O3"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf-g95"

# == LAPACK libraries ==
#LAPACK="-L/usr/local/lib -llapack -lblas" #  The netlib reference LAPACK/BLAS
#LAPACK="-L/usr/local/lib -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas" # ATLAS BLAS (faster)

# == Intel MKL libraries (PARDISO/GMRES solvers) ==
# Intel Math Kernel Library (v11.0)
# MKLPATH="/opt/intel/Compiler/11.0/081/mkl/lib/em64t"
# Intel Math Kernel Library (v11.1)
MKLPATH="/opt/intel/Compiler/11.1/072/mkl/lib/em64t"
LAPACK="-L$MKLPATH -Wl,-rpath,$MKLPATH -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack -lmkl_solver_lp64 -lmkl_core -lguide -lpthread"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# Use -Wl,-rpath <library path> to set a search path for shared libs
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK"

# == Compiler flags ==
# You will need to `make clean' if you change these
FFLAGS="$ARCH $OPTIMIZE $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -DMKL enables MKL solver (PARDISO and GMRES) support
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPFLAGS="-I$MKLPATH/../../include -I$NETCDF/include"
CPPDEFS="-DNETCDF -DMKL -DCLUBB_REAL_TYPE=8"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-fPIC -shared"



# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j5"
