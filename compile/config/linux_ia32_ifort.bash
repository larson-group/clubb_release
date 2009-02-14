# $Id$
# Makefile definitions customized for Linux x86_64 using the Intel Fortran 
# compiler 


# Fortran 95 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

DEBUG="-g -check bounds -check uninit"

# == Warnings ==
WARNINGS="-warn -warn notruncated_source"

# == Machine specific options ==
ARCH="-xW" # This should work on carson/steele (AMD Opteron)
#ARCH="-xT"# Core2 Duo (overlie)

# == Optimization ==
OPTIMIZE="-O3"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf-intel64"

# == LAPACK libraries ==
# Intel Math Kernel Library (v11)
MKLPATH="/opt/intel/Compiler/11.0/081/mkl/lib/em64t"
# Intel Math Kernel Library (v10)
#MKLPATH=/opt/intel/mkl/10.0.5.025/lib/em64t 
LAPACK="-L$MKLPATH -Wl,-rpath,$MKLPATH -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core -lguide -lpthread"
# (v8)
#LAPACK="-L/opt/intel/mkl/8.1/lib/64 -Wl,-rpath,/opt/intel/mkl/8.1/lib/64 -lmkl_lapack64 -lmkl_i2p -lmkl -lmkl_vml_i2p -lmkl_vml -lvml -lguide -lpthread"
# Generic library
#LAPACK="-llapack -lblas -lgfortran"
#LAPACK="-L/usr/lib64 -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas"# ATLAS BLAS (faster then generic library)

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$libdir -lclubb_param -lclubb_bugsrad -lclubb_coamps -L$NETCDF/lib -lnetcdf $LAPACK"

FCFLAGS="$ARCH $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

CPPFLAGS="-DNETCDF -I$NETCDF/include -Dnooverlap -Dradoffline"

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
gmake="make -j 2"

