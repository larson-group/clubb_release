# $Id$
# Makefile definitions customized for Casper using the NVHPC compiler 
# compiler

# Fortran 95 compiler and linker
FC=nvfortran
LD=nvfortran

OPENACC=true

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
# No Debug flags
DEBUG="-g"
# Debugging information and floating-point trapping
#DEBUG="-g -C -Kieee -Ktrap=fp"

# == Machine specific options ==
# The NVHPC Fortran compiler will select the native processor type by default
ARCH=""

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-r8"

if [ $OPENACC == true ]
then
	# == Optimization ==
	OPTIMIZE="-O3 -acc -gpu=cc70,managed -Minfo=accel"
else
	OPTIMIZE="-O3"
fi


# == NetCDF Location ==
# On Casper, the NETCDF enviornment variable is set by our choice of compiler.
# Running 'module load nvhpc' should set NETCDF correctly, so we can simply comment this line out
#NETCDF="$NETCDF"

# == LAPACK libraries ==
# The NVHPC directory contains static versions of LAPACK and BLAS
LAPACK="-llapack -lblas"

if [ $OPENACC == true ]
then
	# == Linking Flags ==
	LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdff $LAPACK -acc -gpu=cc70,managed -lcurand"
else
	LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdff $LAPACK"
fi

FFLAGS="$ARCH $OPTIMIZE $DEBUG -Mbackslash -Mstandard -Kieee"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
# Define include directories. 
# Need location of include and *.mod files for the netcdf library
CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=8"
CPPFLAGS="-I$NETCDF/include"

# == Static library processing ==
AR=ar
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-KPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j5"
