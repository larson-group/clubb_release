# $Id$
# Makefile definitions customized for Linux IA32 using the Portland Group
# compiler

# Fortran 95 compiler and linker
FC=pgfortran
LD=pgfortran

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
# No Debug flags
DEBUG=""
# Debugging information and floating-point trapping
#DEBUG="-g -C -Kieee -Ktrap=fp"

# == Machine specific options ==
# The PGI Fortran compiler will select the native processor type by default
ARCH="-Mcache_align" # -Mcache_align is included for the use of the ACML

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-r8"

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
#NETCDF="$HOME/netcdf-3.6.3"
#Variable defined in larson-group.sh, see here (https://github.com/larson-group/sys_admin/blob/master/set_larson-group_paths/larson-group.sh)
NETCDF="$PGI_NETCDF_FORTRAN"

# == LAPACK libraries ==
# The PGI directory contains static versions of LAPACK and BLAS
#LAPACK="-llapack -lblas"
# This will select the version of ACML that PGI provides, which is generally
# faster than the reference BLAS and LAPACK (above)
#LAPACK="-lacml"

# == Linking Flags ==
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdff $LAPACK -acc -Mcuda"

FFLAGS="$ARCH $OPTIMIZE $DEBUG -Mbackslash -Mstandard -Kieee -acc"

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
