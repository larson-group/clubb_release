# $Id$
# Makefile definitions customized for Linux x86_64 using the Intel Fortran 
# compiler 


# Fortran 95 compiler and linker
FC=ifx
LD=ifx

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# == Debug ==
DEBUG="-debug full -traceback -check bounds -fpe0 -ftz -ftrapuv -init=snan,arrays"

# == Warnings ==
WARNINGS="-warn all -warn notruncated_source"

# These flags are used to disable warnings for code that we don't want to see
# warnings for (like code that isn't ours).
DISABLE_WARNINGS="-warn none"

# == Machine specific options ==
#ARCH="-xHost -qopenmp" # This should work on most modern AMD/Intel computers
# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-real-size 64"

# == Optimization ==
# No optimization
OPTIMIZE="-O0"
#OPTIMIZE="-O3"

# == NetCDF Location ==
#Variable defined in larson-group.sh, see here (https://github.com/larson-group/sys_admin/blob/master/set_larson-group_paths/larson-group.sh)
NETCDF="$(nf-config --prefix)"

# == LAPACK libraries ==
# AMD Core Math Library
#ACML="/opt/acml5.1.0/ifort64/lib"
#LAPACK="-L$ACML -Wl,-rpath,$ACML -lacml"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdff $LAPACK"

FFLAGS="$ARCH $OPTIMIZE $DEBUG"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -DMKL enables MKL solver (PARDISO/GMRES) support
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
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
SHAREDFLAGS="-fPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j5"
