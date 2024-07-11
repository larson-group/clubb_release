# $Id$
# Makefile definitions customized for Linux IA32 using the Portland Group
# compiler

# Fortran 95 compiler and linker
FC=nvfortran
LD=nvfortran

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debugging ==
# Debugging information and floating-point trapping
DEBUG="-traceback -g -C -Kieee -Ktrap=fp -⁠Minfo"

# == Machine specific options ==
# The NVHPC Fortran compiler will select the native processor type by default
ARCH="-Mcache_align" # -Mcache_align is included for the use of the ACML

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-r8"

# == Optimization ==
OPTIMIZE="-O0"

# == NetCDF Location ==
#Variable defined in larson-group.sh, see here (https://github.com/larson-group/sys_admin/blob/master/set_larson-group_paths/larson-group.sh)
NETCDF="$(nf-config --prefix)"

# == Linking Flags ==
LDFLAGS="$ARCH -L$NETCDF/lib -lnetcdff $LAPACK"

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

# gmake command to use and options: '-j #n' enables parallel compilation
# $(echo "scale=0 ; 2+$(nproc) / 2.0" | bc) will compute 2+nproc/2, where nproc is
# the number of logical cores available. This is always a close to optimal choice.
gmake="make -j$(echo "scale=0 ; 2+$(nproc) / 2.0" | bc)"
