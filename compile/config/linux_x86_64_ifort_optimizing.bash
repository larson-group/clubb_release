# $Id$
# Makefile definitions for Linux x86_64 using the Intel Fortran 
# compiler. The compiler options under OPTIMIZE and ARCH are 
# configured to match the compiler options used by NCAR to 
# compile CESM. These definitions are intended for use when
# performance testing.


# Fortran 95 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# == Machine specific options ==
ARCH="-xCORE_AVX2" # Uses Intel CPU's vector units with 256bit length

# == Optimization ==
OPTIMIZE="-assume realloc_lhs -no-fma -O2 -fp-model source"

# == Debug ==
DEBUG="-g -debug all -traceback" # Provides profiling applications such as VTune with extra info

FFLAGS="$ARCH $OPTIMIZE $DEBUG"

# == Warnings ==
WARNINGS="-warn -warn notruncated_source"

# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-real-size 64"

# == NetCDF Location ==
#Variable defined in larson-group.sh, see here (https://github.com/larson-group/sys_admin/blob/master/set_larson-group_paths/larson-group.sh)
NETCDF="$IFORT_NETCDF_FORTRAN"

# == LAPACK libraries ==
# AMD Core Math Library
#ACML="/opt/acml5.1.0/ifort64/lib"
#LAPACK="-L$ACML -Wl,-rpath,$ACML -lacml"
# Intel MKL
LAPACK="-mkl=sequential"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$NETCDF/lib -lnetcdff $LAPACK"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -DMKL enables MKL solver (PARDISO/GMRES) support
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=8 -DMKL"
CPPFLAGS="-I$MKLPATH/../../include -I$NETCDF/include"

# == Static library processing ==
AR=ar
#AR=xiar # Interprocedural optimization archive tool
ARFLAGS=cru
RANLIB=ranlib

# == Shared library processing ==
SHARED=$FC
SHAREDFLAGS="-fPIC -shared"

# Location of 'mkmf' utility
mkmf=$dir/mkmf

# gmake command to use and options: '-j 2' enables parallel compilation
gmake="make -j $(nproc)"
