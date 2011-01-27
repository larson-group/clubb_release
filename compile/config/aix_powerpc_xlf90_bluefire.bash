# $Id$
# Makefile definitions customized for IBM AIX on PowerPC using the XL
# Fortran compiler.


# Fortran 95 compiler and linker
F90=xlf90_r
F77=xlf_r
LD=xlf90_r

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

DEBUG="-g"

# == Warnings ==
WARNINGS=""

# == Machine specific options ==
#ARCH="-DNONSTANDARD_SYSTEM_SUBR -DNATIVE_MASSV -DCLUBB" 
ARCH="-qnoescape" 
#ARCH="-mssse3 -fp-model precise"# Core2 Duo (overlie)

# == Optimization ==
OPTIMIZE="-O3"

# == NetCDF Location ==
NETCDF="/usr/local/netcdf"

# == LAPACK libraries ==
# Intel Math Kernel Library (v11.1)
#MKLPATH="/opt/intel/Compiler/11.1/064/mkl/lib/em64t"
# Intel Math Kernel Library (v10)
#MKLPATH=/opt/intel/mkl/10.0.5.025/lib/em64t 
# Intel Math Kernel Library (v10.2)
#MKLPATH="/opt/intel/mkl/10.2.5.035/lib/em64t"
#LAPACK="-L$MKLPATH -Wl,-rpath,$MKLPATH -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_solver_lp64 -lmkl_core -lguide -lpthread"
# (v8)
#LAPACK="-L/opt/intel/mkl/8.1/lib/64 -Wl,-rpath,/opt/intel/mkl/8.1/lib/64 -lmkl_lapack64 -lmkl_i2p -lmkl -lmkl_vml_i2p -lmkl_vml -lvml -lguide -lpthread"
# Generic library
#LAPACK="-llapack -lblas -lgfortran"
#LAPACK="-L/usr/local/lib -llapack -L/usr/local/atlas/lib -lf77blas -lcblas -latlas" # ATLAS BLAS (faster then generic library)
# Use Bluefire's native IBM XL Fortran version of LAPACK: 
#LAPACK="-qessl"
LAPACK="-L/contrib/lapack/3.2.2/lib"
BLAS="-L/contrib/blas/lib"


# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
# "-lessl" links the IBM XL Fortran version of LAPACK
#LDFLAGS="-L$NETCDF/lib -lnetcdf -lessl"
LDFLAGS="-L$NETCDF/lib -lnetcdf $LAPACK -llapack $BLAS -lblas"

#FFLAGS="$ARCH $DEBUG $LAPACK"
FFLAGS="$ARCH $OPTIMIZE"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -Dradoffline and -Dnooverlap (see bugsrad documentation)
#   -DMKL enables MKL solver (PARDISO/GMRES) support
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

#CPPFLAGS="-DNETCDF -I$MKLPATH/../../include -I$NETCDF/include -Dnooverlap -Dradoffline -DMKL"
#CPPFLAGS="-WF,-DNETCDF,-qmoddir=$NETCDF/include,-Dnooverlap,-Dradoffline"
#CPPFLAGS="-I$NETCDF/include -WF,-DNETCDF,-Dnooverlap,-Dradoffline"
CPPFLAGS="-I$NETCDF/include"
CPPDEFS="-DNETCDF -Dnooverlap -Dradoffline"

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
gmake="make"
