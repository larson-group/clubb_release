# $Id$
# Makefile definitions customized for Linux x86_64 using the Intel Fortran 
# compiler 

#module purge
#module load intel
#module load netcdf
#module load python
#module load intel-mkl 

source /home/software/spack-0.10.1/opt/spack/linux-centos7-x86_64/gcc-4.8.5/lmod-7.4.9-ic63herzfgw5u3na5mdtvp3nwxy6oj2z/lmod/lmod/init/sh;export MODULEPATH=$MODULEPATH:/software/centos7/spack-latest/share/spack/lmod/linux-centos7-x86_64/Core

module purge

module load cmake/3.20.3-vedypwm intel/20.0.4-lednsve intel-mkl/2020.4.304-voqlapk intel-mpi/2019.9.304-i42whlw netcdf-c/4.4.1-blyisdg netcdf-cxx/4.2-gkqc6fq netcdf-fortran/4.4.4-eanrh5t parallel-netcdf/1.11.0-y3nmmej perl/5.30.3-xxmtnqh 

export NETCDF_C_PATH=/gpfs/fs1/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-20.0.4/netcdf-c-4.4.1-blyisdg

export NETCDF_FORTRAN_PATH=/gpfs/fs1/software/centos7/spack-latest/opt/spack/linux-centos7-x86_64/intel-20.0.4/netcdf-fortran-4.4.4-eanrh5t

# Fortran 95 compiler and linker
FC=ifort
LD=ifort

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside


# == Debug ==
DEBUG=""
#DEBUG="-debug full -traceback -check bounds -check uninit -fpe0 -ftz -ftrapuv"

# == Warnings ==
WARNINGS="-warn -warn notruncated_source"

# == Machine specific options ==
ARCH="-xHost" # This should work on most modern AMD/Intel computers
# == Used to promote all real's to double precision ==
DOUBLE_PRECISION="-real-size 64"

# == Optimization ==
OPTIMIZE="-O3"
#OPTIMIZE="-O3 -ipo" # Interprocedural optimization

# == NetCDF Location ==
#NETCDF=$NETCDF_FORTRAN_PATH

# == LAPACK libraries ==
# AMD Core Math Library
#ACML="/opt/acml5.1.0/ifort64/lib"
#LAPACK="-L$ACML -Wl,-rpath,$ACML -lacml"
# Intel MKL
LAPACK="-qmkl=sequential"

# == Linking Flags ==
# Use -s to strip (no debugging); 
# Use -L<library path> -l<lib> to link in an external library
LDFLAGS="-L$NETCDF_FORTRAN_PATH/lib -lnetcdff -L$NETCDF_C_PATH/lib -lnetcdf -L$MKLROOT/lib/intel64 $LAPACK"

FFLAGS="$ARCH $OPTIMIZE $DEBUG -fp-model strict"

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
#   -DMKL enables MKL solver (PARDISO/GMRES) support
# You will need to `make clean' if you change these
# Use -I<include path> to set a module or header file directory
# Define include directories. 
# Need location of include and *.mod files for the netcdf library

CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=8"
CPPFLAGS="-I$MKLPATH/../../include -I$NETCDF_FORTRAN_PATH/include -I$NETCDF_C_PATH"

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
gmake="make -j5"
