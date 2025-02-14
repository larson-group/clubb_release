# $Id$

# Makefile definitions customized for Linux x86_64 using the Cray compiler
# see https://github.com/larson-group/clubb/issues/1138 
# and https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#id4
# and https://cpe.ext.hpe.com/docs/cce/man7/intro_openacc.7.html#compiling
# for more info. This is intended to be used with the following modules 
# module load craype-accel-amd-gfx90a cce/17.0 rocm cray-mpich cray-hdf5-parallel cray-netcdf-hdf5parallel cray-python
# Currently Loaded Modules:
#   1) craype-x86-trento        5) xpmem/2.8.4-1.0_7.3__ga37cbd9.shasta   9) PrgEnv-cray/8.5.0  13) lfs-wrapper/0.0.1        17) cce/15.0.1                   21) cray-netcdf-hdf5parallel/4.9.0.7
#   2) libfabric/1.20.1         6) cray-pmi/6.1.13                       10) Core/24.07         14) DefApps                  18) rocm/5.4.0                   22) cray-python/3.11.5
#   3) craype-network-ofi       7) craype/2.7.31.11                      11) tmux/3.4           15) craype-accel-amd-gfx90a  19) cray-mpich/8.1.27
#   4) perftools-base/23.12.0   8) cray-dsmml/0.2.2                      12) hsi/default        16) cray-libsci/23.09.1.1    20) cray-hdf5-parallel/1.12.2.7

# Inactive Modules:
#   1) darshan-runtime


# Fortran 95 compiler and linker
FC=ftn
LD=ftn

# Define path to directories
dir=`pwd` # dir where this script resides
bindir="$dir/../bin"  # dir for Makefile and executable
objdir="$dir/../obj"  # dir for *.o and *.mod files
libdir="$dir/../lib"  # dir for *.a library files
srcdir="$dir/../src"  # dir where the source files reside

# == Debug ==
#   "-eD" gives debug info with source, great for gdb, but makes collapse() not work
DEBUG=""

# == Optimization ==
OPTIMIZE="-O2"

# == NetCDF Location ==
#NETCDF="/opt/cray/pe/netcdf/default/CRAYCLANG/14.0/"
NETCDF=$(nf-config --prefix)

# == Linking Flags ==
LDFLAGS="-L$NETCDF/lib -lnetcdff -h noomp -h acc"

FFLAGS="$OPTIMIZE $DEBUG -hnopattern -h noomp -h acc -h list=a -h msgs -rm"  

# Preprocessing Directives:
#   -DNETCDF enables netCDF output
CPPDEFS="-DNETCDF -DCLUBB_REAL_TYPE=8 -DCLUBB_GPU"
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

# gmake command that uses "2+num_cores/2" cores for compilation
gmake="make -j$(echo "scale=0 ; 2+$(nproc) / 2.0" | bc)"
#gmake="make -j1"

