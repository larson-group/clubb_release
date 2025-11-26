# Define general dependency information in case we need to build our deps
include(ExternalProject)

find_package(PkgConfig REQUIRED)
pkg_check_modules( NETCDF_FORTRAN netcdf-fortran )

if (NETCDF_FORTRAN_FOUND)
  message(STATUS "netcdf-fortran found via pkg-config")
  message(STATUS "Include dirs: ${NETCDF_FORTRAN_INCLUDE_DIRS}")
  message(STATUS "Libraries:    ${NETCDF_FORTRAN_LIBRARIES}")

  # Example: derive a root path from first include dir
  list(GET NETCDF_FORTRAN_INCLUDE_DIRS 0 first_inc)
  cmake_path(GET first_inc PARENT_PATH NetCDFFortran_ROOT)
  message(STATUS "Derived NetCDFFortran_ROOT = ${NetCDFFortran_ROOT}")
endif()

# Setup RPaths so linked executables will be able to find our external deps
#set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH};${NetCDFFortran_ROOT})
set(BUILD_RPATH ${CMAKE_INSTALL_RPATH})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Netcdf Fortran lib
add_library(netcdf-fortran SHARED IMPORTED)
set_target_properties(netcdf-fortran PROPERTIES
  IMPORTED_LOCATION "${NetCDFFortran_ROOT}/lib/libnetcdff${CMAKE_SHARED_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_FORTRAN_INCLUDE_DIRS}"
)


# Note: This file was adapted from a much more complex version, which was aparently 
#       capable of compiling netcdf from source if it could not be found in the 
#       environment. This didn't seem to work initially, but would certainly be
#       valuable for outsiders compiling CLUBB who don't have netcdf-fortran
#       installed already. 
#
#       This initial version is given below. The way it sets up netcdf-fortran
#       is also incompatible with the new cmake code, so this is meant mainly
#       to be used as reference in case someone wants to take on the task
#       of getting netcdf-fortran to compile here.

# # External dependency locations
# set(NetCDF_C_REPO https://github.com/unidata/netcdf-c.git CACHE STRING "Git repository to fetch NetCDF-C from")
# set(NetCDF_Fortran_REPO https://github.com/unidata/netcdf-fortran.git CACHE STRING "Git repository to fetch NetCDF-Fortran from")
# set(HDF5_REPO https://github.com/HDFGroup/hdf5.git CACHE STRING "Git repository to fetch HDF5 from")
# set(Zlib_REPO https://github.com/madler/zlib.git CACHE STRING "Git repository to fetch zlib from")
# set(LAPACK_URL https://www.netlib.org/lapack/lapack-3.5.0.tgz CACHE STRING "URL to fetch lapack archive from")
# set(CURL_REPO https://github.com/curl/curl.git CACHE STRING "Git repository to fetch Curl from")

# # Specific dep hash/versions/tags to checkout
# set(NetCDF_C_HASH "v4.9.0" CACHE STRING "Version of NetCDF-C to build against (can be a git branch, tag, or hash)")
# set(NetCDF_Fortran_HASH "v4.5.4" CACHE STRING "Version of NetCDF-Fortan to build against (can be a git branch, tag, or hash)")
# set(HDF5_HASH "hdf5-1_12_2" CACHE STRING "Version of HDF5 to build against (can be a git branch, tag, or hash)")
# set(Zlib_HASH "v1.2.13" CACHE STRING "Version of Zlib to build against (can be a git branch, tag, or hash)")
# set(LAPACK_HASH 9ad8f0d3f3fb5521db49f2dd716463b8fb2b6bc9dc386a9956b8c6144f726352 CACHE STRING "Lapack archive sha256 hash")
# set(CURL_HASH "curl-8_2_1" CACHE STRING "Version of Curl to be built (can be a git branch, tag, or hash)")

# include(ExternalProject)
# # Before trying to build anything ourselves, we'll try to detect what's already on the system
# find_package(lapack)
# if (LAPACK_FOUND)
#     # We've found an external lapack, define associated link interface
#     add_link_options(${LAPACK_LINKER_FLAGS})
#     link_libraries(${LAPACK_LIBRARIES})
# else()
#     set(LAPACK_ROOT ${CMAKE_BINARY_DIR}/install_lapack)
#     ExternalProject_Add(
#         clubb_lapack_ext
#         URL ${LAPACK_URL}
#         URL_HASH SHA256=9ad8f0d3f3fb5521db49f2dd716463b8fb2b6bc9dc386a9956b8c6144f726352
#         CMAKE_ARGS
#             -DCMAKE_INSTALL_PREFIX=${LAPACK_ROOT}
#     )
#     add_library(clubb_local_lapack STATIC IMPORTED)
#     set_target_properties(clubb_local_lapack PROPERTIES IMPORTED_LOCATION ${LAPACK_ROOT}/lib/liblapack.a)
#     add_library(clubb_blas STATIC IMPORTED)
#     set_target_properties(clubb_blas PROPERTIES IMPORTED_LOCATION ${LAPACK_ROOT}/lib/libblas.a)
#     link_libraries(clubb_local_lapack clubb_blas)
# endif()

# find_package(ZLIB)
# if(${ZLIB_FOUND})
#     cmake_path(GET ZLIB_LIBRARIES PARENT_PATH result)
#     cmake_path(GET result PARENT_PATH result)
#     # Setting ZLib root will allow any deps we do need to build to
#     # find THIS Zlib so everything is built with a consistent zlib
#     set(ZLIB_ROOT ${result})
# else()
#     # We didn't find Zlib on the system, build it ourselves
#     set(ZLIB_ROOT ${CMAKE_BINARY_DIR}/install_zlib)
#     ExternalProject_Add(
#         clubb_zlib
#         GIT_REPOSITORY    ${Zlib_REPO}
#         GIT_TAG           ${Zlib_HASH}
#         CMAKE_ARGS
#             -DCMAKE_INSTALL_PREFIX=${ZLIB_ROOT}
#     )
#     set(ext_zlib clubb_zlib)
#     set(ZLIB_INCLUDE_DIRS ${ZLIB_ROOT}/include)
# endif()

# find_package(CURL)
# if (CURL_FOUND)
#     cmake_path(GET CURL_INCLUDE_DIRS PARENT_PATH CURL_ROOT)
# else()
#     set(CURL_ROOT ${CMAKE_BINARY_DIR}/install_curl)
#     ExternalProject_Add(
#         clubb_curl
#         GIT_REPOSITORY ${CURL_REPO}
#         GIT_TAG ${CURL_HASH}
#         CMAKE_ARGS
#             -DCMAKE_INSTALL_PREFIX=${CURL_ROOT}
#             -DZLIB_ROOT=${ZLIB_ROOT}
#         DEPENDS ${ext_zlib}
#     )
#     set(ext_curl clubb_curl)
# endif()

# find_package(HDF5)
# if(${HDF5_FOUND})
#     cmake_path(GET HDF5_INCLUDE_DIRS PARENT_PATH HDF5_ROOT)
# else()
#     set(HDF5_ROOT ${CMAKE_BINARY_DIR}/install_hdf5)
#     ExternalProject_Add(
#         clubb_hdf5
#         GIT_REPOSITORY ${HDF5_REPO}
#         GIT_TAG ${HDF5_HASH}
#         CMAKE_ARGS
#             -DCMAKE_INSTALL_PREFIX=${HDF5_ROOT}
#             -DHDF5_ENABLE_Z_LIB_SUPPORT=ON
#             -DHDF5_BUILD_FORTRAN=ON
#             -DZLIB_ROOT=${ZLIB_ROOT} # This allows us to point the newly built hdf5 at our previously built/detected zlib
#         DEPENDS ${ext_zlib}
#     )
#     set(ext_hdf5 clubb_hdf5)
# endif()

# find_package(netCDF)
# if(${NetCDF_FOUND})
#     cmake_path(GET NetCDF_INCLUDE_DIRS PARENT_PATH NetCDFC_ROOT)
# else()
#     set(NetCDFC_ROOT ${CMAKE_BINARY_DIR}/install_NetCDFC)
#     ExternalProject_Add(
#         clubb_netcdfc
#         GIT_REPOSITORY ${NetCDF_C_REPO}
#         GIT_TAG ${NetCDF_C_HASH}
#         CMAKE_ARGS
#             -DCMAKE_INSTALL_PREFIX=${NetCDFC_ROOT}
#             -DENABLE_BYTERANGE=OFF
#             -DENABLE_NETCDF_4=OFF
#             -DENABLE_HDF5=ON
#             -DHDF5_ROOT=${HDF5_ROOT}
#             -DZLIB_ROOT=${ZLIB_ROOT}
#             # -DCURL_ROOT=${CURL_ROOT}
#             -DBUILD_UTILITIES=ON
#             -DENABLE_LARGE_FILE_SUPPORT=ON
#             -DENABLE_DAP=ON
#             -DCMAKE_INSTALL_LIBDIR=${NetCDFC_ROOT}/lib
#         DEPENDS ${ext_hdf5}
#     )
#     set(ext_netcdfc clubb_netcdfc)
# endif()

# find_package(netcdf-fortran)
# if(${NETCDFFortran_FOUND})
#     cmake_path(GET NetCDFFortran_INCLUDE_DIRS PARENT_PATH NetCDFFortan_ROOT)
# else()
#     set(NetCDFFortran_ROOT ${CMAKE_BINARY_DIR}/install_NetCDFFortran)
#     # The CMake args defined below may seem pedantic but are very neccesary due to some
#     # significant inflexibility in the netcdf-fortran CMake code... essentially the project can't
#     # handle being installed anywhere that is not /usr/* without some serious work, i.e. below
#     ExternalProject_Add(
#         clubb_netcdffortran
#         GIT_REPOSITORY ${NetCDF_Fortran_REPO}
#         GIT_TAG ${NetCDF_Fortran_HASH}
#         CMAKE_ARGS
#             -DnetCDF_ROOT=${NetCDFC_ROOT}
#             -DCMAKE_POLICY_DEFAULT_CMP0074:STRING=NEW
#             -DCMAKE_INSTALL_LIBDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_BINDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_INFODIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_DOCDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_INCLUDEDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_DOCDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_DATADIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_DATAROOTDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_LOCALEDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_LOCALSTATEDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_MANDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_OLDINCLUDEDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_RUNSTATEDIR=${NetCDFFortran_ROOT}
#             -DCMAKE_INSTALL_PREFIX=${NetCDFFortran_ROOT}
#             -DENABLE_DAP=ON
#         DEPENDS ${ext_netcdfc}
#     )
# endif()

# # Setup RPaths so linked executables
# # will be able to find our external deps
# set(CMAKE_MACOSX_RPATH TRUE)
# set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH};${NetCDFFortran_ROOT};${NetCDFC_ROOT}/lib)
# set(BUILD_RPATH ${CMAKE_INSTALL_RPATH})
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# # Expose NetCDF-Fortran as a target so we can consume it
# add_library(netcdf-fortran SHARED IMPORTED)
# add_library(netcdf-c SHARED IMPORTED)
# set_target_properties(netcdf-fortran PROPERTIES IMPORTED_LOCATION ${NetCDFFortran_ROOT}/libnetcdff${CMAKE_SHARED_LIBRARY_SUFFIX})
# set_target_properties(netcdf-c PROPERTIES IMPORTED_LOCATION ${NetCDFC_ROOT}/lib/libnetcdf${CMAKE_SHARED_LIBRARY_SUFFIX})
# set(NetCDFFortran_INCLUDE_DIRS ${NetCDFFortran_ROOT})

# add_library(NetCDF::NetCDFC ALIAS netcdf-c)
# add_library(NetCDF::NetCDFFortran ALIAS netcdf-fortran)
# # Ensure all External projects are built before anything from CLUBB
# add_dependencies(netcdf-fortran clubb_netcdffortran)