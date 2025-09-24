# CMake module providing CPPDefinitions

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/MOD)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# enable "ifdef radoffline/nooverlap/COAMPS_MICRO" 
add_compile_definitions(radoffline nooverlap CLUBB COAMPS_MICRO)


if(${USE_NetCDF})
    # enable "ifdef NETCDF"
    add_compile_definitions(NETCDF)
endif()

if(${TUNING})
    # enable "ifdef TUNING"
    add_compile_definitions(TUNER)
endif()

if(${SILHS})
    # enable "ifdef SILHS"
    add_compile_definitions(SILHS)
endif()


# Flags used to define GPU usage. These may neccesary for CPU builds too, as sometimes
# flags are required to explicitly disable GPU usage (e.g. -noacc for nvhpc)
add_compile_options(${GPU_COMPILE_FLAGS})
add_link_options(${GPU_LINK_FLAGS})
add_compile_definitions(${GPU_DEFINITIONS})


# Preprocessor definition for CLUBB_REAL_TYPE, which sets the floating-point precision
set(PRECISION "double" CACHE STRING "Precision for CLUBB_REAL_TYPE (single,double,quad)")
set_property(CACHE PRECISION PROPERTY STRINGS single double quad)

if (PRECISION STREQUAL "single")
  message(STATUS "Compiling with single precision")
  add_compile_definitions(CLUBB_REAL_TYPE=4)

elseif (PRECISION STREQUAL "double")
  message(STATUS "Compiling with double precision")
  add_compile_definitions(CLUBB_REAL_TYPE=8)

elseif (PRECISION STREQUAL "quad")
  message(STATUS "Compiling with quadruple precision")
  add_compile_definitions(CLUBB_REAL_TYPE=16)

else()
  message(FATAL_ERROR "PRECISION must be one of: single, double, quad")
endif()