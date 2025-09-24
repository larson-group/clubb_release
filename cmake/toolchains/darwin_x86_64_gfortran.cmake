set(CMAKE_SYSTEM_NAME Darwin)
set(CMAKE_SYSTEM_PROCESSOR x86_64)
set(CMAKE_SYSTEM_VERSION "${CMAKE_HOST_SYSTEM_VERSION}")
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "CLUBB Fortran Compiler")

add_compile_definitions(__GFORTRAN__)
add_compile_definitions(CLUBB_REAL_TYPE=8)


set(CMAKE_Fortran_FLAGS "-O2")

set(CMAKE_Fortran_FLAGS_DEBUG  "-g -fbounds-check -dynamiclib")

set(WARNINGS "-Wall -pedantic"
)

set(DOUBLE_PRECISION_FLAG "-fdefault-real-8")

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${WARNINGS}")