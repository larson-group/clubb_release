# Try to find OpenMP Fortran support using CMake's built-in module
find_package(OpenMP COMPONENTS Fortran REQUIRED)

if(OpenMP_Fortran_FOUND)
  message(STATUS "OpenMP Fortran support found")
  message(STATUS "Flags:    ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Includes: ${OpenMP_Fortran_INCLUDE_DIRS}")
  message(STATUS "Libraries:${OpenMP_Fortran_LIBRARIES}")
else()
  message(FATAL_ERROR "OpenMP Fortran support not found. Disable ENABLE_OMP or install a compiler with OpenMP support.")
endif()

# Define an imported library for consistency with other dependencies
add_library(openmp-fortran INTERFACE IMPORTED)
set_target_properties(openmp-fortran PROPERTIES
  INTERFACE_COMPILE_OPTIONS "${OpenMP_Fortran_FLAGS}"
  INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_Fortran_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${OpenMP_Fortran_LIBRARIES}"
)