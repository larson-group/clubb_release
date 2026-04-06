# CMake helpers for suppressing warnings on vendor or third-party Fortran code.
#
# Warning-suppression flags are compiler-specific, so we centralize the mapping
# here and keep usage at call sites compiler-agnostic.

function(get_vendor_specific_warning_supression_flags out_var)
  set(flags "")

  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(flags "-w")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(flags "-w")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    set(flags "-w")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(flags "-w")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(flags "-w")
  endif()

  set(${out_var} "${flags}" PARENT_SCOPE)
endfunction()

function(disable_compiler_warnings target)
  get_vendor_specific_warning_supression_flags(vendor_warning_flags)

  if(vendor_warning_flags)
    target_compile_options(${target} PRIVATE ${vendor_warning_flags})
  endif()
endfunction()

function(disable_compiler_warnings_for_sources target)
  get_vendor_specific_warning_supression_flags(vendor_warning_flags)

  if(vendor_warning_flags)
    set_source_files_properties(${ARGN}
      TARGET_DIRECTORY ${target}
      PROPERTIES COMPILE_OPTIONS "${vendor_warning_flags}")
  endif()
endfunction()
