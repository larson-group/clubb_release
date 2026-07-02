# Build integration for F2PY extension(s).
#
# This path is CMake-native and requires USE_NetCDF=ON.
# It builds a wrapper backend library with CMake and then asks numpy.f2py
# to generate the Python extension module that links against that backend.

if(NOT USE_NetCDF)
    message(FATAL_ERROR
        "ENABLE_F2PY requires USE_NetCDF=ON. "
        "Reconfigure without -disable_netcdf (or disable ENABLE_F2PY)."
    )
endif()

# Interpreter selection:
# Prefer the top-level Python_EXECUTABLE (passed from compile.py), so f2py uses
# the same Python environment as the rest of the build. Fallback to find_package.
if(DEFINED Python_EXECUTABLE AND EXISTS "${Python_EXECUTABLE}")
    set(F2PY_PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
else()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    set(F2PY_PYTHON_EXECUTABLE "${Python3_EXECUTABLE}")
endif()

execute_process(
    COMMAND "${F2PY_PYTHON_EXECUTABLE}" -c "import numpy, numpy.f2py"
    RESULT_VARIABLE F2PY_NUMPY_CHECK_RC
)
# Fail early at configure-time if f2py isn't usable in this Python env.
if(NOT F2PY_NUMPY_CHECK_RC EQUAL 0)
    message(FATAL_ERROR "NumPy/F2PY is required for ENABLE_F2PY=ON.")
endif()

# Query Python's extension suffix (e.g. ".cpython-313-x86_64-linux-gnu.so")
# so our generated OUTPUT path matches what f2py/meson actually emits.
execute_process(
    COMMAND "${F2PY_PYTHON_EXECUTABLE}" -c "import importlib.machinery as m; print(m.EXTENSION_SUFFIXES[0])"
    OUTPUT_VARIABLE F2PY_EXT_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT F2PY_EXT_SUFFIX)
    set(F2PY_EXT_SUFFIX ".so")
endif()

message(STATUS "ENABLE_F2PY=ON with USE_NetCDF=ON: reusing CMake-built CLUBB libraries")

# Wrapper layout:
# - `derived_type_storage.F90` defines shared Fortran-derived-type storage.
# - Every callable wrapper source uses `*_wrapper.F90` naming and is picked up here.
file(GLOB F2PY_MODULE_WRAPPERS "${CMAKE_CURRENT_SOURCE_DIR}/f2py_fortran_wrappers/*_wrapper.F90")

set(F2PY_WRAPPER_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/f2py_fortran_wrappers/derived_type_storage.F90"
    ${F2PY_MODULE_WRAPPERS}
)

# Artifact paths:
# - F2PY_WORK_DIR is a transient build-tree directory for invoking f2py.
# - F2PY_RUNTIME_DIR stores the build-specific extension and backend library.
# - F2PY_MODULE_IN_BUILD is what f2py produces before installation.
# - F2PY_MODULE_RUNTIME is the build-specific module imported by selected/latest installs.
set(F2PY_WORK_DIR "${CMAKE_CURRENT_BINARY_DIR}/f2py_work")
set(F2PY_RUNTIME_DIR "${CMAKE_CURRENT_BINARY_DIR}/f2py_runtime")
set(F2PY_INSTALL_SUBDIR "python")
set(F2PY_MODULE_IN_BUILD "${F2PY_WORK_DIR}/clubb_f2py${F2PY_EXT_SUFFIX}")
set(F2PY_MODULE_RUNTIME "${F2PY_RUNTIME_DIR}/clubb_f2py${F2PY_EXT_SUFFIX}")
set(F2PY_PYF_FILE "${CMAKE_CURRENT_SOURCE_DIR}/clubb_f2py.pyf")
set(F2PY_PYF_IN_BUILD "${F2PY_WORK_DIR}/clubb_f2py.pyf")

add_library(clubb_f2py_backend SHARED ${F2PY_WRAPPER_SOURCES})
# Link to clubb_driver_lib so this backend can reuse the full CLUBB dependency
# graph already built by CMake (including clubb_core_lib and friends).
target_link_libraries(clubb_f2py_backend PRIVATE clubb_driver_lib)

# Keep backend next to the Python extension for import-time loading.
set_target_properties(clubb_f2py_backend PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${F2PY_RUNTIME_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${F2PY_RUNTIME_DIR}"
)

set(F2PY_LDFLAG_LIST "-Wl,-rpath,\$ORIGIN")
foreach(_f2py_link_dir IN LISTS CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES)
    if(_f2py_link_dir AND EXISTS "${_f2py_link_dir}")
        list(APPEND F2PY_LDFLAG_LIST
            "-L${_f2py_link_dir}"
            "-Wl,-rpath,${_f2py_link_dir}"
        )
    endif()
endforeach()
list(JOIN F2PY_LDFLAG_LIST " " F2PY_LDFLAGS)

add_custom_command(
    OUTPUT "${F2PY_MODULE_IN_BUILD}"
    # Start from a clean transient f2py work directory to avoid stale artifacts.
    COMMAND "${CMAKE_COMMAND}" -E rm -rf "${F2PY_WORK_DIR}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${F2PY_WORK_DIR}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${F2PY_RUNTIME_DIR}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${F2PY_PYF_FILE}" "${F2PY_PYF_IN_BUILD}"
    COMMAND "${CMAKE_COMMAND}" -E chdir "${F2PY_WORK_DIR}"
        "${CMAKE_COMMAND}" -E env
            # Ensure meson/f2py uses the same Fortran compiler as the CMake build.
            "FC=${CMAKE_Fortran_COMPILER}"
            "CC=${CMAKE_C_COMPILER}"
            "CXX=${CMAKE_CXX_COMPILER}"
            # Ensure generated extension can locate libclubb_f2py_backend.so in
            # the same directory at runtime. Also carry Fortran runtime link paths
            # detected by CMake, since f2py/meson links the extension separately.
            "LDFLAGS=${F2PY_LDFLAGS}"
            "${F2PY_PYTHON_EXECUTABLE}" -m numpy.f2py
            -c "clubb_f2py.pyf"
            # Link against the backend shared library we just built via CMake.
            "-L${F2PY_RUNTIME_DIR}" "-lclubb_f2py_backend"
    DEPENDS clubb_f2py_backend "${F2PY_PYF_FILE}"
    COMMENT "Building clubb_f2py extension via numpy.f2py (reusing CMake-built CLUBB libs)"
    VERBATIM
)

add_custom_command(
    OUTPUT "${F2PY_MODULE_RUNTIME}"
    # Copy extension next to the backend library so the $ORIGIN rpath resolves.
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${F2PY_MODULE_IN_BUILD}" "${F2PY_MODULE_RUNTIME}"
    DEPENDS "${F2PY_MODULE_IN_BUILD}"
    COMMENT "Copying clubb_f2py extension into the build runtime directory"
    VERBATIM
)

# Public target:
# building `clubb_f2py` builds backend + extension and leaves artifacts in the build tree.
# Marked ALL so it runs during normal `cmake --build ...` when ENABLE_F2PY=ON.
add_custom_target(clubb_f2py ALL DEPENDS "${F2PY_MODULE_RUNTIME}")

install(TARGETS clubb_f2py_backend
    LIBRARY DESTINATION "${F2PY_INSTALL_SUBDIR}"
    RUNTIME DESTINATION "${F2PY_INSTALL_SUBDIR}"
)
install(FILES "${F2PY_MODULE_RUNTIME}" DESTINATION "${F2PY_INSTALL_SUBDIR}")
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/clubb_python"
    DESTINATION "${F2PY_INSTALL_SUBDIR}"
    FILES_MATCHING PATTERN "*.py"
    PATTERN "__pycache__" EXCLUDE
)
