# -------------------------------------------------------------------------
# GPTL Dependency Handler for CLUBB
# -------------------------------------------------------------------------
# - Searches for existing GPTL installation first
# - If not found, downloads and builds GPTL v8.1.1 from GitHub
# - Supports Fortran-enabled build
# - Ensures GPTL is built *before* any CLUBB targets link to it
# -------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)

include(ExternalProject)

# -----------------------------------------------------------------------------
# 1. Try to find GPTL in environment or provided root
# -----------------------------------------------------------------------------
find_library(GPTLF_LIBRARY NAMES gptlf HINTS ${GPTL_ROOT}/lib $ENV{GPTL_ROOT}/lib)
find_library(GPTL_LIBRARY  NAMES gptl  HINTS ${GPTL_ROOT}/lib $ENV{GPTL_ROOT}/lib)

find_path(GPTL_INCLUDE_DIR
  NAMES gptl.mod gptl.inc GPTL.h
  HINTS
    ${GPTL_ROOT}/include $ENV{GPTL_ROOT}/include
    ${GPTL_ROOT}/mod $ENV{GPTL_ROOT}/mod
  PATH_SUFFIXES include mod include/nvfortran mod/nvfortran
)

set(GPTL_LIBRARIES "")
if(GPTLF_LIBRARY)
  list(APPEND GPTL_LIBRARIES ${GPTLF_LIBRARY})
endif()
if(GPTL_LIBRARY)
  list(APPEND GPTL_LIBRARIES ${GPTL_LIBRARY})
endif()

# -----------------------------------------------------------------------------
# 2. If GPTL is found, create imported target and return
# -----------------------------------------------------------------------------
if(GPTL_LIBRARIES AND GPTL_INCLUDE_DIR)
  message(STATUS "Found GPTL installation:")
  message(STATUS "  Includes: ${GPTL_INCLUDE_DIR}")
  message(STATUS "  Libraries: ${GPTL_LIBRARIES}")

  add_library(gptl INTERFACE IMPORTED)
  set(GPTL_ALL_LIBS ${GPTL_LIBRARIES})

  # Try to link libunwind if it exists
  find_library(UNWIND_LIBRARY NAMES unwind unwind-x86_64)
  if(UNWIND_LIBRARY)
    list(APPEND GPTL_ALL_LIBS ${UNWIND_LIBRARY})
  endif()

  set_target_properties(gptl PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${GPTL_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${GPTL_ALL_LIBS}"
  )
  return()
endif()

# -----------------------------------------------------------------------------
# 3. Otherwise, download and build GPTL v8.1.1
# -----------------------------------------------------------------------------
message(WARNING "GPTL not found — downloading and building GPTL v8.1.1...")

set(GPTL_PREFIX ${CMAKE_BINARY_DIR}/gptl_install)
set(GPTL_SRC    ${CMAKE_BINARY_DIR}/gptl_src)
file(MAKE_DIRECTORY ${GPTL_PREFIX}/include)
file(MAKE_DIRECTORY ${GPTL_PREFIX}/lib)

set(GPTL_ENV_VARS
  "CC=${CMAKE_C_COMPILER}"
  "FC=${CMAKE_Fortran_COMPILER}"
  "CFLAGS=${CMAKE_C_FLAGS}"
  "FFLAGS=${CMAKE_Fortran_FLAGS}"
)

ExternalProject_Add(GPTL_project
  URL https://github.com/jmrosinski/GPTL/releases/download/v8.1.1/gptl-8.1.1.tar.gz
  PREFIX ${CMAKE_BINARY_DIR}/gptl_build
  SOURCE_DIR ${GPTL_SRC}
  CONFIGURE_COMMAND /bin/sh -c "
    # Prevent autoheader/autoconf from running if timestamps confuse 'missing'
    touch aclocal.m4 configure Makefile.in fortran/Makefile.in config.h.in &&
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} \
    ./configure --prefix=${GPTL_PREFIX} --enable-fortran"
  BUILD_COMMAND /bin/sh -c "
    # Remove problematic instrumentation flag and build
    find . -name Makefile | xargs sed -i 's/-finstrument-functions//g';
    make -j"
  INSTALL_COMMAND make install
  BUILD_IN_SOURCE 1
  LOG_CONFIGURE 1
  LOG_BUILD 1
  LOG_INSTALL 1
)

# -----------------------------------------------------------------------------
# 4. Create a "stamp" target so Ninja knows these files are generated
# -----------------------------------------------------------------------------
add_custom_target(gptl_built
  COMMAND ${CMAKE_COMMAND} -E echo "GPTL static libraries ready."
  DEPENDS GPTL_project
  BYPRODUCTS
    ${GPTL_PREFIX}/lib/libgptl.a
    ${GPTL_PREFIX}/lib/libgptlf.a
  COMMENT "Waiting for GPTL static libraries to be built"
)

# -----------------------------------------------------------------------------
# 5. Create imported GPTL target for CLUBB to link to
# -----------------------------------------------------------------------------
add_library(gptl INTERFACE IMPORTED)
add_dependencies(gptl gptl_built)

set(GPTL_INCLUDE_DIR ${GPTL_PREFIX}/include)
set(GPTL_LIBRARIES
  ${GPTL_PREFIX}/lib/libgptlf.a
  ${GPTL_PREFIX}/lib/libgptl.a
)

set_target_properties(gptl PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${GPTL_INCLUDE_DIR}"
  INTERFACE_LINK_LIBRARIES "${GPTL_LIBRARIES}"
)

# -----------------------------------------------------------------------------
# 6. Link libunwind if available (used by GPTL for callstack tracking)
# -----------------------------------------------------------------------------
find_package(PkgConfig QUIET)
pkg_check_modules(UNWIND QUIET libunwind)

if(UNWIND_FOUND)
  message(STATUS "Found libunwind via pkg-config: ${UNWIND_LIBRARIES}")
  set(UNWIND_LIBRARY ${UNWIND_LIBRARIES})
else()
  message(STATUS "pkg-config could not locate libunwind — trying manual search...")
  find_library(UNWIND_LIBRARY
    NAMES unwind unwind-x86_64
    HINTS
      /usr/lib
      /usr/local/lib
      /usr/lib64
      /usr/local/lib64
      /usr/lib/x86_64-linux-gnu
      /opt/local/lib
  )
endif()

if(UNWIND_LIBRARY)
  message(STATUS "Found libunwind: ${UNWIND_LIBRARY}")
  set_property(TARGET gptl APPEND PROPERTY
    INTERFACE_LINK_LIBRARIES "${UNWIND_LIBRARY}"
  )
else()
  message(WARNING "libunwind not found — GPTL callstack functionality may be disabled")
endif()
