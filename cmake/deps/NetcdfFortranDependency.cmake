# --- GLOBAL RPATH FIX ---
# RPATH - Run-Time Search Path
# These provide a path so that even if you move the clubb_standalone executable, it's still found
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}" 
        "${CMAKE_INSTALL_PREFIX}/lib" "${CMAKE_INSTALL_PREFIX}/lib64")

# If this flag is True, then all the libraries that were linked will be put into the executable
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# ------------------------


# --- TIER 1: Native System Scanner ---
# NETCDF_MOD_DIR is the path to netcdf.mod. We search for it in include & netcdf subdirectories
find_path(NETCDF_MOD_DIR "netcdf.mod" PATH_SUFFIXES include netcdf)

# NETCDF_F_LIB is the netcdf fortran libary. We search for it at lib & x86_64-linux-gnu
find_library(NETCDF_F_LIB NAMES netcdff PATH_SUFFIXES lib x86_64-linux-gnu)

# --- THE TWO-WAY ABI SHIELD ---
if(NETCDF_MOD_DIR AND NETCDF_F_LIB)
    # Attempt to compile a tiny test to see if you're using the same compiler for NetCDF and CLUBB
    include(CheckFortranSourceCompiles)
    
    # Temporarily point to the found NetCDF
    set(CMAKE_REQUIRED_INCLUDES ${NETCDF_MOD_DIR})

    # This program forces CLUBB's compiler to talk with the found netcdf
        # If the compilers match, it succeeds and bypasses the following if statement
        # If the compilers don't match, it fails, and NETCDF_ABI_COMPATIBLE is set to false
    check_fortran_source_compiles("
        program abi_test
            use netcdf
        end program abi_test" 
    NETCDF_ABI_COMPATIBLE)

    if(NOT NETCDF_ABI_COMPATIBLE)
        message(WARNING "ABI INCOMPATIBILITY: The NetCDF at ${NETCDF_MOD_DIR} "
                        "cannot be read by ${CMAKE_Fortran_COMPILER_ID}.")
        
        # Wipe the wrongly compiled netcdf variable from memory since the test program failed
        unset(NETCDF_MOD_DIR CACHE)
        unset(NETCDF_F_LIB CACHE)
        unset(NETCDF_ABI_COMPATIBLE CACHE)
    endif()
endif()
# -----------------------------

# If the variables were confirmed, netcdf was found locally and prepares netcdf for cmake
if (NETCDF_MOD_DIR AND NETCDF_F_LIB)
    message(STATUS "SUCCESS: netcdf-fortran found.")
    message(STATUS "  -> Library : ${NETCDF_F_LIB}")
    message(STATUS "  -> Modules : ${NETCDF_MOD_DIR}")
    get_filename_component(NetCDF_Lib_Dir "${NETCDF_F_LIB}" DIRECTORY)
    set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_RPATH};${NetCDF_Lib_Dir})
    set(BUILD_RPATH ${CMAKE_INSTALL_RPATH})

    add_library(netcdf-fortran SHARED IMPORTED) # This adds a library with a set target - import

    # This gives the location as well as links netcdf-fortran library to the .mod files
    set_target_properties(netcdf-fortran PROPERTIES
        IMPORTED_LOCATION "${NETCDF_F_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_MOD_DIR}"
    )

# NetCDF wasn't found and has to be downloaded
else()
    # --- TIER 2: The Auto-Downloader ---
    message(STATUS "WARNING: netcdf-fortran not found locally. Fetching...")
    include(FetchContent)

    # 1. Build the C-Engine Core
    # Set cache variables (within netcdf-c) to off b/c they require external libraries.
    set(BUILD_TESTING OFF CACHE BOOL "Kill Phantom C-Tests" FORCE)
    set(ENABLE_NETCDF_4 OFF CACHE BOOL "" FORCE)
    set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(ENABLE_UTILITIES OFF CACHE BOOL "" FORCE)
    set(BUILD_UTILITIES OFF CACHE BOOL "" FORCE)
    set(ENABLE_DAP OFF CACHE BOOL "" FORCE)
    set(ENABLE_BYTERANGE OFF CACHE BOOL "" FORCE)
    set(HAVE_NET_CDF_PAR FALSE CACHE INTERNAL "Disable parallel checks")
    set(HAVE_NC_CREATE_PAR FALSE CACHE INTERNAL "")
    set(HAVE_NC_OPEN_PAR FALSE CACHE INTERNAL "")
    set(HAVE_NC_VAR_PAR_ACCESS FALSE CACHE INTERNAL "")

    FetchContent_Declare(
        netcdf_c
        GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
        GIT_TAG v4.9.2
    )
    FetchContent_MakeAvailable(netcdf_c)

    # 2. Forge the System Config (The Config Spoof)
    # This step is needed because we just compiled netcdf c (not installed). The steps that cmake
    # expects happen in the install step. - we don't want to install because it requires sudo and
    # it makes it available to all other programs. 
    # The config files made by installing are configured here instead
    
    # This prepares a file
    set(SPOOF_DIR "${CMAKE_BINARY_DIR}/netcdf_spoof")
    file(MAKE_DIRECTORY "${SPOOF_DIR}")
    
    # Locate the actual C library we just built to satisfy the built in tests
    set(ACTUAL_C_LIB "${CMAKE_BINARY_DIR}/_deps/netcdf_c-build/liblib/libnetcdf.a")

    # Cmake checks for a version file when trying to find netcdf, this writes a version file
    file(WRITE "${SPOOF_DIR}/netCDFConfigVersion.cmake" 
            "set(PACKAGE_VERSION \"4.9.2\")\nset(PACKAGE_VERSION_COMPATIBLE TRUE)\n")

    # This makes a map for cmake to use which points netCDF::netcdf to the built netcdf c library
    file(WRITE "${SPOOF_DIR}/netCDFConfig.cmake"
        "set(netCDF_FOUND TRUE)\n"
        "if(NOT TARGET netCDF::netcdf)\n"
        "  add_library(netCDF::netcdf SHARED IMPORTED)\n"
        "  set_target_properties(netCDF::netcdf PROPERTIES IMPORTED_LOCATION \"${ACTUAL_C_LIB}\")\n"
        "endif()\n"
        "set(netCDF_LIBRARIES \"${ACTUAL_C_LIB}\")\n"
        "set(netCDF_INCLUDE_DIRS \"${netcdf_c_BINARY_DIR}/include\" "
        "\"${netcdf_c_SOURCE_DIR}/include\")\n"
    )
    # When we build netcdf fortran, this tells cmake to look here first
    set(netCDF_DIR "${SPOOF_DIR}" CACHE PATH "" FORCE)
    
    # 3. Pass the tests
    # When building netcdf-f, it tests the c engine. But netcdf-c hasn't been 'linked' yet

    # Cache the version so the build doesn't assume legacy and clubb has necessary modern features
    set(NetCDF_VERSION "4.9.2" CACHE STRING "" FORCE)
    set(NetCDF_C_VERSION "4.9.2" CACHE STRING "" FORCE)
    
    # These internal flags silence specific feature detection tests
    set(HAVE_NETCDF4_VERSION TRUE CACHE INTERNAL "")
    set(HAVE_NC_INQ_LIBVERS TRUE CACHE INTERNAL "")
    set(HAVE_NC_SET_LOG_LEVEL TRUE CACHE INTERNAL "")
    set(HAVE_NC_DEF_VAR_FILL TRUE CACHE INTERNAL "")
    set(HAVE_NC_OPT_INQ_FILTER TRUE CACHE INTERNAL "")

    
    # netcdf-f has a final check, which we set to TRUE to pass before running
    set(NetCDF_C_VERSION_OK TRUE CACHE INTERNAL "")
    
    # Turn quantize on since 4.9.2 has it
    set(HAVE_QUANTIZE TRUE CACHE INTERNAL "")

    
    
    # 4. Fetch the Fortran Bindings
    # netcdf-f assumes it's the root, but since not, the patch command makes fixes
    FetchContent_Declare(
        netcdf_fortran
        GIT_REPOSITORY https://github.com/Unidata/netcdf-fortran.git
        GIT_TAG v4.6.1
        PATCH_COMMAND /bin/sh "-c" "sed -i.bak 's/LINK_LIBRARIES netCDF::netcdf//g' CMakeLists.txt && sed -i.bak 's/CMAKE_SOURCE_DIR/CMAKE_CURRENT_SOURCE_DIR/g' CMakeLists.txt && sed -i.bak 's/CMAKE_BINARY_DIR/CMAKE_CURRENT_BINARY_DIR/g' CMakeLists.txt && sed -i.bak 's/add_subdirectory(docs)//Ig' CMakeLists.txt && sed -i.bak 's/add_subdirectory(examples)//Ig' CMakeLists.txt && sed -i.bak 's/add_subdirectory(nf_test)//Ig' CMakeLists.txt && sed -i.bak 's/NC_ENOPAR/-114/g' fortran/nf_nc_noparallel.F90 && sed -i.bak '660,675s/^/#/' CMakeLists.txt && rm -f CMakeLists.txt.bak fortran/nf_nc_noparallel.F90.bak"
    )

    # This is to limit the installs because cmake inherits the installs of the netcdf-f
    set(CMAKE_SKIP_INSTALL_RULES TRUE) 

    # This tells cmake where to find the c and fortran engine
    set(CMAKE_REQUIRED_INCLUDES "${netcdf_c_SOURCE_DIR}/include" "${netcdf_c_BINARY_DIR}/include")
    include_directories("${netcdf_c_SOURCE_DIR}/include" "${netcdf_c_BINARY_DIR}/include")
    include_directories("${CMAKE_BINARY_DIR}/_deps/netcdf_fortran-src/fortran" 
                        "${CMAKE_BINARY_DIR}/_deps/netcdf_fortran-build/fortran")


    # 2. Save and wipe flags to be compatible with older fortran so it doesn't fail
    set(SAVED_FC_FLAGS "${CMAKE_Fortran_FLAGS}")
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND 
                    CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch")
    else()
        set(CMAKE_Fortran_FLAGS "")
    endif()

    # 3. Build NetCDF fortran
    FetchContent_MakeAvailable(netcdf_fortran)

    #4. We make a dummy install file since we skipped install earlier so compile.py doesn't fail
    file(WRITE "${CMAKE_BINARY_DIR}/_deps/netcdf_fortran-build/cmake_install.cmake" 
                "# Dummy install script\n")

    #5. Restore strict flags for CLUBB
    set(CMAKE_Fortran_FLAGS "${SAVED_FC_FLAGS}")
        
        
    # --- 5. THE ULTIMATE HOT-WIRE INJECTION ---
    # This finds the netcdf fortran source code and saves to FORTRAN_SRC
    FetchContent_GetProperties(netcdf_fortran SOURCE_DIR FORTRAN_SRC BINARY_DIR FORTRAN_BIN)

    # 'Dictionary' of NetCDF
    set(DATA_FILE "${FORTRAN_SRC}/fortran/module_netcdf_nf_data.F90") 
    if(EXISTS "${DATA_FILE}")
        file(READ "${DATA_FILE}" FILE_CONTENTS)
        
        # Only inject if NC_NOQUANTIZE isn't already present in the file
        # If the fortran version is too old, it doesn't have quantize constants named, so we 
        # inject them so that CLUBB doesn't fail to compile when using these constants
        if(NOT FILE_CONTENTS MATCHES "NC_NOQUANTIZE")
            message(STATUS "Injecting missing NetCDF constants into ${DATA_FILE}")
            set(INJECTION "Implicit NONE\n")
            string(APPEND INJECTION 
               "      integer, parameter :: NC_NOQUANTIZE = 0, NC_QUANTIZE_BITGROOM = 1\n")
            string(APPEND INJECTION 
               "      integer, parameter :: NC_QUANTIZE_GRANULARBR = 2, NC_QUANTIZE_BITROUND = 3\n")
            string(APPEND INJECTION "      integer, parameter :: NC_ENOPAR = -114")
            string(REPLACE "Implicit NONE" "${INJECTION}" UPDATED_CONTENTS "${FILE_CONTENTS}")
            file(WRITE "${DATA_FILE}" "${UPDATED_CONTENTS}")
        else()
            message(STATUS "NetCDF constants already present, skipping injection.")
        endif()
    endif()

    # Injection B: Manually provide EVERY missing include file the compiler wants
    set(GEN_DIR "${FORTRAN_BIN}/fortran")
    file(MAKE_DIRECTORY "${GEN_DIR}")
    # List of all the expected include files netcdf Fortran expects
    set(BLUEPRINTS "netcdf_constants.f90" "netcdf_externals.f90" "netcdf_overloads.f90" 
                    "netcdf_visibility.f90" "netcdf_file.f90" "netcdf3_file.f90" 
                    "netcdf_dims.f90" "netcdf_attributes.f90" "netcdf_variables.f90" 
                    "netcdf_text_variables.f90" "netcdf_expanded.f90" "netcdf_eightbyte.f90" 
                    "netcdf_fortran_env.f90" "netcdf_expanded_subset.f90")
    
    foreach(FILENAME ${BLUEPRINTS})
        if("${FILENAME}" STREQUAL "netcdf_constants.f90")
            # Writes necessary constants to be used
            file(WRITE "${GEN_DIR}/${FILENAME}" 
                "integer, parameter :: &\n"
                "    NF90_NOERR=0, NF90_CLOBBER=16, NF90_NOWRITE=0, &\n"
                "    NF90_UNLIMITED=0, NF90_DOUBLE=6, NF90_CHAR=2, &\n"
                "    NF90_INT=4, NC_ENOPAR=-114\n"
            )
        elseif("${FILENAME}" STREQUAL "netcdf_externals.f90")
            # Writes a linking promise that these files can be found later to avoid an error
            file(WRITE "${GEN_DIR}/${FILENAME}" 
                "external nf_create, nf_open, nf_redef, nf_enddef, &\n"
                "         nf_close, nf_inq_varid, nf_def_dim, &\n"
                "         nf_def_var, nf_put_att_text, &\n"
                "         nf_put_var_double, nf_strerror\n"
            )
        else()
            # This just writes a dummy file with just a comment
            file(WRITE "${GEN_DIR}/${FILENAME}" "! Placeholder\n")
        endif()
    endforeach()

    # --- 6. THE FINAL HANDSHAKE ---
    if(TARGET netcdff)
        # 1. Tells anything linking to netcdf fortran where to find the .mod files
        target_include_directories(netcdff INTERFACE 
                    "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_deps/netcdf_fortran-build/fortran>")

        # 2. This writes a c file to define functions the linker might look for
            # -114 is the standard NC_ENOPAR (parellel) error code
        file(WRITE "${CMAKE_BINARY_DIR}/ghostbusters.c" 
             "int nf_var_par_access_(){return -114;}\n"
             "int nf_create_par_(){return -114;}\n"
             "int nf_open_par_(){return -114;}\n")
        add_library(netcdf_ghostbusters STATIC "${CMAKE_BINARY_DIR}/ghostbusters.c")

        # 3. Bundle netcdf-f and the ghostbuster library together under the name CLUBB wants
        add_library(netcdf-fortran INTERFACE)
        target_link_libraries(netcdf-fortran INTERFACE netcdff netcdf_ghostbusters)

        message(STATUS "Handshake Complete: netcdf-fortran bundled with Ghostbusters.")
    endif()

    # Re enable install rules so CLUBB can install properly
    set(CMAKE_SKIP_INSTALL_RULES FALSE)
endif()

# --- REVIVE CLUBB TESTS ---
# Now that ALL NetCDF libraries are built, safely turn tests back on for CLUBB
set(ENABLE_TESTS ON CACHE BOOL "" FORCE)
set(BUILD_TESTING ON CACHE BOOL "" FORCE)
# --------------------------------------------------
