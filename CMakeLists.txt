# This is the top level CMakeLists.txt in clubb
# Initial issue: https://github.com/larson-group/clubb/issues/1117
# Introduced in PR: https://github.com/larson-group/clubb/pull/1253
# WIKI page: https://github.com/larson-group/clubb/wiki/CMake

cmake_minimum_required(VERSION 3.18)

# Define Build system options
include(CMakeDependentOption)
option(USE_NetCDF   "Utilize netcdf for model output data format"       ON  )   # netcdf is required for output
option(TUNING       "Compile additional code useful for TUNING runs"    OFF )   # enables some extra runtime checks
option(SILHS        "Enable SILHS"                                      ON  )   # silhs is default in some cases
option(ENABLE_OMP   "Enable CLUBB to run with OpenMP"                   OFF )   # enables openmp
option(ENABLE_TESTS "Enable creation of ctests"                         ON  )   # enable quick run tests
option(USE_GPTL     "Enables use of the gptl timing library"            OFF )   # a library with nice timing features

# Create the project
project(CLUBB Fortran C CXX)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/deps;${CMAKE_CURRENT_LIST_DIR}/cmake/utils;${CMAKE_MODULE_PATH})
set(CLUBB_ROOT ${CMAKE_CURRENT_LIST_DIR})


# Handle our CPPDefinitions such as CLUBB_REAL_TYPE
include(CPPDefinitions)

if(${USE_NetCDF})
    # netCDF is an external library, cmake/deps/NetcdfFortranDependency.cmake handles making netcdf available
    # and enables linking via 'netcdf-fortran'
    include(NetcdfFortranDependency)
endif()

if(${ENABLE_OMP})
    # cmake/deps/OpenMPDependency.cmake finds openmp and makes it available via 'openmp-fortran'
    include(OpenMPDependency)
endif()

if(${USE_GPTL})
    # cmake/deps/GPTLDependency.cmake finds gptl and makes it available via 'gptl'
    include(GPTLDependency)
endif()

# Add libraries/executables
add_subdirectory(src)


# These define ctests, which as currenlty implemented are intended as a quick set of 
# tests just to ensure the code runs with different options and cases.
# These tests can can be run by navigating to the build directory and running "ctest"
# or these can be ran automatically by calling "./compile.py -run_tests"
if(${ENABLE_TESTS})

    include(CTest)

    set(TEST_OUT_DIR ${CMAKE_BINARY_DIR}/clubb_ctest_output)

    # One-time setup before any tests that require the fixture
    add_test(NAME _init_output_dir
            COMMAND bash -c
            "${CMAKE_COMMAND} -E rm -rf ${TEST_OUT_DIR} && ${CMAKE_COMMAND} -E make_directory ${TEST_OUT_DIR}")
    set_tests_properties(_init_output_dir PROPERTIES FIXTURES_SETUP CleanOutputDir)

    # Function to add a case without wiping the whole output dir each time
    function(run_clubb_case TEST_NAME)
        set(RUN_COMMAND
            "${CMAKE_COMMAND} -E make_directory ${TEST_OUT_DIR}/${TEST_NAME} && \
            ${Python_EXECUTABLE} ${CLUBB_ROOT}/run_scripts/run_scm.py \
            -out_dir ${TEST_OUT_DIR}/${TEST_NAME} ${ARGN}"
        )

        add_test(NAME ${TEST_NAME}
                COMMAND bash -c "${RUN_COMMAND}"
                WORKING_DIRECTORY ${CLUBB_ROOT}/run_scripts)

        # Ensure the one-time init runs before any case
        set_tests_properties(${TEST_NAME} PROPERTIES FIXTURES_REQUIRED CleanOutputDir)
    endfunction()


    # Run 100 iterations of a handful of common casse
    run_clubb_case( bomex_100iter "-max_iters 100 bomex" )
    run_clubb_case( wangara_100iter "-max_iters 100 wangara" )
    run_clubb_case( fire_100iter "-max_iters 100 fire" )
    run_clubb_case( atex_100iter "-max_iters 100 atex" )

    # rico_silhs to test silhs_lib
    run_clubb_case( rico_silhs_100iter "-max_iters 100 -override l_lh_straight_mc=.true. rico_silhs" )

    # arm_97 to test microphys and SILHS
    run_clubb_case( arm_97_50iter "-max_iters 50 -override l_lh_straight_mc=.true. arm_97" )

    # cobra to test scalars and ice code
    run_clubb_case( cobra_100iter "-max_iters 100 cobra" )

    # test for r8029 parameter/flag set
    set(R8029_CONFIG ${CLUBB_ROOT}/input/tunable_parameters_compatible_r8029/)
    run_clubb_case( arm_r8029_100iter "-max_iters 100 -config ${R8029_CONFIG} arm" )

    # multicol test and disable normal stats output
    run_clubb_case( arm_100iter_8col "-max_iters 100 -ngrdcol 8 -tout 0 arm" )

    # test with small timestep to check that it doesn't break anything
    run_clubb_case( arm_1000iter_5dt "-max_iters 1000 -dt_main 5 arm" )

endif()

