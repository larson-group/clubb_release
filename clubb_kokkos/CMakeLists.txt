cmake_minimum_required(VERSION 3.18)

project(CLUBB_Kokkos LANGUAGES C CXX Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build configuration" FORCE)
endif()

# Keep Native comparison targets and optional PDF matrix calls in this checkout.
get_filename_component(CLUBB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(CMAKE_MODULE_PATH
    "${CLUBB_ROOT}/cmake/deps;${CLUBB_ROOT}/cmake/utils;${CMAKE_MODULE_PATH}")

option(TUNING "Compile additional tuning checks" OFF)
option(ENABLE_OMP "Enable Fortran OpenMP" OFF)
option(USE_GPTL "Enable GPTL timing" OFF)
option(SILHS "Enable SILHS" ON)
option(KOKKOS_BFB "Use Fortran/f2k SILHS matrix operations for BFB comparisons" OFF)

include(CPPDefinitions)
include(CompilerWarnings)
include(NetcdfFortranDependency)
if(ENABLE_OMP)
    include(OpenMPDependency)
endif()
if(USE_GPTL)
    include(GPTLDependency)
endif()

add_subdirectory("${CLUBB_ROOT}/src" "${CMAKE_BINARY_DIR}/fortran" EXCLUDE_FROM_ALL)

find_package(Kokkos CONFIG REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(NETCDF_C REQUIRED IMPORTED_TARGET netcdf)

add_library(clubb_kokkos_error_code_lib STATIC src/CLUBB_core/error_code.cpp)
target_compile_features(clubb_kokkos_error_code_lib PUBLIC cxx_std_17)
target_include_directories(clubb_kokkos_error_code_lib PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/CLUBB_core")
target_link_libraries(clubb_kokkos_error_code_lib PUBLIC Kokkos::kokkos)

add_library(clubb_kokkos_timer_lib STATIC src/CLUBB_core/code_timer_module.cpp)
target_compile_features(clubb_kokkos_timer_lib PUBLIC cxx_std_17)
target_include_directories(clubb_kokkos_timer_lib PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/CLUBB_core")
if(USE_GPTL)
    if(TARGET GPTL_project)
        # The existing dependency installs this C header before compiling timers.
        set(CLUBB_KOKKOS_GPTL_C_INCLUDE_DIR "${GPTL_INCLUDE_DIR}")
    else()
        find_path(CLUBB_KOKKOS_GPTL_C_INCLUDE_DIR NAMES gptl.h
            HINTS "${GPTL_INCLUDE_DIR}" "${GPTL_ROOT}/include" "$ENV{GPTL_ROOT}/include"
                  "${GPTL_INCLUDE_DIR}/..")
        if(NOT CLUBB_KOKKOS_GPTL_C_INCLUDE_DIR)
            message(FATAL_ERROR "C++ GPTL timers require gptl.h; set GPTL_ROOT to its installation")
        endif()
    endif()
    target_include_directories(clubb_kokkos_timer_lib PRIVATE "${CLUBB_KOKKOS_GPTL_C_INCLUDE_DIR}")
    target_link_libraries(clubb_kokkos_timer_lib PRIVATE gptl)
endif()

add_library(clubb_kokkos_input_io_lib STATIC
    src/fortran_io.cpp
    src/namelist.cpp
    src/netcdf.cpp
    src/CLUBB_core/file_functions.cpp)
target_compile_features(clubb_kokkos_input_io_lib PUBLIC cxx_std_20)
target_include_directories(clubb_kokkos_input_io_lib PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/CLUBB_core")
target_link_libraries(clubb_kokkos_input_io_lib PUBLIC Kokkos::kokkos PkgConfig::NETCDF_C)

# The remaining production bridge is only the optional PDF/LAPACK ABI.
if(KOKKOS_BFB)
    add_library(clubb_f2k_lib STATIC src/f2k/lapack_bridge.F90)
    add_dependencies(clubb_f2k_lib clubb_core_lib)
    target_link_libraries(clubb_f2k_lib PRIVATE lapack_lib)
endif()

add_library(clubb_kokkos_core_lib STATIC
    src/CLUBB_core/advance_clubb_core_module.cpp
    src/CLUBB_core/stats_netcdf.cpp
    src/CLUBB_core/stats_clubb_utilities.cpp
    src/CLUBB_core/corr_varnce_module.cpp
    src/CLUBB_core/precipitation_fraction.cpp
    src/CLUBB_core/setup_clubb_pdf_params.cpp
    src/CLUBB_core/advance_xp3_module.cpp
    src/CLUBB_core/advance_xm_wpxp_module.cpp
    src/CLUBB_core/advance_xp2_xpyp_module.cpp
    src/CLUBB_core/advance_wp2_wp3_module.cpp
    src/CLUBB_core/advance_windm_edsclrm_module.cpp
    src/CLUBB_core/diffusion.cpp
    src/CLUBB_core/mean_adv.cpp
    src/CLUBB_core/turbulent_adv_pdf.cpp
    src/CLUBB_core/pos_definite_module.cpp
    src/CLUBB_core/clip_explicit.cpp
    src/CLUBB_core/mono_flux_limiter.cpp
    src/CLUBB_core/fill_holes.cpp
    src/CLUBB_core/index_mapping.cpp
    src/CLUBB_core/Skx_module.cpp
    src/CLUBB_core/sigma_sqd_w_module.cpp
    src/CLUBB_core/tridiag_lu_solver.cpp
    src/CLUBB_core/penta_lu_solver.cpp
    src/CLUBB_core/lapack_wrap.cpp
    src/CLUBB_core/matrix_solver_wrapper.cpp
    src/CLUBB_core/err_info_type_module.cpp
    src/CLUBB_core/hydromet_pdf_parameter_module.cpp
    src/CLUBB_core/pdf_parameter_module.cpp
    src/CLUBB_core/calendar.cpp
    src/CLUBB_core/remapping_module.cpp
    src/CLUBB_core/model_flags.cpp
    src/CLUBB_core/parameters_tunable.cpp
    src/CLUBB_core/grid_class.cpp
    src/CLUBB_core/mixing_length.cpp
    src/CLUBB_core/sfc_varnce_module.cpp
    src/CLUBB_core/numerical_check.cpp
    src/CLUBB_core/pdf_closure_module.cpp
    src/CLUBB_core/adg1_adg2_3d_luhar_pdf.cpp
    src/CLUBB_core/new_pdf_main.cpp
    src/CLUBB_core/new_hybrid_pdf_main.cpp)
target_compile_features(clubb_kokkos_core_lib PRIVATE cxx_std_17)
target_link_libraries(clubb_kokkos_core_lib PUBLIC
    clubb_kokkos_input_io_lib clubb_kokkos_timer_lib clubb_kokkos_error_code_lib)

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/KokkosBfb.cmake")
if(KOKKOS_BFB_ENABLED)
    target_link_libraries(clubb_kokkos_core_lib PRIVATE clubb_f2k_lib)
    target_compile_definitions(clubb_kokkos_core_lib PUBLIC KOKKOS_BFB=1)
endif()

# Keep execution backends and requested/effective BFB modes with test results.
list(JOIN Kokkos_DEVICES "\", \"" clubb_kokkos_devices_json)
set(clubb_kokkos_bfb_requested_json false)
set(clubb_kokkos_bfb_json false)
if(KOKKOS_BFB)
    set(clubb_kokkos_bfb_requested_json true)
endif()
if(KOKKOS_BFB_ENABLED)
    set(clubb_kokkos_bfb_json true)
endif()
file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/kokkos_backends.json"
    CONTENT "{\n  \"devices\": [\"${clubb_kokkos_devices_json}\"],\n  \"kokkos_bfb_requested\": ${clubb_kokkos_bfb_requested_json},\n  \"kokkos_bfb\": ${clubb_kokkos_bfb_json}\n}\n")

add_library(clubb_kokkos_driver_lib STATIC
    src/Radiation/radiation_io.cpp
    src/clubb_driver.cpp
    src/Input_fields/hydrostatic_module.cpp
    src/Input_fields/input_fields.cpp
    src/Input_fields/input_interpret.cpp
    src/Input_fields/input_netcdf.cpp
    src/Input_fields/input_reader.cpp
    src/Input_fields/sounding.cpp
    src/Input_fields/stat_file_utils.cpp
    src/Benchmark_cases/mpace_a.cpp
    src/Benchmark_cases/prescribe_forcings.cpp
    src/Benchmark_cases/time_dependent_input.cpp
    src/Radiation/cos_solar_zen_module.cpp
    src/Radiation/extended_atmosphere_module.cpp
    src/Radiation/parameters_radiation.cpp
    src/Radiation/rad_lwsw_module.cpp
    src/Radiation/radiation_module.cpp
    src/Radiation/radiation_variables_module.cpp
    src/Radiation/simple_rad_module.cpp
    src/Radiation/soil_vegetation.cpp
    src/SILHS/est_kessler_microphys_module.cpp
    src/SILHS/generate_uniform_sample_module.cpp
    src/SILHS/latin_hypercube_arrays.cpp
    src/SILHS/latin_hypercube_driver_module.cpp
    src/SILHS/lh_microphys_var_covar_module.cpp
    src/SILHS/math_utilities.cpp
    src/SILHS/mt95.cpp
    src/SILHS/parameters_silhs.cpp
    src/SILHS/silhs_api_module.cpp
    src/SILHS/silhs_importance_sample_module.cpp
    src/SILHS/transform_to_pdf_module.cpp
    src/Microphys/cloud_sed_module.cpp
    src/readiopdata_core/readiopdata_core_module.cpp
    src/readiopdata_module.cpp)
target_compile_features(clubb_kokkos_driver_lib PUBLIC cxx_std_20)
target_link_libraries(clubb_kokkos_driver_lib PUBLIC clubb_kokkos_core_lib)

# src/ also defines the native clubb_standalone target; keep the public filename.
add_executable(clubb_kokkos_standalone src/clubb_standalone.cpp)
set_target_properties(clubb_kokkos_standalone PROPERTIES OUTPUT_NAME clubb_standalone)
target_compile_features(clubb_kokkos_standalone PRIVATE cxx_std_17)
target_link_libraries(clubb_kokkos_standalone PRIVATE clubb_kokkos_driver_lib)

include(CTest)
if(BUILD_TESTING)
    find_package(Threads REQUIRED)
    # Native derived-type conversion and numerical oracles are test-only.
    add_library(clubb_kokkos_native_test_lib STATIC
        src/f2k/derived_type_bridge.cpp
        src/f2k/derived_type_bridge.F90)
    target_compile_features(clubb_kokkos_native_test_lib PRIVATE cxx_std_17)
    target_include_directories(clubb_kokkos_native_test_lib PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(clubb_kokkos_native_test_lib PUBLIC Kokkos::kokkos PRIVATE
        clubb_core_lib input_fields_lib radiation_lib silhs_lib
        benchmark_cases_lib microphys_lib text_writer_lib netcdf-fortran)
    add_executable(clubb_kokkos_debug_level
        auto_llm_generated_tests/debug_level.cpp
        auto_llm_generated_tests/debug_level_verify.F90)
    target_link_libraries(clubb_kokkos_debug_level PRIVATE
        clubb_kokkos_error_code_lib clubb_core_lib Threads::Threads)
    add_test(NAME clubb_kokkos_debug_level COMMAND clubb_kokkos_debug_level)

    add_executable(clubb_kokkos_code_timer auto_llm_generated_tests/code_timer.cpp)
    target_link_libraries(clubb_kokkos_code_timer PRIVATE clubb_kokkos_timer_lib)
    add_test(NAME clubb_kokkos_code_timer COMMAND clubb_kokkos_code_timer
        "${CMAKE_CURRENT_BINARY_DIR}/code_timer_test_data")

    add_executable(clubb_kokkos_core_grid_operators
        auto_llm_generated_tests/core_grid_operators.cpp
        auto_llm_generated_tests/core_grid_operators_verify.F90)
    target_compile_features(clubb_kokkos_core_grid_operators PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_grid_operators PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_core_grid_operators COMMAND clubb_kokkos_core_grid_operators)

    add_executable(clubb_kokkos_core_interpolation
        auto_llm_generated_tests/core_interpolation.cpp
        auto_llm_generated_tests/core_interpolation_verify.F90)
    target_compile_features(clubb_kokkos_core_interpolation PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_interpolation PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_core_interpolation COMMAND clubb_kokkos_core_interpolation)

    add_executable(clubb_kokkos_core_matrix_wrapper
        auto_llm_generated_tests/core_matrix_wrapper.cpp
        auto_llm_generated_tests/core_matrix_wrapper_verify.F90)
    target_compile_features(clubb_kokkos_core_matrix_wrapper PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_matrix_wrapper PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_core_matrix_wrapper COMMAND clubb_kokkos_core_matrix_wrapper)

    add_executable(clubb_kokkos_core_lapack_diagnostics
        auto_llm_generated_tests/core_lapack_diagnostics.cpp)
    target_compile_features(clubb_kokkos_core_lapack_diagnostics PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_lapack_diagnostics PRIVATE clubb_kokkos_core_lib)
    add_test(NAME clubb_kokkos_core_lapack_diagnostics
        COMMAND clubb_kokkos_core_lapack_diagnostics)

    add_executable(clubb_kokkos_core_cholesky
        auto_llm_generated_tests/core_cholesky.cpp)
    target_compile_features(clubb_kokkos_core_cholesky PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_cholesky PRIVATE clubb_kokkos_core_lib)
    add_test(NAME clubb_kokkos_core_cholesky COMMAND clubb_kokkos_core_cholesky)

    add_executable(clubb_kokkos_namelist auto_llm_generated_tests/namelist.cpp)
    target_link_libraries(clubb_kokkos_namelist PRIVATE clubb_kokkos_input_io_lib)
    add_test(NAME clubb_kokkos_namelist COMMAND clubb_kokkos_namelist)

    add_executable(clubb_kokkos_netcdf_input auto_llm_generated_tests/netcdf_input.cpp)
    target_link_libraries(clubb_kokkos_netcdf_input PRIVATE clubb_kokkos_driver_lib)
    add_test(NAME clubb_kokkos_netcdf_input COMMAND clubb_kokkos_netcdf_input)

    add_executable(clubb_kokkos_input_text_io auto_llm_generated_tests/input_text_io.cpp)
    target_link_libraries(clubb_kokkos_input_text_io PRIVATE clubb_kokkos_driver_lib)
    add_test(NAME clubb_kokkos_input_text_io COMMAND clubb_kokkos_input_text_io
        "${CMAKE_CURRENT_BINARY_DIR}/input_text_io_data")

    add_executable(clubb_kokkos_core_lu_solvers
        auto_llm_generated_tests/core_lu_solvers.cpp
        auto_llm_generated_tests/core_lu_solvers_verify.F90)
    target_compile_features(clubb_kokkos_core_lu_solvers PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_lu_solvers PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_core_lu_solvers COMMAND clubb_kokkos_core_lu_solvers)

    add_executable(clubb_kokkos_core_pdf_kernels
        auto_llm_generated_tests/core_pdf_kernels.cpp
        auto_llm_generated_tests/core_pdf_kernels_verify.F90)
    target_compile_features(clubb_kokkos_core_pdf_kernels PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_pdf_kernels PRIVATE
        clubb_kokkos_core_lib clubb_core_lib)
    add_test(NAME clubb_kokkos_core_pdf_kernels COMMAND clubb_kokkos_core_pdf_kernels)

    add_executable(clubb_kokkos_radiation_bridge
        auto_llm_generated_tests/radiation_bridge.cpp)
    target_compile_features(clubb_kokkos_radiation_bridge PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_radiation_bridge PRIVATE clubb_kokkos_driver_lib)
    add_test(NAME clubb_kokkos_radiation_bridge COMMAND clubb_kokkos_radiation_bridge)

    add_executable(clubb_kokkos_core_error_return
        auto_llm_generated_tests/core_error_return.cpp)
    target_compile_features(clubb_kokkos_core_error_return PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_error_return PRIVATE clubb_kokkos_core_lib)
    add_test(NAME clubb_kokkos_core_error_return COMMAND clubb_kokkos_core_error_return)

    add_executable(clubb_kokkos_config_error_output
        auto_llm_generated_tests/config_error_output.cpp
        auto_llm_generated_tests/config_error_output_verify.F90)
    target_link_libraries(clubb_kokkos_config_error_output PRIVATE
        clubb_kokkos_driver_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_config_error_output COMMAND clubb_kokkos_config_error_output)

    add_executable(clubb_kokkos_core_calendar
        auto_llm_generated_tests/core_calendar.cpp)
    target_compile_features(clubb_kokkos_core_calendar PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_calendar PRIVATE clubb_kokkos_core_lib)
    add_test(NAME clubb_kokkos_core_calendar COMMAND clubb_kokkos_core_calendar)

    add_executable(clubb_kokkos_core_parameter_validation
        auto_llm_generated_tests/core_parameter_validation.cpp
        auto_llm_generated_tests/core_parameter_validation_verify.F90)
    target_compile_features(clubb_kokkos_core_parameter_validation PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_parameter_validation PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_core_parameter_validation
        COMMAND clubb_kokkos_core_parameter_validation)

    add_executable(clubb_kokkos_core_array_transfer
        auto_llm_generated_tests/core_array_transfer.cpp
        auto_llm_generated_tests/core_array_transfer_verify.F90)
    target_compile_features(clubb_kokkos_core_array_transfer PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_core_array_transfer PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_core_array_transfer COMMAND clubb_kokkos_core_array_transfer)

    add_executable(clubb_kokkos_derived_type_roundtrip
        auto_llm_generated_tests/derived_type_roundtrip.cpp)
    target_compile_features(clubb_kokkos_derived_type_roundtrip PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_derived_type_roundtrip PRIVATE
        clubb_kokkos_core_lib clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_derived_type_roundtrip
        COMMAND clubb_kokkos_derived_type_roundtrip)

    add_executable(clubb_kokkos_hydromet_pdf_parameter_collection
        auto_llm_generated_tests/hydromet_pdf_parameter_collection.cpp)
    target_compile_features(clubb_kokkos_hydromet_pdf_parameter_collection PRIVATE cxx_std_17)
    target_include_directories(clubb_kokkos_hydromet_pdf_parameter_collection PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(clubb_kokkos_hydromet_pdf_parameter_collection PRIVATE clubb_kokkos_native_test_lib)
    add_test(NAME clubb_kokkos_hydromet_pdf_parameter_collection
        COMMAND clubb_kokkos_hydromet_pdf_parameter_collection)

    # This case test accepts a generated aggregate namelist through run_scm.py -exe.
    add_executable(clubb_kokkos_driver_sequence
        auto_llm_generated_tests/clubb_driver_sequence.cpp)
    target_compile_features(clubb_kokkos_driver_sequence PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_driver_sequence PRIVATE clubb_kokkos_driver_lib)

    add_executable(clubb_kokkos_driver_data
        auto_llm_generated_tests/driver_data.cpp)
    target_compile_features(clubb_kokkos_driver_data PRIVATE cxx_std_17)
    target_include_directories(clubb_kokkos_driver_data PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(clubb_kokkos_driver_data PRIVATE Kokkos::kokkos)
    add_test(NAME clubb_kokkos_driver_data COMMAND clubb_kokkos_driver_data)

    add_executable(clubb_kokkos_scratch_workspace
        auto_llm_generated_tests/scratch_workspace.cpp)
    target_compile_features(clubb_kokkos_scratch_workspace PRIVATE cxx_std_20)
    target_include_directories(clubb_kokkos_scratch_workspace PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/src")
    target_link_libraries(clubb_kokkos_scratch_workspace PRIVATE Kokkos::kokkos)
    add_test(NAME clubb_kokkos_scratch_workspace COMMAND clubb_kokkos_scratch_workspace)

    add_executable(clubb_kokkos_stats_netcdf
        auto_llm_generated_tests/stats_netcdf.cpp
        auto_llm_generated_tests/stats_netcdf_verify.F90)
    target_compile_features(clubb_kokkos_stats_netcdf PRIVATE cxx_std_17)
    target_link_libraries(clubb_kokkos_stats_netcdf PRIVATE clubb_kokkos_core_lib netcdf-fortran)
    add_test(NAME clubb_kokkos_stats_netcdf COMMAND clubb_kokkos_stats_netcdf)
endif()
