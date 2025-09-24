# Establish libraries from this directory
set(CLUBB_DRIVER_SOURCES
        advance_microphys_module.F90
        clubb_driver.F90
        clubb_model_settings.F90
        bugsrad_driver.F90
        coamps_microphys_driver_module.F90
        cos_solar_zen_module.F90
        ice_dfsn_module.F90
        morrison_microphys_module.F90
        microphys_driver.F90
        microphys_init_cleanup.F90
        microphys_interface.inc
        mixed_moment_PDF_integrals.F90
        quicksort.F90
        simple_rad_module.F90
        soil_vegetation.F90
        sounding.F90
        stat_file_utils.F90
        time_dependent_input.F90
        extended_atmosphere_module.F90
        input_netcdf.F90
        input_interpret.F90
        input_fields.F90
        input_grads.F90
        gfdl_activation.F90
        extrapolation.F90
        hydrostatic_module.F90
        estimate_scm_microphys_module.F90
        lh_microphys_driver_module.F90
        parameters_radiation.F90
        parameters_microphys.F90
        KK_microphys_module.F90
        silhs_category_variance_module.F90
        enhanced_simann.F90
        error.F90
        pdf_hydromet_microphys_wrapper.F90
        generalized_grid_test.F90
        text_writer.F90
)

add_library(clubb_src_lib ${CLUBB_DRIVER_SOURCES})

if(${USE_NetCDF})
    target_link_libraries(clubb_src_lib PRIVATE netcdf-fortran)
endif()

if(${ENABLE_OMP})
    target_link_libraries(clubb_src_lib PRIVATE openmp-fortran)
endif()

# Add subdirectories we need
add_subdirectory(KK_microphys)          # provides:  KK_microphys_lib
add_subdirectory(COAMPS_microphys)      # provides:  coamps_lib
add_subdirectory(Morrison_microphys)    # provides:  morrison_microphys_lib
add_subdirectory(SCM_Activation)        # provides:  scm_activation_lib
add_subdirectory(SILHS)                 # provides:  silhs_lib
add_subdirectory(Numerical_recipes)     # provides:  numerical_recipes_lib
add_subdirectory(G_unit_test_types)     # provides:  G_unit_test_lib
add_subdirectory(Microphys_utils)       # provides:  microphys_utils_lib
add_subdirectory(BUGSrad)               # provides:  bugsrad_lib
add_subdirectory(CLUBB_core)            # provides:  clubb_core_lib
add_subdirectory(Lapack)                # provides:  lapack_lib
add_subdirectory(Benchmark_cases)       # provides:  additional sources added to clubb_src_lib

# add source file to BUGGSrad lib
target_sources(bugsrad_lib PRIVATE variables_radiation_module.F90)

if(${ENABLE_OMP})
    target_link_libraries(bugsrad_lib PRIVATE openmp-fortran)
endif()

# Create libclubb_parabolic and link in clubb_core_lib 
add_library(parabolic_lib AiryFunction.f90 Parabolic_constants.f90 Parabolic.f90)
target_link_libraries(parabolic_lib PRIVATE clubb_core_lib)

# Define link interface for clubb_src_lib
target_link_libraries(clubb_src_lib PRIVATE clubb_core_lib
                                            bugsrad_lib
                                            KK_microphys_lib
                                            coamps_lib
                                            morrison_microphys_lib
                                            scm_activation_lib
                                            silhs_lib
                                        )

# Create clubb_standalone executable
add_executable(clubb_standalone clubb_standalone.F90)
target_link_libraries(clubb_standalone PRIVATE  bugsrad_lib
                                                clubb_core_lib
                                                KK_microphys_lib
                                                coamps_lib
                                                morrison_microphys_lib
                                                clubb_src_lib
                                                scm_activation_lib
                                                silhs_lib
                                            )

# Create clubb_thread_test exectuable
add_executable(clubb_thread_test clubb_thread_test.F90)
target_link_libraries(clubb_thread_test PRIVATE bugsrad_lib
                                                clubb_core_lib
                                                KK_microphys_lib
                                                coamps_lib
                                                morrison_microphys_lib
                                                clubb_src_lib
                                                scm_activation_lib
                                                silhs_lib
                                            )
if(${ENABLE_OMP})
    target_link_libraries(clubb_thread_test PRIVATE openmp-fortran)
endif()

# Create clibb_tuner executable
add_executable(clubb_tuner clubb_tuner.F90)
target_link_libraries(clubb_tuner PRIVATE   bugsrad_lib
                                            clubb_core_lib
                                            KK_microphys_lib
                                            coamps_lib
                                            morrison_microphys_lib
                                            clubb_src_lib
                                            scm_activation_lib
                                            silhs_lib
                                            numerical_recipes_lib
                                        )

# Create jacobian executable
add_executable(jacobian jacobian.F90)
target_link_libraries(jacobian PRIVATE  bugsrad_lib
                                        clubb_core_lib
                                        KK_microphys_lib
                                        coamps_lib
                                        morrison_microphys_lib
                                        clubb_src_lib
                                        scm_activation_lib
                                        silhs_lib
                                    )

# Create G_unit_tests executable
add_executable(G_unit_tests G_unit_tests.F90)
target_link_libraries(G_unit_tests PRIVATE  bugsrad_lib
                                            clubb_core_lib
                                            KK_microphys_lib
                                            coamps_lib
                                            morrison_microphys_lib
                                            clubb_src_lib
                                            scm_activation_lib
                                            silhs_lib
                                            G_unit_test_lib
                                        )

# Create int2text executable
add_executable(int2text int2txt.F90)
target_link_libraries(int2text PRIVATE  bugsrad_lib
                                        clubb_core_lib
                                        KK_microphys_lib
                                        coamps_lib
                                        morrison_microphys_lib
                                        clubb_src_lib
                                    )

install(TARGETS clubb_core_lib
                parabolic_lib
                microphys_utils_lib
                KK_microphys_lib
                silhs_lib
                morrison_microphys_lib
                scm_activation_lib
                coamps_lib
                bugsrad_lib
                clubb_src_lib
                clubb_tuner
                G_unit_tests
                clubb_thread_test
                jacobian
                clubb_standalone
                int2text
        DESTINATION ${CMAKE_INSTALL_PREFIX})
