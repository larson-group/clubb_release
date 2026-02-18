# Add subdirectories we need
add_subdirectory(Radiation)             # provides:  radiation_lib
add_subdirectory(Microphys)             # provides:  microphys_lib
add_subdirectory(SILHS)                 # provides:  silhs_lib
add_subdirectory(Numerical_recipes)     # provides:  numerical_recipes_lib
add_subdirectory(G_unit_test_types)     # provides:  G_unit_test_lib
add_subdirectory(CLUBB_core)            # provides:  clubb_core_lib
add_subdirectory(Lapack)                # provides:  lapack_lib
add_subdirectory(Input_fields)          # provides:  input_fields_lib
add_subdirectory(Benchmark_cases)       # provides:  benchmark_cases_lib

# create text_writer_lib used in various files
add_library(text_writer_lib text_writer.F90)
target_link_libraries(text_writer_lib PRIVATE clubb_core_lib)

# Establish libraries from this directory
set(CLUBB_DRIVER_SOURCES
        clubb_driver.F90
        soil_vegetation.F90
        generalized_grid_test.F90  
)

add_library(clubb_driver_lib ${CLUBB_DRIVER_SOURCES})

# Define link interface for clubb_driver_lib
target_link_libraries(clubb_driver_lib PRIVATE  clubb_core_lib
                                                input_fields_lib
                                                radiation_lib
                                                silhs_lib
                                                benchmark_cases_lib
                                                microphys_lib
                                                text_writer_lib
                                            )

if(${USE_NetCDF})
    target_link_libraries(clubb_driver_lib PUBLIC netcdf-fortran)
endif()

if(${ENABLE_OMP})
    target_link_libraries(clubb_driver_lib PRIVATE openmp-fortran)
endif()

if(${USE_GPTL})
    target_link_libraries(clubb_driver_lib PRIVATE gptl)
endif()

# Create clubb_standalone executable
add_executable(clubb_standalone clubb_standalone.F90)
target_link_libraries(clubb_standalone PRIVATE clubb_driver_lib)

# Create clubb_thread_test exectuable
add_executable(clubb_thread_test clubb_thread_test.F90)
target_link_libraries(clubb_thread_test PRIVATE clubb_driver_lib)

if(${ENABLE_OMP})
    target_link_libraries(clubb_thread_test PRIVATE openmp-fortran)
endif()

# Create clubb_tuner executable
add_executable(clubb_tuner clubb_tuner.F90 error.F90)
target_link_libraries(clubb_tuner PRIVATE clubb_driver_lib 
                                          input_fields_lib
                                          numerical_recipes_lib
                                        )

# Create G_unit_tests executable
add_executable(G_unit_tests G_unit_tests.F90)
target_link_libraries(G_unit_tests PRIVATE G_unit_test_lib )

# Create clubb_driver_test executable
add_executable(clubb_driver_test clubb_driver_test.F90)
target_link_libraries(clubb_driver_test PRIVATE clubb_driver_lib )

# Commented out, since this appears unused as of 2025
# # Create jacobian executable
# add_executable(jacobian jacobian.F90)
# target_link_libraries(jacobian PRIVATE clubb_driver_lib)


# Commented out, since this appears unused as of 2025
# Create int2txt executable
# add_executable(int2txt int2txt.F90)
# target_link_libraries(int2txt PRIVATE clubb_driver_lib)

install(TARGETS clubb_driver_lib
                clubb_tuner
                G_unit_tests
                clubb_thread_test
                clubb_standalone
                clubb_driver_test
                #jacobian
                #int2txt
        DESTINATION ${CMAKE_INSTALL_PREFIX})
