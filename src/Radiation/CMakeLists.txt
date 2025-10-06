add_subdirectory(BUGSrad)               # provides:  bugsrad_lib

set(RADIATION_SOURCES parameters_radiation.F90 
                      radiation_module.F90 
                      simple_rad_module.F90
                      bugsrad_driver.F90
                      cos_solar_zen_module.F90
                      rad_lwsw_module.F90
                      extended_atmosphere_module.F90
                    )

add_library(radiation_lib ${RADIATION_SOURCES})

target_link_libraries(radiation_lib PRIVATE clubb_core_lib 
                                            bugsrad_lib 
                                            silhs_lib
                                            text_writer_lib 
                                          )

if(${ENABLE_OMP})
  target_link_libraries(radiation_lib PRIVATE openmp-fortran)
endif()