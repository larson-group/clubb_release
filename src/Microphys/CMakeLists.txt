add_subdirectory(KK_microphys)          # provides:  KK_microphys_lib
add_subdirectory(COAMPS_microphys)      # provides:  coamps_lib
add_subdirectory(SCM_Activation)        # provides:  scm_activation_lib
add_subdirectory(Microphys_utils)       # provides:  microphys_utils_lib
add_subdirectory(Morrison_microphys)       # provides:  morrison_microphys_lib

set(MICROPHYS_SOURCES advance_microphys_module.F90
                      coamps_microphys_driver_module.F90  
                      microphys_driver.F90                
                      cloud_sed_module.F90
                      morrison_microphys_module.F90
                      KK_microphys_module.F90  
                      microphys_init_cleanup.F90
                      estimate_scm_microphys_module.F90
                      lh_microphys_driver_module.F90
                      parameters_microphys.F90
                      gfdl_activation.F90
                      ice_dfsn_module.F90
                      silhs_category_variance_module.F90
                      pdf_hydromet_microphys_wrapper.F90
                      mixed_moment_PDF_integrals.F90
                      microphys_interface.inc
                    )

add_library(microphys_lib ${MICROPHYS_SOURCES})

target_link_libraries(microphys_lib PRIVATE clubb_core_lib
                                            input_fields_lib
                                            KK_microphys_lib
                                            coamps_lib
                                            silhs_lib
                                            scm_activation_lib
                                            microphys_utils_lib
                                            morrison_microphys_lib
                                          )

if(${ENABLE_OMP})
  target_link_libraries(microphys_lib PRIVATE openmp-fortran)
endif()
