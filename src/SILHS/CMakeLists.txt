set(SILHS_SOURCES
                    est_kessler_microphys_module.F90
                    generate_uniform_sample_module.F90
                    latin_hypercube_arrays.F90
                    latin_hypercube_driver_module.F90
                    lh_microphys_var_covar_module.F90
                    math_utilities.F90
                    output_2D_samples_module.F90
                    parameters_silhs.F90
                    silhs_api_module.F90
                    silhs_importance_sample_module.F90
                    transform_to_pdf_module.F90
)

add_library(silhs_lib ${SILHS_SOURCES})
target_link_libraries(silhs_lib PRIVATE clubb_core_lib KK_microphys_lib microphys_utils_lib)
