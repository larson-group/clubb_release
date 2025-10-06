set(KK_MICROPHYS_SOURCES
        KK_utilities.F90
        KK_Nrm_tendencies.F90
        KK_local_means.F90
        KK_upscaled_means.F90
        KK_upscaled_turbulent_sed.F90
        KK_upscaled_covariances.F90
        KK_upscaled_variances.F90
        PDF_integrals_means.F90
        PDF_integrals_covar.F90
        PDF_integrals_all_MM.F90
        parameters_KK.F90
        Parabolic.f90
        Parabolic_constants.f90
        AiryFunction.f90
)

add_library(KK_microphys_lib ${KK_MICROPHYS_SOURCES})

target_link_libraries(KK_microphys_lib PRIVATE clubb_core_lib microphys_utils_lib)