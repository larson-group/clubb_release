set(BUGSrad_SOURCES
        bugs_lwr.F
        bugs_rad.F
        bugsrad_physconst.F90
        bugsrad_planck.F90
        bugs_swr.F
        cloud_correlate.F90
        cloudg.F
        comscp1.F
        comscp2.F
        gascon.F
        gases_ckd_data.h
        gases_ckd.F90
        kinds.F90
        newexp.F90
        rayle.F
        two_rt_lw.F
        two_rt_lw_iter.F
        two_rt_lw_ocastrndm.F90
        two_rt_lw_sel.F
        two_rt_sw_bs.F
        two_rt_sw.F
        two_rt_sw_iter.F
        two_rt_sw_ocastrndm.F90
        variables_radiation_module.F90
)

add_library(bugsrad_lib ${BUGSrad_SOURCES})

target_link_libraries(bugsrad_lib PRIVATE clubb_core_lib lapack_lib)

if(${ENABLE_OMP})
    target_link_libraries(bugsrad_lib PRIVATE openmp-fortran)
endif()
