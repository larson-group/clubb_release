!$Id$
module input_names
!
!  Description: This module contains all of the strings used to define the
!  headers for input_reader.F90 compatable files.
!
!---------------------------------------------------------------------------------------------------
  implicit none
  ! Column identifiers
  character(len=*), public, parameter :: &
    z_name = 'z[m]'

  character(len=*), public, parameter :: &
    pressure_name = 'Press[Pa]', &
    press_mb_name = "Press[mb]"

  character(len=*), public, parameter :: &
    temperature_name = 'T[K]', &
    theta_name = 'thm[K]', &
    thetal_name = 'thlm[K]'

  character(len=*), public, parameter :: &
    temperature_f_name = 'T_f[K\s]', &
    thetal_f_name = 'thlm_f[K\s]', &
    theta_f_name = 'thm_f[K\s]'

  character(len=*), public, parameter :: &
    rt_name = 'rt[kg\kg]', &
    sp_humidity_name = "sp_hmdty[kg\kg]"

  character(len=*), public, parameter :: &
    rt_f_name = 'rtm_f[kg\kg\s]', &
    sp_humidity_f_name = 'sp_hmdty_f[kg\kg\s]'

  character(len=*), public, parameter :: &
    um_name = 'u[m\s]', &
    vm_name = 'v[m\s]'

  character(len=*), public, parameter :: &
    ug_name = 'ug[m\s]', &
    vg_name = 'vg[m\s]'

  character(len=*), public, parameter :: &
    um_ref_name = 'um_ref[m\s]', &
    vm_ref_name = 'vm_ref[m\s]'

  character(len=*), public, parameter :: &
    um_f_name = 'um_f[m\s^2]', &
    vm_f_name = 'vm_f[m\s^2]'

  character(len=*), public, parameter :: &
    wm_name = 'w[m\s]', &
    omega_name = 'omega[Pa\s]', &
    omega_mb_hr_name = 'omega[mb\hr]'

  character(len=*), public, parameter :: &
    CO2_name = 'CO2[ppmv]', &
    ozone_name = "o3[kg\kg]"

  character(len=*), public, parameter :: &
    time_name = 'Time[s]'

  character(len=*), public, parameter :: &
    LH_name = 'LH[W\m^2]', &
    SH_name = 'SH[W\m^2]'

  private ! Default Scope

end module input_names
