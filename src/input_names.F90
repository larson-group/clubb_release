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
  z_name = 'z[m]', &
  pressure_name = 'Press[Pa]', &
  press_mb_name = "Press[mb]", &
  temperature_name = 'T[K]', &
  theta_name = 'thm[K]', &
  thetal_name = 'thlm[K]', &
  rt_name = 'rt[kg\kg]', &
  sp_humidity_name = "sp_hmdty[kg\kg]",&
  um_name = 'u[m\s]', &
  vm_name = 'v[m\s]', &
  ug_name = 'ug[m\s]', &
  vg_name = 'vg[m\s]', &
  um_f_name = 'um_f[m\s^2]', &
  vm_f_name = 'vm_f[m\s^2]', &
  wm_name = 'w[m\s]', &
  omega_name = 'omega[Pa\s]', &
  CO2_name = 'CO2[ppmv]', &
  time_name = 'Time[s]', &
  LH_name = 'LH[W\m^2]', &
  SH_name = 'SH[W\m^2]', &
  thetal_f_name = 'thlm_f[K\s]', &
  rt_f_name = 'rtm_f[kg\kg\s]', &
  um_ref_name = 'um_ref[m\s]', &
  vm_ref_name = 'vm_ref[m\s]', &
  ozone_name = "o3[kg\kg]"
  
private ! Default Scope

end module input_names
