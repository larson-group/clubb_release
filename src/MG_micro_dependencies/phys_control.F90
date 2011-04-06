! $Id$
module phys_control
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  implicit none

  private

  public :: phys_getopts
  
  contains
  
!================================================================================================
  subroutine phys_getopts &
             (deep_scheme_out, shallow_scheme_out, eddy_scheme_out, microp_scheme_out, &
             macrop_scheme_out, atm_dep_flux_out, history_aerosol_out, history_microphysics_out, &
             history_budget_out, history_budget_histfile_num_out, do_tms_out, do_iss_out, &
             tms_orocnst_out, tms_z0fac_out, conv_water_in_rad_out )
    !
    !  Description: The original subroutine sets options for determining
    !               turbulent diffusion/microphysics history output
    !               In our case, we don't use this, so this subroutine does nothing.
    !
    !---------------------------------------------------------------------------------
  
       character(len=16), intent(out), optional :: deep_scheme_out
       character(len=16), intent(out), optional :: shallow_scheme_out
       character(len=16), intent(out), optional :: eddy_scheme_out
       character(len=16), intent(out), optional :: microp_scheme_out
       character(len=16), intent(out), optional :: macrop_scheme_out
       logical,           intent(out), optional :: atm_dep_flux_out
       logical,           intent(out), optional :: history_aerosol_out
       logical,           intent(out), optional :: history_microphysics_out
       logical,           intent(out), optional :: history_budget_out
       integer,           intent(out), optional :: history_budget_histfile_num_out
       logical,           intent(out), optional :: do_tms_out
       logical,           intent(out), optional :: do_iss_out
       real,              intent(out), optional :: tms_orocnst_out
       real,              intent(out), optional :: tms_z0fac_out
       integer,           intent(out), optional :: conv_water_in_rad_out
       
       ! Avoid compiler warnings
       deep_scheme_out = ''
       shallow_scheme_out = ''
       eddy_scheme_out = ''
       microp_scheme_out = ''
       macrop_scheme_out = ''
       atm_dep_flux_out = .false.
       history_aerosol_out = .false.
       history_microphysics_out = .false.
       history_budget_out = .false.
       history_budget_histfile_num_out = 0
       do_tms_out = .false.
       do_iss_out = .false.
       tms_orocnst_out = 0
       tms_z0fac_out = 0
       conv_water_in_rad_out = 0
      
  end subroutine phys_getopts
  
!================================================================================================
    

end module phys_control
