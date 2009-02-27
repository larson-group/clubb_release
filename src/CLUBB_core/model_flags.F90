!===============================================================================
! $Id$

module model_flags

! Description:
!   Various model options that can be toggled off and on as desired.

! References:
!   None
!-------------------------------------------------------------------------------

  implicit none

  public :: setup_model_flags
  private ! Default Scope

  logical, parameter, public ::  & 
    l_LH_on              = .false., & ! Latin hypercube calculation
    l_local_kk           = .false., & ! Local drizzle for rain microphysics
    l_hyper_dfsn         = .false., & ! 4th-order hyper-diffusion
    l_pos_def            = .false., & ! Flux limiting pos. def. scheme on rtm
    l_hole_fill          = .true.,  & ! Hole filling pos. def. scheme on wp2,up2,rtp2,etc
    l_clip_semi_implicit = .false., & ! Semi-implicit clipping scheme on wpthlp and wprtp
    l_3pt_sqd_dfsn       = .false., & ! Three-point squared diffusion coefficient
    l_clip_turb_adv      = .false.    ! Corrects thlm/rtm when w'th_l'/w'r_t' is clipped

  logical, parameter, public :: &
    l_standard_term_ta = .false.    ! Use the standard discretization for the
  ! turbulent advection terms.  Setting to
  ! .false. means that a_1 and a_3 are pulled
  ! outside of the derivative in advance_wp2_wp3_mod.F90
  ! and in advance_xp2_xpyp_module.F90.

  logical, parameter, public :: & 
    l_single_C2_Skw = .false.,  & ! Use a single Skw dependent value for C2
    l_gamma_Skw     = .true.,   & ! Use a Skw dependent gamma parameter
    l_byteswap_io   = .false.     ! Swap byte order in GrADS output

  logical, public :: & 
    l_bugsrad,       & ! BUGSrad interactive radiation scheme
    l_uv_nudge,      & ! For wind speed nudging. - Michael Falk
    l_soil_veg,      & ! Simple surface scheme - Joshua Fasching
    l_tke_aniso        ! For anisotropic turbulent kinetic energy,
                       !   i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

! OpenMP directives. These cannot be indented.
!$omp threadprivate(l_bugsrad, l_uv_nudge, l_tke_aniso)

  contains

!===============================================================================
  subroutine setup_model_flags & 
             ( l_bugsrad_in, l_soil_veg_in, & 
               l_uv_nudge_in, l_tke_aniso_in )

! Description:
!   Setup model flags

! References:
!   None
!-------------------------------------------------------------------------------
    use constants, only:  & 
      fstderr ! Variable(s)

    implicit none

        ! Input Variables
    logical, intent(in) ::  & 
      l_bugsrad_in, l_soil_veg_in, & 
      l_uv_nudge_in, & 
      l_tke_aniso_in

    !---- Begin Code ----
    
    l_soil_veg     = l_soil_veg_in
    l_bugsrad      = l_bugsrad_in
    l_uv_nudge     = l_uv_nudge_in
    l_tke_aniso    = l_tke_aniso_in


    return
  end subroutine setup_model_flags

end module model_flags
