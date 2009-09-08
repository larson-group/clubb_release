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
    l_uv_nudge,      & ! For wind speed nudging. - Michael Falk
    l_soil_veg,      & ! Simple surface scheme - Joshua Fasching
    l_tke_aniso        ! For anisotropic turbulent kinetic energy,
                       !   i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

  ! Use to determine whether a host model has already applied the surface flux,
  ! to avoid double counting.
  logical, public :: &
    l_host_applies_sfc_fluxes 

  character(len=6), public :: &
    saturation_formula ! "bolton" approx. or "flatau" approx.

! OpenMP directives. These cannot be indented.
!$omp threadprivate(l_uv_nudge, l_tke_aniso, l_host_applies_sfc_fluxes, &
!$omp   saturation_formula)

  contains

!===============================================================================
  subroutine setup_model_flags & 
             ( l_soil_veg_in, l_host_applies_sfc_fluxes_in, & 
               l_uv_nudge_in, l_tke_aniso_in, saturation_formula_in )

! Description:
!   Setup model flags

! References:
!   None
!-------------------------------------------------------------------------------
    use constants, only:  & 
      fstderr ! Variable(s)

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    logical, intent(in) ::  & 
      l_soil_veg_in, & 
      l_host_applies_sfc_fluxes_in, &
      l_uv_nudge_in, & 
      l_tke_aniso_in

    character(len=*), intent(in) :: &
      saturation_formula_in

    !---- Begin Code ----

    ! Logicals
    l_soil_veg     = l_soil_veg_in
    l_uv_nudge     = l_uv_nudge_in
    l_tke_aniso    = l_tke_aniso_in

    l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes_in

    ! String
    saturation_formula = trim( saturation_formula_in )

    return
  end subroutine setup_model_flags

end module model_flags
