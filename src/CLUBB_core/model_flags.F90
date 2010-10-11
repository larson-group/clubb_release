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
    l_clip_turb_adv      = .false., & ! Corrects thlm/rtm when w'th_l'/w'r_t' is clipped
    l_gmres              = .false., & ! Use GMRES iterative solver rather than LAPACK
    l_sat_mixrat_lookup  = .false.    ! Use a lookup table for mixing length
                                      ! saturation vapor pressure calculations

  logical, parameter, public :: &
    l_standard_term_ta = .false.    ! Use the standard discretization for the
  ! turbulent advection terms.  Setting to
  ! .false. means that a_1 and a_3 are pulled
  ! outside of the derivative in advance_wp2_wp3_module.F90
  ! and in advance_xp2_xpyp_module.F90.

  logical, parameter, public :: &
    l_single_C2_Skw = .false.,  & ! Use a single Skw dependent value for C2
#ifdef BYTESWAP_IO
    l_byteswap_io   = .true.,  & ! Don't use the native byte ordering in GrADS output
#else
    l_byteswap_io   = .false., & ! Use the native byte ordering in GrADS output
#endif
    l_gamma_Skw     = .true.      ! Use a Skw dependent gamma parameter

  logical, public :: & 
    l_uv_nudge,      & ! For wind speed nudging. - Michael Falk
    l_soil_veg,      & ! Simple surface scheme - Joshua Fasching
    l_tke_aniso        ! For anisotropic turbulent kinetic energy,
                       !   i.e. TKE = 1/2 (u'^2 + v'^2 + w'^2)

! OpenMP directives. These cannot be indented.
!$omp threadprivate(l_uv_nudge, l_tke_aniso, l_soil_veg)

  logical, parameter, public :: &
    l_use_boussinesq = .false.  ! Flag to use the Boussinesq form of the
                                ! predictive equations.  The predictive
                                ! equations are anelastic by default.

  ! Use to determine whether a host model has already applied the surface flux,
  ! to avoid double counting.
  logical, public :: &
    l_host_applies_sfc_fluxes 

!$omp threadprivate(l_host_applies_sfc_fluxes)

  integer, public :: &
    saturation_formula ! Integer that stores the saturation formula to be used

!$omp threadprivate(saturation_formula)

  ! These are the integer constants that represent the various saturation
  ! formulas. To add a new formula, add an additional constant here,
  ! add the logic to check the strings for the new formula in clubb_core and
  ! this module, and add logic in saturation to call the proper function--
  ! the control logic will be based on these named constants.

  integer, parameter, public :: &
    saturation_bolton = 1, & ! Constant for Bolton approximations of saturation
#ifdef GFDL
    saturation_gfdl   = 2, & ! Constant for the GFDL approximation of saturation
#endif
    saturation_flatau = 3    ! Constant for Flatau approximations of saturation

 
#ifdef GFDL
  logical, public :: &
     I_sat_sphum       ! h1g, 2010-06-15
#endif

  contains

!===============================================================================
  subroutine setup_model_flags & 
             ( l_soil_veg_in, l_host_applies_sfc_fluxes_in, & 
               l_uv_nudge_in, l_tke_aniso_in, saturation_formula_in &
#ifdef GFDL
               ,  I_sat_sphum_in   &  ! h1g, 2010-06-15
#endif
                )

! Description:
!   Setup model flags

! References:
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only:  & 
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

#ifdef GFDL
         logical, intent(in) ::  & 
         I_sat_sphum_in           ! h1g, 2010-06-15
#endif

    !---- Begin Code ----

    ! Logicals
    l_soil_veg     = l_soil_veg_in
    l_uv_nudge     = l_uv_nudge_in
    l_tke_aniso    = l_tke_aniso_in

    l_host_applies_sfc_fluxes = l_host_applies_sfc_fluxes_in

!    ! String
!    saturation_formula = trim( saturation_formula_in )

    ! Set up the saturation formula value
    select case ( trim( saturation_formula_in ) )
    case ( "bolton", "Bolton" )
      saturation_formula = saturation_bolton
    case ( "flatau", "Flatau" )
      saturation_formula = saturation_flatau
#ifdef GDFL
    case ( "gfdl", "GFDL" )
      saturation_formula = saturation_gfdl
#endif
    ! Add new saturation formulas after this.
    end select

#ifdef GFDL
      I_sat_sphum = I_sat_sphum_in  ! h1g, 2010-06-15
#endif

    return
  end subroutine setup_model_flags

end module model_flags
