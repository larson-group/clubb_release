!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module sponge_layer_damping

  ! Description:
  ! This module is used for damping variables in upper altitudes of the grid.
  !
  ! References:
  !   None
  !-------------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  public :: sponge_damp_xm,             & ! Procedure(s)
            initialize_tau_sponge_damp, &
            finalize_tau_sponge_damp,   &
            sponge_damp_settings,       & ! Variable type(s)
            sponge_damp_profile


  type sponge_damp_settings

    real( kind = core_rknd ) :: &
      tau_sponge_damp_min, & ! Minimum damping time scale (model top)        [s]
      tau_sponge_damp_max, & ! Maximum damping time scale (damp layer base)  [s]
      sponge_damp_depth      ! Damping depth as a fraction of domain height  [-]

    logical :: &
      l_sponge_damping       ! True if damping is being used

  end type sponge_damp_settings

  type sponge_damp_profile

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      tau_sponge_damp    ! Damping time scale    [1/s]

    real( kind = core_rknd ) :: &
      sponge_layer_depth    ! Depth of sponge damping layer  [m]

  end type sponge_damp_profile


  type(sponge_damp_settings), public :: &
    thlm_sponge_damp_settings, & ! Variable(s)
    rtm_sponge_damp_settings,  &
    uv_sponge_damp_settings
!$omp threadprivate( thlm_sponge_damp_settings, rtm_sponge_damp_settings, &
!$omp                uv_sponge_damp_settings )

  type(sponge_damp_profile), public :: &
    thlm_sponge_damp_profile, & ! Variable(s)
    rtm_sponge_damp_profile,  &
    uv_sponge_damp_profile
!$omp threadprivate( thlm_sponge_damp_profile, rtm_sponge_damp_profile, &
!$omp                uv_sponge_damp_profile )


  private

  contains

  !=============================================================================
  function sponge_damp_xm( dt, z, xm_ref, xm, damping_profile ) result( xm_p )

    ! Description:
    ! Damps specified mean field toward a reference profile.  The module must be
    ! initialized for this function to work.  Otherwise a stop is issued.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    !  "Sponge"-layer damping at the domain top region

    use grid_class, only: &
        gr    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! External
    intrinsic :: allocated

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z,      & ! Height of model grid levels                [m]
      xm_ref, & ! Reference profile of x to damp xm towards  [units vary]
      xm        ! Mean field being damped                    [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xm_p   ! Damped value of xm  [units_vary]

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of timestep to damping timescale  [-]

    integer :: k

    ! ---- Begin Code ----

    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       xm_p = xm
     
       do k = gr%nz, 1, -1

          ! The height of the model top is gr%zm(gr%nz).
          if ( gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth ) then

             ! Vince Larson used implicit discretization in order to 
             ! reduce noise in rtm in cloud_feedback_s12 (CGILS) 
             !xm_p(k) = xm(k) - real( ( ( xm(k) - xm_ref(k) ) / & 
             !                  damping_profile%tau_sponge_damp(k) ) * dt )
             dt_on_tau = dt / damping_profile%tau_sponge_damp(k)

             ! Really, we should be using xm_ref at time n+1 rather than n.
             ! However, for steady profiles of xm_ref, it won't matter.        
             xm_p(k) = ( xm(k) + dt_on_tau * xm_ref(k) ) / &
                             ( 1.0_core_rknd + dt_on_tau )
             ! End Vince Larson's change

          else ! gr%zm(gr%nz) - z(k) >= damping_profile%sponge_layer_depth

             ! Below sponge damping layer; exit loop.
             exit

          endif ! gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth


       enddo ! k = gr%nz, 1, -1

    else

       stop "tau_sponge_damp in damping used before initialization"

    endif


    return

  end function sponge_damp_xm

  !=============================================================================
  subroutine initialize_tau_sponge_damp( dt, z, settings, damping_profile )

    ! Description:
    ! Initialize time scale, tau_sponge_damp, used for damping.  The time scale
    ! attains its maximum value, tau_sponge_damp_max, at the bottom of the
    ! "sponge" damping layer, which results in minimal damping.  Likewise, the
    ! time scale attains its minimum value, tau_sponge_damp_min, at the top of
    ! the model, which results in maximum damping.  At levels in-between the top
    ! of the model and the base of the sponge damping layer, the value of
    ! tau_sponge_damp is in-between tau_sponge_damp_min and tau_sponge_damp_max,
    ! as calculated by an interpolation formula.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd    ! Variable(s)
    
    use constants_clubb, only: &
        two,     & ! Constant(s)
        fstderr

    use grid_class, only: &
        gr    ! Variable(s)

!    use interpolation, only: &
!        lin_interpolate_two_points    ! Procedure(s)

    implicit none

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep    [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z    ! Height of model grid levels    [m]

    type(sponge_damp_settings), intent(in) :: &
      settings

    ! Output Variable(s)
    type(sponge_damp_profile), intent(out) :: &
      damping_profile

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      tau_sponge_damp_exponent  ! Exponent in calculation of tau_sponge_damp [-]

    integer :: &
      k    ! Loop index

    ! ---- Begin Code ----

    ! Allocate the damping time scale.
    allocate( damping_profile%tau_sponge_damp(1:gr%nz) )

    ! Calculate the depth of the sponge layer.
    ! The height of the model top is gr%zm(gr%nz).
    damping_profile%sponge_layer_depth &
    = settings%sponge_damp_depth * gr%zm(gr%nz)

    ! Check the value of tau_sponge_damp_min.
    if ( settings%tau_sponge_damp_min < two * dt ) then
       write(fstderr,*) "Error:  tau_sponge_damp_min is too small!"
       write(fstderr,*) "It must be at least 2.0 * dt"
       stop
    endif

    ! Calculate the value of the damping time scale, tau_sponge_damp, at levels
    ! that are within the sponge damping layer.
    do k = gr%nz, 1, -1

       ! The height of the model top is gr%zm(gr%nz).
       if ( gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth ) then

          ! Vince Larson added code to use standard linear interpolation.
          ! Brian Griffin reverted the linear interpolation in order to use code
          ! that is similar to what is found in SAM LES.

          tau_sponge_damp_exponent &
          = ( gr%zm(gr%nz) - z(k) ) / damping_profile%sponge_layer_depth

          damping_profile%tau_sponge_damp(k) &
          = settings%tau_sponge_damp_min &
            * ( settings%tau_sponge_damp_max &
                / settings%tau_sponge_damp_min )**tau_sponge_damp_exponent

          !damping_profile%tau_sponge_damp(k) &
          != lin_interpolate_two_points( z(k), gr%zm(gr%nz), &
          !                              gr%zm(gr%nz) &
          !                              - damping_profile%sponge_layer_depth, &
          !                              settings%tau_sponge_damp_min, &
          !                              settings%tau_sponge_damp_max )

          ! End Vince Larson's change
          ! End Brian Griffin's rebellious reversion.

       else ! gr%zm(gr%nz) - z(k) >= damping_profile%sponge_layer_depth

          ! Below sponge damping layer; exit loop.
          exit

       endif ! gr%zm(gr%nz) - z(k) < damping_profile%sponge_layer_depth

    enddo ! k = gr%nz, 1, -1


    return

  end subroutine initialize_tau_sponge_damp

  !=============================================================================
  subroutine finalize_tau_sponge_damp( damping_profile )

    ! Description:
    ! Frees memory allocated in initialize_tau_sponge_damp
    ! 
    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input/Output Variable(s)
    type(sponge_damp_profile), intent(inout) :: &
      damping_profile ! Information for damping the profile

    ! ---- Begin Code ----

    ! Deallocate the damping time scale.
    deallocate( damping_profile%tau_sponge_damp )


    return

  end subroutine finalize_tau_sponge_damp

!===============================================================================

end module sponge_layer_damping
