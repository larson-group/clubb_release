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

  public :: sponge_damp_xm,                 & ! Procedure(s)
            sponge_damp_xp2,                &
            sponge_damp_xp3,                &
            initialize_tau_sponge_damp_api, &
            finalize_tau_sponge_damp_api,   &
            sponge_damp_settings,           & ! Variable type(s)
            sponge_damp_profile


  type sponge_damp_settings

    real( kind = core_rknd ) :: &
      tau_sponge_damp_min, & ! Minimum damping time scale (model top)        [s]
      tau_sponge_damp_max, & ! Maximum damping time scale (damp layer base)  [s]
      sponge_damp_depth      ! Damping depth as a fraction of domain height  [-]

    logical :: &
      l_sponge_damping = .false.       ! True if damping is being used, default false

  end type sponge_damp_settings

  type sponge_damp_profile

    real( kind = core_rknd ), allocatable, dimension(:,:) :: &
      tau_sponge_damp    ! Damping time scale by column and level [s]

    real( kind = core_rknd ), allocatable, dimension(:) :: &
      sponge_layer_depth    ! Depth of sponge damping layer  [m]

  end type sponge_damp_profile


  type(sponge_damp_settings), public :: &
    thlm_sponge_damp_settings,    & ! Variable(s)
    rtm_sponge_damp_settings,     &
    uv_sponge_damp_settings,      &
    wp2_sponge_damp_settings,     &
    wp3_sponge_damp_settings,     &
    up2_vp2_sponge_damp_settings

  type(sponge_damp_profile), public :: &
    thlm_sponge_damp_profile,    & ! Variable(s)
    rtm_sponge_damp_profile,     &
    uv_sponge_damp_profile,      &
    wp2_sponge_damp_profile,     &
    wp3_sponge_damp_profile,     &
    up2_vp2_sponge_damp_profile


  private

  !$omp threadprivate( thlm_sponge_damp_settings, rtm_sponge_damp_settings )
  !$omp threadprivate( uv_sponge_damp_settings, wp2_sponge_damp_settings )
  !$omp threadprivate( wp3_sponge_damp_settings, up2_vp2_sponge_damp_settings )
  !$omp threadprivate( thlm_sponge_damp_profile, rtm_sponge_damp_profile )
  !$omp threadprivate( uv_sponge_damp_profile, wp2_sponge_damp_profile )
  !$omp threadprivate( wp3_sponge_damp_profile, up2_vp2_sponge_damp_profile )

  contains

  !=============================================================================
  function sponge_damp_xm( ngrdcol, nzm, nzt, gr, dt, zt, zm, &
                           xm_ref, xm, damping_profile ) result( xm_p )

    ! Description:
    ! Damps each column of a mean field toward its reference profile.
    ! The damping profile must be initialized before calling this function;
    ! otherwise a stop is issued.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    !  "Sponge"-layer damping at the domain top region
    use grid_class, only: &
            grid    ! Type(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! External
    intrinsic :: allocated

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol, & ! Number of columns [-]
      nzm,     & ! Number of momentum levels [-]
      nzt        ! Number of thermodynamic levels [-]

    type( grid ), intent(in) :: &
      gr    ! Grid bounds and vertical direction

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      zt    ! Thermodynamic grid heights [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      zm    ! Momentum grid heights [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      xm_ref, & ! Reference mean field [units vary]
      xm        ! Mean field before damping [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile    ! Column time scales and layer depths

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      xm_p   ! Damped value of xm  [units_vary]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of timestep to damping timescale  [-]

    integer :: i, k    ! Column and vertical loop indices [-]


    !--------------------------- Begin Code ---------------------------

    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       xm_p = xm
     
       do k = gr%k_ub_zt, gr%k_lb_zt, -gr%grid_dir_indx
          do i = 1, ngrdcol

             ! The height of the model top is zm(i,gr%k_ub_zm).
             if ( zm(i,gr%k_ub_zm) - zt(i,k) < damping_profile%sponge_layer_depth(i) ) then

                ! Vince Larson used implicit discretization in order to
                ! reduce noise in rtm in cloud_feedback_s12 (CGILS)
                dt_on_tau = dt / damping_profile%tau_sponge_damp(i,k)

                ! Really, we should be using xm_ref at time n+1 rather than n.
                ! However, for steady profiles of xm_ref, it won't matter.
                xm_p(i,k) = ( xm(i,k) + dt_on_tau * xm_ref(i,k) ) / &
                                ( 1.0_core_rknd + dt_on_tau )
                ! End Vince Larson's change

             endif ! zm(i,gr%k_ub_zm) - zt(i,k) < damping_profile%sponge_layer_depth(i)

          end do ! i
       end do ! k

    else

       error stop "tau_sponge_damp in sponge_damp_xm used before initialization"

    endif


    return

  end function sponge_damp_xm

  !=============================================================================
  function sponge_damp_xp2( ngrdcol, nzm, gr, dt, zm, xp2, x_tol_sqd, &
                            damping_profile ) &
  result( xp2_damped )

    ! Description:
    ! Calculates "sponge"-layer damping of the variance xp2 in each column.
    !
    ! Sponge damping on a local value of x is given by the equation:
    !
    ! x_d = x - ( delta_t / tau ) * ( x - <x> ),
    !
    ! where x is the local value prior to damping, x_d is the damped local value
    ! of x, <x> is the grid-level mean value of x, delta_t is the model time
    ! step duration, and tau is the damping time scale.  Since delta_t / tau has
    ! the same value everywhere at a grid level, the grid-level mean of x does
    ! not change as a result of damping.
    !
    ! Subtracting <x> from both sides:
    !
    ! x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),
    !
    ! which results in:
    !
    ! x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).
    !
    ! Squaring both sides:
    !
    ! ( x_d - <x> )^2 = ( 1 - delta_t / tau )^2 * ( x - <x> )^2.
    !
    ! After averaging both sides, the damped value of xp2 is:
    !
    ! < x_d'^2 > = ( 1 - delta_t / tau )^2 * < x'^2 >.
    !
    ! Any sponge damping is applied to (predictive) xp2 after the field has been
    ! advanced in time.  This allows sponge damping to be applied in an implicit
    ! manner.  The damped value of xp2 must also be limited at a minimum value
    ! of x_tol^2.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid

    use constants_clubb, only: &
        one    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol, & ! Number of columns [-]
      nzm        ! Number of momentum levels [-]

    type( grid ), intent(in) :: &
      gr    ! Grid bounds and vertical direction

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      zm,  & ! Momentum grid heights [m]
      xp2    ! Variance before damping [units vary]

    real( kind = core_rknd ), intent(in) :: &
      x_tol_sqd    ! Minimum variance [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile    ! Column time scales and layer depths

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      xp2_damped    ! Variance of x, <x'^2>, after damping   [units vary]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of model time step duration to damping timescale  [-]

    integer :: &
      i, k    ! Column and vertical loop indices [-]


    !--------------------------- Begin Code ---------------------------

    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       ! Set the entire profile of <x'^2> after damping to the profile of <x'^2>
       ! before damping.  The values of <x'^2> after damping will be overwritten
       ! at any levels where "sponge"-layer damping occurs.
       xp2_damped = xp2
     
       do k = gr%k_ub_zm, gr%k_lb_zm, -gr%grid_dir_indx
          do i = 1, ngrdcol

             ! The height of the model top is zm(i,gr%k_ub_zm).
             if ( zm(i,gr%k_ub_zm) - zm(i,k) < damping_profile%sponge_layer_depth(i) ) then

                ! Calculate the value of delta_t / tau at the grid level.
                dt_on_tau = dt / damping_profile%tau_sponge_damp(i,k)

                ! Calculate the damped value of <x'^2>.
                xp2_damped(i,k) = ( one - dt_on_tau )**2 * xp2(i,k)

                ! The damped value of <x'^2> needs to be greater than or equal to
                ! x_tol^2.
                xp2_damped(i,k) = max( xp2_damped(i,k), x_tol_sqd )

             endif ! zm(i,gr%k_ub_zm) - zm(i,k) < damping_profile%sponge_layer_depth(i)

          end do ! i
       end do ! k

    else

       error stop "tau_sponge_damp in sponge_damp_xp2 used before initialization"

    endif


    return

  end function sponge_damp_xp2

  !=============================================================================
  function sponge_damp_xp3( ngrdcol, nzm, nzt, gr, dt, z, zm, xp3, &
                            damping_profile ) &
  result( xp3_damped )

    ! Description:
    ! Calculates "sponge"-layer damping of the third moment xp3 in each column.
    !
    ! Sponge damping on a local value of x is given by the equation:
    !
    ! x_d = x - ( delta_t / tau ) * ( x - <x> ),
    !
    ! where x is the local value prior to damping, x_d is the damped local value
    ! of x, <x> is the grid-level mean value of x, delta_t is the model time
    ! step duration, and tau is the damping time scale.  Since delta_t / tau has
    ! the same value everywhere at a grid level, the grid-level mean of x does
    ! not change as a result of damping.
    !
    ! Subtracting <x> from both sides:
    !
    ! x_d - <x> = ( x - <x> ) - ( delta_t / tau ) * ( x - <x> ),
    !
    ! which results in:
    !
    ! x_d - <x> = ( 1 - delta_t / tau ) * ( x - <x> ).
    !
    ! Taking both sides to the third power:
    !
    ! ( x_d - <x> )^3 = ( 1 - delta_t / tau )^3 * ( x - <x> )^3.
    !
    ! After averaging both sides, the damped value of xp3 is:
    !
    ! < x_d'^3 > = ( 1 - delta_t / tau )^3 * < x'^3 >.
    !
    ! Any sponge damping is applied to (predictive) xp3 after the field has been
    ! advanced in time.  This allows sponge damping to be applied in an implicit
    ! manner.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid  ! Type(s)

    use constants_clubb, only: &
        one    ! Constant(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol, & ! Number of columns [-]
      nzm,     & ! Number of momentum levels [-]
      nzt        ! Number of thermodynamic levels [-]

    type( grid ), intent(in) :: &
      gr    ! Grid bounds and vertical direction

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep  [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      z    ! Thermodynamic grid heights [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      zm    ! Momentum grid heights [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      xp3    ! Third moment before damping [units vary]

    type(sponge_damp_profile), intent(in) :: &
      damping_profile    ! Column time scales and layer depths

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      xp3_damped    ! <x'^3> after damping   [units vary]

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ) :: &
      dt_on_tau    ! Ratio of model time step duration to damping timescale  [-]

    integer :: &
      i, k    ! Column and vertical loop indices [-]


    !--------------------------- Begin Code ---------------------------

    if ( allocated( damping_profile%tau_sponge_damp ) ) then

       ! Set the entire profile of <x'^3> after damping to the profile of <x'^3>
       ! before damping.  The values of <x'^3> after damping will be overwritten
       ! at any levels where "sponge"-layer damping occurs.
       xp3_damped = xp3
     
       do k = gr%k_ub_zt, gr%k_lb_zt, -gr%grid_dir_indx
          do i = 1, ngrdcol

             ! The height of the model top is zm(i,gr%k_ub_zm).
             if ( zm(i,gr%k_ub_zm) - z(i,k) < damping_profile%sponge_layer_depth(i) ) then

                ! Calculate the value of delta_t / tau at the grid level.
                dt_on_tau = dt / damping_profile%tau_sponge_damp(i,k)

                ! Calculate the damped value of <x'^3>.
                xp3_damped(i,k) = ( one - dt_on_tau )**3 * xp3(i,k)

             endif ! zm(i,gr%k_ub_zm) - z(i,k) < damping_profile%sponge_layer_depth(i)

          end do ! i
       end do ! k

    else

       error stop "tau_sponge_damp in sponge_damp_xp3 used before initialization"

    endif


    return

  end function sponge_damp_xp3

  !=============================================================================
  subroutine initialize_tau_sponge_damp_api( ngrdcol, nz, gr, dt, z, & ! In
                                             settings, & ! In
                                             damping_profile ) ! Out

    ! Description:
    ! Initialize the damping time scale independently in each column.  The time scale
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
        grid ! Type

!    use interpolation, only: &
!        lin_interpolate_two_points    ! Procedure(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol, & ! Number of columns [-]
      nz         ! Number of vertical levels [-]

    type (grid), intent(in) :: gr    ! Model grid

    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep    [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      z     ! Height of model grid levels    [m]

    type(sponge_damp_settings), intent(in) :: &
      settings    ! Shared damping configuration

    !--------------------------- Output Variables ---------------------------
    type(sponge_damp_profile), intent(out) :: &
      damping_profile    ! Column time scales and layer depths

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ) :: &
      tau_sponge_damp_exponent  ! Exponent in calculation of tau_sponge_damp [-]

    integer :: &
      i,           & ! Column index [-]
      k,           & ! Vertical loop index [-]
      top_indx,    & ! Model-top index [-]
      bottom_indx, & ! Lowest-level index [-]
      incr           ! Top-down loop increment [-]


    !--------------------------- Begin Code ---------------------------

    allocate( damping_profile%tau_sponge_damp(ngrdcol,nz), &
              damping_profile%sponge_layer_depth(ngrdcol) )

    damping_profile%tau_sponge_damp = settings%tau_sponge_damp_max
    damping_profile%sponge_layer_depth = settings%sponge_damp_depth * gr%zm(:,gr%k_ub_zm)

    ! Check the value of tau_sponge_damp_min.
    if ( settings%tau_sponge_damp_min < two * dt ) then
       write(fstderr,*) "Error:  tau_sponge_damp_min is too small!"
       write(fstderr,*) "It must be at least 2.0 * dt"
       error stop
    endif

    if ( gr%grid_dir_indx > 0 ) then
       ! Ascending grid
       top_indx = nz ! top level
       bottom_indx = 1 ! lowest level
       incr = -1 ! loop increment
    else ! gr%grid_dir_indx < 0
       ! Descending grid
       top_indx = 1 ! top level
       bottom_indx = nz ! lowest level
       incr = 1
    endif

    ! Each column's layer depth and time scales follow its current heights.
    do k = top_indx, bottom_indx, incr
      do i = 1, ngrdcol

        if ( gr%zm(i,gr%k_ub_zm) - z(i,k) < damping_profile%sponge_layer_depth(i) ) then

          tau_sponge_damp_exponent &
            = ( gr%zm(i,gr%k_ub_zm) - z(i,k) ) / damping_profile%sponge_layer_depth(i)

          damping_profile%tau_sponge_damp(i,k) &
            = settings%tau_sponge_damp_min &
              * ( settings%tau_sponge_damp_max &
                  / settings%tau_sponge_damp_min )**tau_sponge_damp_exponent

        end if

      end do
    end do


    return

  end subroutine initialize_tau_sponge_damp_api

  !=============================================================================
  subroutine finalize_tau_sponge_damp_api( damping_profile )

    ! Description:
    ! Frees memory allocated in initialize_tau_sponge_damp_api
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
    deallocate( damping_profile%tau_sponge_damp, damping_profile%sponge_layer_depth )


    return

  end subroutine finalize_tau_sponge_damp_api

!===============================================================================

end module sponge_layer_damping
