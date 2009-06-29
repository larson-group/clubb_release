!$Id$
module sponge_layer_damping
!
! This module is used for damping variables in upper altitudes of the grid.
!
!---------------------------------------------------------------------------------------------------
  implicit none

  public :: sponge_damp_xm, initialize_tau_sponge_damp, finalize_tau_sponge_damp, &
            sponge_damp_settings, sponge_damp_profile


  type sponge_damp_settings

    real :: &
      tau_sponge_damp_min, & ! Minimum damping time-scale (at the top) [s]
      tau_sponge_damp_max, & ! Maximum damping time-scale (base of damping layer) [s]
      sponge_damp_depth      ! damping depth as a fraction of domain height [-]

    logical :: &
      l_sponge_damping       ! True if damping is being used

  end type sponge_damp_settings

  type sponge_damp_profile
    real, allocatable, dimension(:) :: &
      tau_sponge_damp ! Damping factor

    integer :: &
      n_sponge_damp ! Number of levels damped

  end type sponge_damp_profile


  type(sponge_damp_settings), public :: &
    thlm_sponge_damp_settings, &
    rtm_sponge_damp_settings, &
    uv_sponge_damp_settings

  type(sponge_damp_profile), public :: &
    thlm_sponge_damp_profile, &
    rtm_sponge_damp_profile, &
    uv_sponge_damp_profile


  private

  contains

  !---------------------------------------------------------------------------------------------
  function sponge_damp_xm( dt, xm_ref, xm, damping_profile ) result( xm_p )
    !
    !  Description: Damps specified variable. The module must be initialized for
    !  this function to work. Otherwise a stop is issued.
    !
    !-------------------------------------------------------------------------------------------

    !  "Sponge"-layer damping at the domain top region

    use grid_class, only: gr ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    implicit none

    ! Input Variable(s)
    real(kind=time_precision), intent(in) :: dt ! Model Timestep

    real, dimension(gr%nnzp), intent(in) :: &
      xm_ref ! Reference to damp to [-]

    real, dimension(gr%nnzp), intent(in) :: &
      xm ! Variable being damped [-]

    type(sponge_damp_profile), intent(in) :: &
        damping_profile

    ! Output Variable(s)
    real, dimension(gr%nnzp) :: xm_p ! Variable damped [-]

    integer k

    if( allocated( damping_profile%tau_sponge_damp ) ) then

      xm_p = xm

      do k = gr%nnzp, gr%nnzp-damping_profile%n_sponge_damp, -1

        xm_p(k) = xm(k) - real( ( ( xm(k) - xm_ref(k) ) / & 
                        damping_profile%tau_sponge_damp(k) ) * dt )

      end do ! k

    else

      stop "tau_sponge_damp in damping used before initialization"

    end if

  end function sponge_damp_xm

  !---------------------------------------------------------------------------------------------
  subroutine initialize_tau_sponge_damp( dt, settings, damping_profile )
    !
    !  Description:
    !    Initialize tau_sponge_damp used for damping
    !
    !
    !-------------------------------------------------------------------------------------------
    use stats_precision, only: time_precision ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    implicit none

    ! Input Variable(s)
    real(kind=time_precision), intent(in) :: dt ! Model Timestep [s}

    type(sponge_damp_settings), intent(in) :: &
        settings

    type(sponge_damp_profile), intent(out) :: &
      damping_profile

    integer k

    allocate( damping_profile%tau_sponge_damp(1:gr%nnzp))

    if( settings%tau_sponge_damp_min < 2 * dt) then
      print*,'Error: in damping() tau_sponge_damp_min is too small!'
      stop
    end if

    do k=gr%nnzp,1,-1
      if(gr%zt(gr%nnzp)-gr%zt(k) < settings%sponge_damp_depth*gr%zt(gr%nnzp)) then
        damping_profile%n_sponge_damp=gr%nnzp-k+1
      endif
    end do

    do k=gr%nnzp,gr%nnzp-damping_profile%n_sponge_damp,-1
      damping_profile%tau_sponge_damp(k) = settings%tau_sponge_damp_min *&
        (settings%tau_sponge_damp_max/settings%tau_sponge_damp_min)** &
        ( ( gr%zt(gr%nnzp)-gr%zt(k) ) / &
          (gr%zt(gr%nnzp) - gr%zt( gr%nnzp-damping_profile%n_sponge_damp ) ) )
    end do

  end subroutine initialize_tau_sponge_damp

  !---------------------------------------------------------------------------------------------
  subroutine finalize_tau_sponge_damp( damping_profile )
    !
    !  Description:
    !    Frees memory allocated in initialize_tau_sponge_damp
    !
    !-------------------------------------------------------------------------------------------
    implicit none

    type(sponge_damp_profile), intent(inout) :: &
        damping_profile

    deallocate( damping_profile%tau_sponge_damp )

  end subroutine finalize_tau_sponge_damp


end module sponge_layer_damping
