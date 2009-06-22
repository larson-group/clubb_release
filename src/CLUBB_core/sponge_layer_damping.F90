!$Id$
module sponge_layer_damping
!
! This module is used for damping variables in upper altitudes of the grid.
!
!---------------------------------------------------------------------------------------------------
  implicit none

  public :: sponge_damp_xm, initialize_tau_sponge_damp, finalize_tau_sponge_damp

  real, public :: &
    tau_sponge_damp_min, & ! Minimum damping time-scale (at the top)
    tau_sponge_damp_max, & ! maxim damping time-scale (base of damping layer)
    sponge_damp_depth      ! damping depth as a fraction of domain height

  logical, public :: &
    l_sponge_damping       ! True if damping is being used


  real, private, allocatable, dimension(:) :: &
    tau_sponge_damp ! Damping factor

  integer, private :: &
    n_sponge_damp ! Number of levels damped

  private

  contains

  !-------------------------------------------------------------------------------------------------
  function sponge_damp_xm( dt, xm_ref, xm ) result( xm_p )
    !
    !  Description: Damps specified variable. The module must be initialized for
    !  this function to work. Otherwise a stop is issued.
    !
    !-----------------------------------------------------------------------------------------------

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

    ! Output Variable(s)
    real, dimension(gr%nnzp) :: xm_p ! Variable damped [-]

    integer k

    if( allocated( tau_sponge_damp ) ) then

      xm_p = xm

      do k = gr%nnzp, gr%nnzp-n_sponge_damp, -1

        xm_p(k) = xm(k) - real( ( ( xm(k) - xm_ref(k) ) / tau_sponge_damp(k) ) * dt )
  
      end do ! k

    else

      stop "tau_sponge_damp in damping used before initialization"

    end if

  end function sponge_damp_xm

  !-------------------------------------------------------------------------------------------------
  subroutine initialize_tau_sponge_damp( dt )
    !
    !  Description:
    !    Initialize tau_sponge_damp used for damping
    !
    !
    !-----------------------------------------------------------------------------------------------
    use stats_precision, only: time_precision ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    implicit none

    ! Input Variable(s)
    real(kind=time_precision), intent(in) :: dt ! Model Timestep [s}

    integer k

    allocate(tau_sponge_damp(1:gr%nnzp))

    if( tau_sponge_damp_min < 2 * dt) then
      print*,'Error: in damping() tau_sponge_damp_min is too small!'
      stop
    end if

    do k=gr%nnzp,1,-1
      if(gr%zt(gr%nnzp)-gr%zt(k) < sponge_damp_depth*gr%zt(gr%nnzp)) then
        n_sponge_damp=gr%nnzp-k+1
      endif
    end do

    do k=gr%nnzp,gr%nnzp-n_sponge_damp,-1
      tau_sponge_damp(k) = tau_sponge_damp_min *(tau_sponge_damp_max/tau_sponge_damp_min)** &
                 ((gr%zt(gr%nnzp)-gr%zt(k))/(gr%zt(gr%nnzp)-gr%zt(gr%nnzp-n_sponge_damp)))
    end do

  end subroutine initialize_tau_sponge_damp

  !-------------------------------------------------------------------------------------------------
  subroutine finalize_tau_sponge_damp()
    !
    !  Description:
    !    Frees memory allocated in initialize_tau_sponge_damp
    !
    !-----------------------------------------------------------------------------------------------
    implicit none

    deallocate( tau_sponge_damp )

  end subroutine finalize_tau_sponge_damp


end module sponge_layer_damping
