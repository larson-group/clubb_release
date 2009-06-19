$Id$
module damping
!
! This module is used for damping variables in upper altitudes of the grid.
!
!---------------------------------------------------------------------------------------------------
  implicit none

  public :: damp_xm, initialize_tau_damp, finalize_tau_damp

  real, public :: &
    tau_damp_min, & ! Minimum damping time-scale (at the top)
    tau_damp_max, & ! maxim damping time-scale (base of damping layer)
    damp_depth      ! damping depth as a fraction of domain height

  logical, public :: &
    l_damping       ! True if damping is being used


  real, private, allocatable, dimension(:) :: &
    tau ! Damping factor

  integer, private :: &
    n_damp ! Number of levels damped

  private

  contains

  !-------------------------------------------------------------------------------------------------
  function damp_xm( dt, xm_ref, xm ) result( xm_p )
    !
    !  Description: Damps specified variable. The module must be initialized for
    !  this function to work. Otherwise a stop is issued.
    !
    !-----------------------------------------------------------------------------------------------

  !  "Spange"-layer damping at the domain top region

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

    if( allocated( tau ) ) then

      xm_p = xm

      do k = gr%nnzp, gr%nnzp-n_damp, -1

        xm_p(k) = xm(k) - (( xm(k) - xm_ref(k) ) / tau(k) ) * dt
  
      end do ! k

    else

      stop "tau in damping used before initialization"

    end if

  end function damp_xm

  !-------------------------------------------------------------------------------------------------
  subroutine initialize_tau_damp( dt )
    !
    !  Description:
    !    Initialize tau used for damping
    !
    !
    !-----------------------------------------------------------------------------------------------
    use stats_precision, only: time_precision ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    implicit none

    ! Input Variable(s)
    real(kind=time_precision), intent(in) :: dt ! Model Timestep [s}

    integer k

    allocate(tau(1:gr%nnzp))

    if( tau_damp_min < 2 * dt) then
      print*,'Error: in damping() tau_damp_min is too small!'
      stop
    end if

    do k=gr%nnzp,1,-1
      if(gr%zt(gr%nnzp)-gr%zt(k).lt.damp_depth*gr%zt(gr%nnzp)) then
        n_damp=gr%nnzp-k+1
      endif
    end do

    do k=gr%nnzp,gr%nnzp-n_damp,-1
      tau(k) = tau_damp_min *(tau_damp_max/tau_damp_min)** &
                 ((gr%zt(gr%nnzp)-gr%zt(k))/(gr%zt(gr%nnzp)-gr%zt(gr%nnzp-n_damp)))
    end do

  end subroutine initialize_tau_damp

  !-------------------------------------------------------------------------------------------------
  subroutine finalize_tau_damp()
    !
    !  Description:
    !    Frees memory allocated in initialize_tau_damp
    !
    !-----------------------------------------------------------------------------------------------
    implicit none

    deallocate( tau )

  end subroutine finalize_tau_damp


end module damping
