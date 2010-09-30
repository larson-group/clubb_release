!----------------------------------------------------------------------
! $Id$
module jun25

!       Description:
!       Contains subroutines for the June 11 Altocumulous case.
!----------------------------------------------------------------------

  implicit none

  public :: jun25_altocu_tndcy

  private ! Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine jun25_altocu_tndcy( time, time_initial, & 
                                 wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
                                 sclrm_forcing, edsclrm_forcing )

! Description:
!   Computes subsidence, radiation, and LS tendencies for the June
!   25th Altocumulous case.

! References:
!   None
!-----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use constants_clubb, only: pi, Cp, Lv, zero_threshold, fstderr ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use interpolation, only: linear_interpolation ! Procedure(s)

    use error_code, only: clubb_at_least_debug_level ! Procedure(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use interpolation, only: lin_int ! Procedure(s)

    implicit none

    ! Constant parameters

    ! Input variables
    real(kind=time_precision), intent(in) :: & 
    time,          & ! Time of simulation since start  [s]
    time_initial     ! Initial time of simulation      [s]

    ! Output variables
    real, dimension(gr%nnzp), intent(inout) ::  & 
    wm_zt,           & ! Vertical ascent/descent on therm. grid      [m/s]
    wm_zm,           & ! Vertical ascent/descent on moment. grid     [m/s]
    thlm_forcing,    & ! Change in liq. water potential temperature 
                     ! due to radiative heating and ice diffusion  [K/s]
    rtm_forcing        ! Change in total water due to ice diffusion  [kg/kg/s]

    real, dimension(gr%nnzp,sclr_dim),intent(out) ::  & 
    sclrm_forcing ! Large-scale tendency for passive scalars      [units/s]

    real, dimension(gr%nnzp,edsclr_dim),intent(out) ::  & 
    edsclrm_forcing ! Large-scale tendency for passive scalars    [units/s]

    !---------------------------------------------------------------
    ! Working arrays for subsidence interpolation
    !---------------------------------------------------------------
    real, dimension(5) :: &
    zsubs ! Heights at which wm_zt data is supplied (used for subsidence interpolation) [m]

    real, dimension(6) :: & 
      tsubs ! times after initialization at which wm_zt data is supplied [s]

    real, dimension(5) :: & 
    wt1, wt2, wt3, wt4, wt5, wt6 ! wtX(Y) vertical velocity specified at height Y and time X [m/s]

    real, dimension(gr%nnzp) :: & 
      w1, w2 ! vertical velocity before (w1) and after (w2) 
        !      the current time at the specified level [m/s]

    integer :: k

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    zsubs(1) = gr%zm(1)
    zsubs(2) = gr%zm(1) + 360.
    zsubs(3) = gr%zm(1) + 1090.
    zsubs(4) = gr%zm(1) + 1890.
    zsubs(5) = gr%zm(1) + 2500.

    tsubs(1) = 0
    tsubs(2) = 10800
    tsubs(3) = 28800
    tsubs(4) = 36000
    tsubs(5) = 36000
    tsubs(6) = 36000

    wt1(1) = 0.
    wt1(2) = .004
    wt1(3) = .004
    wt1(4) = .004
    wt1(5) = 0.

    wt2(1) = 0.
    wt2(2) = .004
    wt2(3) = .004
    wt2(4) = .004
    wt2(5) = 0.

    wt3(1) = 0.
    wt3(2) = -.003
    wt3(3) = -.003
    wt3(4) = -.003
    wt3(5) = 0.

    wt4(1) = 0.
    wt4(2) = -.003
    wt4(3) = -.003
    wt4(4) = -.003
    wt4(5) = 0.

    wt5(1) = 0.
    wt5(2) = -.003
    wt5(3) = -.003
    wt5(4) = -.003
    wt5(5) = 0.

    wt6(1) = 0.
    wt6(2) = -.003
    wt6(3) = -.003
    wt6(4) = -.003
    wt6(5) = 0.

    !---------------------------------------------------------------
    ! Using linear interpolation to calculate subsidence
    ! Original code by Michael Falk
    ! Added for Jun.25 case by Adam Smith, 13 April 2006
    !---------------------------------------------------------------
    if ( (time - time_initial) < tsubs(1) ) then
      do k=1,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt1,gr%zt(k),wm_zt(k))
        else
          wm_zt(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
      end do

    else if ( (time - time_initial) < tsubs(2)) then
      do k=2,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt1,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt2,gr%zt(k),w2(k))
        else
          w1(k) = 0.0
          w2(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
        !wm_zt(k) = &
        !  real(((time-time_initial)-tsubs(1)) &
        !         /(tsubs(2)-tsubs(1))*(w2(k)-w1(k))+w1(k))
        wm_zt(k) = lin_int( real(time-time_initial), tsubs(2), tsubs(1), w2(k) ,w1(k) )
      end do

    else if ( (time - time_initial) < tsubs(3)) then
      do k=2,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt2,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt3,gr%zt(k),w2(k))
        else
          w1(k) = 0.0
          w2(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
        !wm_zt(k) =  &
        !  real(((time-time_initial)-tsubs(2)) &
        !         /(tsubs(3)-tsubs(2))*(w2(k)-w1(k))+w1(k))
        wm_zt(k) = lin_int(real(time-time_initial), tsubs(3), tsubs(2), w2(k), w1(k) )
      end do

    else if ( (time - time_initial) < tsubs(4)) then
      do k=2,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt3,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt4,gr%zt(k),w2(k))
        else
          w1(k) = 0.0
          w2(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
        !wm_zt(k) =  &
        !  real(((time-time_initial)-tsubs(3)) &
        !         /(tsubs(4)-tsubs(3))*(w2(k)-w1(k))+w1(k))
        wm_zt(k) = lin_int( real(time-time_initial), tsubs(4), tsubs(3), w2(k), w1(k) )
      end do

    else if ( (time - time_initial) < tsubs(5)) then
      do k=2,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt4,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt5,gr%zt(k),w2(k))
        else
          w1(k) = 0.0
          w2(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
        !wm_zt(k) =  &
        !  real(((time-time_initial)-tsubs(4)) &
        !         /(tsubs(5)-tsubs(4))*(w2(k)-w1(k))+w1(k))
        wm_zt = lin_int( real(time-time_initial), tsubs(5), tsubs(4), w2(k), w1(k) )
      end do

    else if ( (time - time_initial) < tsubs(6)) then
      do k=2,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt5,gr%zt(k),w1(k))
          call linear_interpolation(5,zsubs,wt6,gr%zt(k),w2(k))
        else
          w1(k) = 0.0
          w2(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
        !wm_zt(k) = &
        !  real(((time-time_initial)-tsubs(5)) &
        !         /(tsubs(6)-tsubs(5))*(w2(k)-w1(k))+w1(k))
        wm_zt = lin_int( real(time-time_initial), tsubs(6), tsubs(5), w2(k), w1(k) )
      end do

    else if ( (time - time_initial) >= tsubs(6)) then
      do k=2,gr%nnzp
        if ( gr%zt(k) <= zsubs(5) ) then
          call linear_interpolation(5,zsubs,wt6,gr%zt(k),wm_zt(k))
        else
          wm_zt(k) = 0.0
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
                             gr%zt(k), "m. is above the highest level ",  &
                             "specified in the subsidence sounding, which ",  &
                             "is at height", zsubs(5), "m."
            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
                             "this altitude."
          endif
        endif
      end do
    end if

    wm_zt(1) = wm_zt(2)

    wm_zm = zt2zm(wm_zt)

    rtm_forcing(1:gr%nnzp) = 0.
    thlm_forcing(1:gr%nnzp) = 0.

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine jun25_altocu_tndcy

end module jun25
