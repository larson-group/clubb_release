!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module mean_adv

! Description:
! Module mean_adv computes the mean advection terms for all of the
! time-tendency (prognostic) equations in the CLUBB parameterization.  All of
! the mean advection terms are solved for completely implicitly, and therefore
! become part of the left-hand side of their respective equations.
!
! Function term_ma_zt_lhs handles the mean advection terms for the variables
! located at thermodynamic grid levels.  These variables are:  rtm, thlm, wp3,
! all hydrometeor species, and sclrm.
!
! Function term_ma_zm_lhs handles the mean advection terms for the variables
! located at momentum grid levels.  The variables are:  wprtp, wpthlp, wp2,
! rtp2, thlp2, rtpthlp, up2, vp2, wpsclrp, sclrprtp, sclrpthlp, and sclrp2.

  implicit none

  private ! Default scope

  public :: term_ma_zt_lhs, & 
            term_ma_zm_lhs

  contains

!===============================================================================
  pure function term_ma_zt_lhs( wm_zt, dzt, level ) & 
  result( lhs )

! Description:
! Mean advection of var_zt:  implicit portion of the code.
!
! The variable "var_zt" stands for a variable that is located at thermodynamic
! grid levels.
!
! The d(var_zt)/dt equation contains a mean advection term:
!
! - w * d(var_zt)/dz.
!
! This term is solved for completely implicitly, such that:
!
! - w * d( var_zt(t+1) )/dz.
!
! Note:  When the term is brought over to the left-hand side, the sign is
!        reversed and the leading "-" in front of the term is changed to a "+".
!
! The timestep index (t+1) means that the value of var_zt being used is from the
! next timestep, which is being advanced to in solving the d(var_zt)/dt
! equation.
!
! This term is discretized as follows:
!
! The values of var_zt are found on the thermodynamic levels, as are the values
! of wm_zt (mean vertical velocity on thermodynamic levels).  The variable
! var_zt is interpolated to the intermediate momentum levels.  The derivative of
! the interpolated values is taken over the central thermodynamic level.  The
! derivative is multiplied by wm_zt at the central thermodynamic level to get
! the desired result.
!
! -----var_ztp1-------------------------------------------- t(k+1)
!
! ================var_zt(interp)=========================== m(k)
!
! -----var_zt---------------------d(var_zt)/dz-----wm_zt--- t(k)
!
! ================var_zt(interp)=========================== m(k-1)
!
! -----var_ztm1-------------------------------------------- t(k-1)
!
! The vertical indices t(k+1), m(k), t(k), m(k-1), and t(k-1) correspond with
! altitudes zt(k+1), zm(k), zt(k), zm(k-1), and zt(k-1), respectively.  The
! letter "t" is used for thermodynamic levels and the letter "m" is used for
! momentum levels.
!
! dzt(k) = 1 / ( zm(k) - zm(k-1) )

! References:
!   None
!-----------------------------------------------------------------------

    use grid_class, only: & 
        gr ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1,    & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2,    & ! Thermodynamic main diagonal index.
      km1_tdiag = 3       ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1,    & ! Index for upper thermodynamic level grid weight.
      t_below = 2       ! Index for lower thermodynamic level grid weight.

    ! Input Variables
    real, intent(in) :: & 
      wm_zt,   & ! wm_zt(k)                        [m/s]
      dzt        ! Inverse of grid spacing (k)   [1/m]

    integer, intent(in) :: & 
      level ! Central thermodynamic level (on which calculation occurs).

    ! Return Variable
    real, dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      mk,    & ! Momentum level directly above central thermodynamic level.
      mkm1     ! Momentum level directly below central thermodynamic level.

    ! Momentum level (k) is between thermodynamic level (k+1)
    ! and thermodynamic level (k).
    mk = level

    ! Momentum level (k-1) is between thermodynamic level (k)
    ! and thermodynamic level (k-1).
    mkm1 = level - 1

    if ( level == 1 ) then

      ! k = 1 (bottom level); lower boundary level.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) & 
      = 0.0

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
      lhs(k_tdiag) & 
      = 0.0

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) & 
      = 0.0


    elseif ( level > 1 .and. level < gr%nnzp ) then

      ! Most of the interior model; normal conditions.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) & 
      = + wm_zt * dzt * gr%weights_zt2zm(t_above,mk)

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
      lhs(k_tdiag) & 
      = + wm_zt * dzt * (   gr%weights_zt2zm(t_below,mk) & 
                        - gr%weights_zt2zm(t_above,mkm1)   )

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) & 
      = - wm_zt * dzt * gr%weights_zt2zm(t_below,mkm1)


    elseif ( level == gr%nnzp ) then

      ! k = gr%nnzp (top level); upper boundary level.

      ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
      lhs(kp1_tdiag) & 
      = 0.0

      ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
      lhs(k_tdiag) & 
      = 0.0

      ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
      lhs(km1_tdiag) & 
      = 0.0


    endif


    return
  end function term_ma_zt_lhs

!===============================================================================
  pure function term_ma_zm_lhs( wm_zm, dzm, level ) & 
  result( lhs )

! Description:
! Mean advection of var_zm:  implicit portion of the code.
!
! The variable "var_zm" stands for a variable that is located at momentum grid
! levels.
!
! The d(var_zm)/dt equation contains a mean advection term:
!
! - w * d(var_zm)/dz.
!
! This term is solved for completely implicitly, such that:
!
! - w * d( var_zm(t+1) )/dz.
!
! Note:  When the term is brought over to the left-hand side, the sign is
!        reversed and the leading "-" in front of the term is changed to a "+".
!
! The timestep index (t+1) means that the value of var_zm being used is from the
! next timestep, which is being advanced to in solving the d(var_zm)/dt
! equation.
!
! This term is discretized as follows:
!
! The values of var_zm are found on the momentum levels, as are the values of
! wm_zm (mean vertical velocity on momentum levels).  The variable var_zm is
! interpolated to the intermediate thermodynamic levels.  The derivative of the
! interpolated values is taken over the central momentum level.  The derivative
! is multiplied by wm_zm at the central momentum level to get the desired
! result.
!
! =====var_zmp1============================================ m(k+1)
!
! ----------------var_zm(interp)--------------------------- t(k+1)
!
! =====var_zm=====================d(var_zm)/dz=====wm_zm=== m(k)
!
! ----------------var_zm(interp)--------------------------- t(k)
!
! =====var_zmm1============================================ m(k-1)
!
! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond with
! altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.  The
! letter "t" is used for thermodynamic levels and the letter "m" is used for
! momentum levels.
!
! dzm(k) = 1 / ( zt(k+1) - zt(k) )

! References:
!-----------------------------------------------------------------------

    use grid_class, only: & 
        gr

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1,    & ! Momentum superdiagonal index.
      k_mdiag   = 2,    & ! Momentum main diagonal index.
      km1_mdiag = 3       ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1,    & ! Index for upper momentum level grid weight.
      m_below = 2       ! Index for lower momentum level grid weight.

    ! Input Variables
    real, intent(in) :: & 
      wm_zm,   & ! wm_zm(k)                        [m/s]
      dzm        ! Inverse of grid spacing (k)     [1/m]

    integer, intent(in) :: & 
      level ! Central momentum level (on which calculation occurs).

    ! Return Variable
    real, dimension(3) :: lhs

    ! Local Variables
    integer :: & 
      tkp1,  & ! Thermodynamic level directly above central momentum level.
      tk       ! Thermodynamic level directly below central momentum level.

    ! Thermodynamic level (k+1) is between momentum level (k+1)
    ! and momentum level (k).
    tkp1 = level + 1

    ! Thermodynamic level (k) is between momentum level (k)
    ! and momentum level (k-1).
    tk = level

    if ( level == 1 ) then

      ! k = 1; lower boundery level at surface.

      ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
      lhs(kp1_mdiag) & 
      = 0.0

      ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
      lhs(k_mdiag) & 
      = 0.0

      ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
      lhs(km1_mdiag) & 
      = 0.0


    elseif ( level > 1 .and. level < gr%nnzp ) then

      ! Most of the interior model; normal conditions.

      ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
      lhs(kp1_mdiag) & 
      = + wm_zm * dzm * gr%weights_zm2zt(m_above,tkp1)

      ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
      lhs(k_mdiag) & 
      = + wm_zm * dzm * (   gr%weights_zm2zt(m_below,tkp1) & 
                        - gr%weights_zm2zt(m_above,tk)  )

      ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
      lhs(km1_mdiag) & 
      = - wm_zm * dzm * gr%weights_zm2zt(m_below,tk)


    elseif ( level == gr%nnzp ) then

      ! k = gr%nnzp (top level); upper boundary level.

      ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
      lhs(kp1_mdiag) & 
      = 0.0

      ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
      lhs(k_mdiag) & 
      = 0.0

      ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
      lhs(km1_mdiag) & 
      = 0.0


    endif

    return
  end function term_ma_zm_lhs

!===============================================================================

end module mean_adv
