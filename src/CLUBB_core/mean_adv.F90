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

#ifdef GPTL
  use gptl
#endif

  implicit none

  private ! Default scope

  public :: term_ma_zt_lhs, & 
            term_ma_zm_lhs

  integer, parameter :: &
    ndiags3 = 3

  integer :: ret_code

  contains

  !=============================================================================
  subroutine term_ma_zt_lhs( nzm, nzt, ngrdcol, wm_zt, weights_zt2zm, & ! Intent(in)
                             invrs_dzt, invrs_dzm,     & ! Intent(in)
                             l_upwind_xm_ma,           & ! Intent(in)
                             lhs_ma )                    ! Intent(out)

    ! Description:
    ! Mean advection of var_zt:  implicit portion of the code.
    !
    ! The variable "var_zt" stands for a variable that is located at
    ! thermodynamic grid levels.
    !
    ! The d(var_zt)/dt equation contains a mean advection term:
    !
    ! - w * d(var_zt)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - w * d( var_zt(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! The timestep index (t+1) means that the value of var_zt being used is from
    ! the next timestep, which is being advanced to in solving the d(var_zt)/dt
    ! equation.
    !
    ! This term is discretized as follows:
    !
    ! The values of var_zt are found on the thermodynamic levels, as are the
    ! values of wm_zt (mean vertical velocity on thermodynamic levels).  The
    ! variable var_zt is interpolated to the intermediate momentum levels.  The
    ! derivative of the interpolated values is taken over the central
    ! thermodynamic level.  The derivative is multiplied by wm_zt at the central
    ! thermodynamic level to get the desired result.
    !
    ! -----var_zt---------------------------------------------- t(k+1)
    !
    ! =================var_zt(interp)========================== m(k+1)
    !
    ! -----var_zt---------------------d(var_zt)/dz-----wm_zt--- t(k)
    !
    ! =================var_zt(interp)========================== m(k)
    !
    ! -----var_zt---------------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+1), m(k+1), t(k), m(k), and t(k-1) correspond
    ! with altitudes zt(k+1), zm(k+1), zt(k), zm(k), and zt(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) )
    !
    !
    ! Special discretization for upper and lower boundary levels:
    !
    ! Zero derivative method: the derivative d(var_zt)/dz is set to 0 over the
    ! model top and the model bottom.
    !
    ! This method corresponds with the "zero-flux" boundary condition option
    ! for eddy diffusion, where d(var_zt)/dz is set to 0 across the upper
    ! boundary.
    !
    ! In order to discretize the upper boundary condition, consider a new level
    ! outside the model (thermodynamic level gr%nzt+1) just above the upper
    ! boundary level (thermodynamic level gr%nzt).  The value of var_zt at the
    ! level just outside the model is defined to be the same as the value of
    ! var_zt at thermodynamic level gr%nzt.  Therefore, the value of
    ! d(var_zt)/dz between the level just outside the model and the uppermost
    ! thermodynamic level is 0, staying consistent with the zero-flux boundary
    ! condition option for the eddy diffusion portion of the code.  Therefore,
    ! the value of var_zt at momentum level gr%nzm, which is the upper boundary
    ! of the model, would be the same as the value of var_zt at the uppermost
    ! thermodynamic level.
    !
    ! The values of var_zt are found on the thermodynamic levels, as are the
    ! values of wm_zt (mean vertical velocity on the thermodynamic levels).  The
    ! variable var_zt is interpolated to momentum level gr%nzm-1, based on
    ! the values of var_zt at thermodynamic levels gr%nzt and gr%nzt-1.  The
    ! value of var_zt at momentum level gr%nzm is set equal to the value of
    ! var_zt at thermodynamic level gr%nzt, as described above.  The derivative
    ! of the set and interpolated values, d(var_zt)/dz, is taken over the
    ! central thermodynamic level.  The derivative is multiplied by wm_zt at the
    ! central thermodynamic level to get the desired result.
    !
    ! For the following diagram, k = gr%nzm, which is the uppermost level of
    ! the model:
    !
    ! --[var_zt(k+1) = var_zt(k)]----(level outside model)----- t(k+1)
    !
    ! ==[var_zt(top) = var_zt(k)]===[d(var_zt)/dz|_(top) = 0]== m(k+1) Boundary
    !
    ! -----var_zt(k)------------------d(var_zt)/dz-----wm_zt--- t(k)
    !
    ! =================var_zt(interp)========================== m(k)
    !
    ! -----var_zt(k-1)----------------------------------------- t(k-1)
    !
    ! where (top) stands for the grid index of momentum level k = gr%nzm, which
    ! is the upper boundary of the model.
    !
    ! An analogous discretization occurs at the lower boundary.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_tdiag = 1, & ! Thermodynamic superdiagonal index.
      k_tdiag   = 2, & ! Thermodynamic main diagonal index.
      km1_tdiag = 3    ! Thermodynamic subdiagonal index.

    integer, parameter :: & 
      t_above = 1, & ! Index for upper thermodynamic level grid weight.
      t_below = 2    ! Index for lower thermodynamic level grid weight.

    ! -------------------------- Input Variables --------------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: & 
      wm_zt,     & ! wm_zt                        [m/s]
      invrs_dzt    ! Inverse of grid spacing      [1/m]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: & 
      invrs_dzm    ! Inverse of grid spacing      [1/m]
      
    real( kind = core_rknd ), dimension(ngrdcol,nzm,t_above:t_below), intent(in) :: & 
      weights_zt2zm

    logical, intent(in) :: &
      l_upwind_xm_ma    ! This flag determines whether we want to use an upwind
                        ! differencing approximation rather than a centered
                        ! differencing for turbulent or mean advection terms.
                        ! It affects rtm, thlm, sclrm, um and vm.

    ! -------------------------- Return Variable --------------------------
    real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzt), intent(out) :: &
      lhs_ma    ! Mean advection contributions to lhs    [1/s]

    ! -------------------------- Local Variables --------------------------

    integer :: i, k    ! Vertical level index

    !-------------------------- Begin Code --------------------------

#ifdef GPTL
    ret_code = GPTLstart('ik_loops')
#endif

    if ( .not. l_upwind_xm_ma ) then  ! Use centered differencing

      ! Most of the interior model; normal conditions.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 2, nzt-1, 1
        do i = 1, ngrdcol

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,i,k) = + wm_zt(i,k) * invrs_dzt(i,k) * weights_zt2zm(i,k+1,t_above)

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,i,k) = + wm_zt(i,k) * invrs_dzt(i,k) &
                                  * ( weights_zt2zm(i,k+1,t_below) - weights_zt2zm(i,k,t_above) )

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,i,k) = - wm_zt(i,k) * invrs_dzt(i,k) * weights_zt2zm(i,k,t_below)
        end do
      end do ! k = 2, nzt-1, 1
      !$acc end parallel loop

      ! Upper Boundary

      ! Special discretization for zero derivative method, where the
      ! derivative d(var_zt)/dz over the model top is set to 0, in order
      ! to stay consistent with the zero-flux boundary condition option
      ! in the eddy diffusion code.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        
        ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
        lhs_ma(kp1_tdiag,i,nzt) = zero
        
        ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
        lhs_ma(k_tdiag,i,nzt) = + wm_zt(i,nzt) * invrs_dzt(i,nzt) &
                                  * ( one - weights_zt2zm(i,nzm-1,t_above) )

        ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
        lhs_ma(km1_tdiag,i,nzt) = - wm_zt(i,nzt) * invrs_dzt(i,nzt) * weights_zt2zm(i,nzm-1,t_below)

      end do
      !$acc end parallel loop

      ! Lower Boundary

      ! Special discretization for zero derivative method, where the
      ! derivative d(var_zt)/dz over the model bottom is set to 0.
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        
        ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
        lhs_ma(kp1_tdiag,i,1) = + wm_zt(i,1) * invrs_dzt(i,1) * weights_zt2zm(i,2,t_above)
        
        ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
        lhs_ma(k_tdiag,i,1) = - wm_zt(i,1) * invrs_dzt(i,1) &
                                * ( one - weights_zt2zm(i,2,t_below) )

        ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
        lhs_ma(km1_tdiag,i,1) = zero

      end do
      !$acc end parallel loop

    else ! l_upwind_xm_ma == .true.; use "upwind" differencing

      ! Most of the interior model; normal conditions.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 2, nzt-1, 1
        do i = 1, ngrdcol
          if ( wm_zt(i,k) >= zero ) then  ! Mean wind is in upward direction

             ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
             lhs_ma(kp1_tdiag,i,k) = zero

             ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
             lhs_ma(k_tdiag,i,k) = + wm_zt(i,k) * invrs_dzm(i,k)

             ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
             lhs_ma(km1_tdiag,i,k) = - wm_zt(i,k) * invrs_dzm(i,k)
             
          else  ! wm_zt < 0; Mean wind is in downward direction

             ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
             lhs_ma(kp1_tdiag,i,k) = + wm_zt(i,k) * invrs_dzm(i,k+1)

             ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
             lhs_ma(k_tdiag,i,k) = - wm_zt(i,k) * invrs_dzm(i,k+1)

             ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
             lhs_ma(km1_tdiag,i,k) = zero

          endif ! wm_zt > 0
          
        end do
      end do ! k = 2, nzt-1, 1
      !$acc end parallel loop

      ! Upper Boundary
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        if ( wm_zt(i,nzt) >= zero ) then  ! Mean wind is in upward direction

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,i,nzt) = zero

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,i,nzt) = + wm_zt(i,nzt) * invrs_dzm(i,nzm-1)

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,i,nzt) = - wm_zt(i,nzt) * invrs_dzm(i,nzm-1)

        else  ! wm_zt < 0; Mean wind is in downward direction

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,i,nzt) = zero

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,i,nzt) = zero

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,i,nzt) = zero

        end if ! wm_zt > 0
      end do
      !$acc end parallel loop

      ! Lower Boundary
      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        if ( wm_zt(i,1) >= zero ) then  ! Mean wind is in upward direction

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,i,1) = zero

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,i,1) = zero

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,i,1) = zero

        else  ! wm_zt < 0; Mean wind is in downward direction

          ! Thermodynamic superdiagonal: [ x var_zt(k+1,<t+1>) ]
          lhs_ma(kp1_tdiag,i,1) = + wm_zt(i,1) * invrs_dzm(i,2)

          ! Thermodynamic main diagonal: [ x var_zt(k,<t+1>) ]
          lhs_ma(k_tdiag,i,1) = - wm_zt(i,1) * invrs_dzm(i,2)

          ! Thermodynamic subdiagonal: [ x var_zt(k-1,<t+1>) ]
          lhs_ma(km1_tdiag,i,1) = zero

        end if ! wm_zt > 0
      end do
      !$acc end parallel loop

    endif ! l_upwind_xm_ma

#ifdef GPTL
    !$acc wait
    ret_code = GPTLstop('ik_loops')
#endif

    return

  end subroutine term_ma_zt_lhs

  !=============================================================================
  subroutine term_ma_zm_lhs( nzm, nzt, ngrdcol, wm_zm, &
                             invrs_dzm, weights_zm2zt, & 
                             lhs_ma )

    ! Description:
    ! Mean advection of var_zm:  implicit portion of the code.
    !
    ! The variable "var_zm" stands for a variable that is located at momentum
    ! grid levels.
    !
    ! The d(var_zm)/dt equation contains a mean advection term:
    !
    ! - w * d(var_zm)/dz.
    !
    ! This term is solved for completely implicitly, such that:
    !
    ! - w * d( var_zm(t+1) )/dz.
    !
    ! Note:  When the term is brought over to the left-hand side, the sign
    !        is reversed and the leading "-" in front of the term is changed to
    !        a "+".
    !
    ! The timestep index (t+1) means that the value of var_zm being used is from
    ! the next timestep, which is being advanced to in solving the d(var_zm)/dt
    ! equation.
    !
    ! This term is discretized as follows:
    !
    ! The values of var_zm are found on the momentum levels, as are the values
    ! of wm_zm (mean vertical velocity on momentum levels).  The variable var_zm
    ! is interpolated to the intermediate thermodynamic levels.  The derivative
    ! of the interpolated values is taken over the central momentum level.  The
    ! derivative is multiplied by wm_zm at the central momentum level to get the
    ! desired result.
    !
    ! =====var_zm============================================== m(k+1)
    !
    ! -----------------var_zm(interp)-------------------------- t(k)
    !
    ! =====var_zm=====================d(var_zm)/dz=====wm_zm=== m(k)
    !
    ! -----------------var_zm(interp)-------------------------- t(k-1)
    !
    ! =====var_zm============================================== m(k-1)
    !
    ! The vertical indices m(k+1), t(k), m(k), t(k-1), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k), zm(k), zt(k-1), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k) - zt(k-1) )

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

    use constants_clubb, only: &
        zero ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: & 
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    integer, parameter :: & 
      m_above = 1, & ! Index for upper momentum level grid weight.
      m_below = 2    ! Index for lower momentum level grid weight.

    ! -------------------------- Input Variables --------------------------
    integer, intent(in) :: &
      nzm, &
      nzt, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: & 
      wm_zm,     & ! wm_zm                        [m/s]
      invrs_dzm    ! Inverse of grid spacing      [1/m]
      
    real( kind = core_rknd ), dimension(ngrdcol,nzt,m_above:m_below), intent(in) :: & 
      weights_zm2zt

    ! -------------------------- Return Variable --------------------------
    real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nzm), intent(out) :: &
      lhs_ma    ! Mean advection contributions to lhs  [1/s]

    ! -------------------------- Local Variables
    integer :: i, k, b    ! Vertical level index

    ! -------------------------- Begin Code --------------------------

    ! Set lower boundary array to 0
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do b = 1, ndiags3
        lhs_ma(b,i,1) = zero
      end do
    end do
    !$acc end parallel loop

    ! Most of the interior model; normal conditions.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 2, nzm-1, 1
      do i = 1, ngrdcol
        
        ! Momentum superdiagonal: [ x var_zm(k+1,<t+1>) ]
        lhs_ma(kp1_mdiag,i,k) = + wm_zm(i,k) * invrs_dzm(i,k) * weights_zm2zt(i,k,m_above)

        ! Momentum main diagonal: [ x var_zm(k,<t+1>) ]
        lhs_ma(k_mdiag,i,k) = + wm_zm(i,k) * invrs_dzm(i,k) * ( weights_zm2zt(i,k,m_below) & 
                                       - weights_zm2zt(i,k-1,m_above) )

        ! Momentum subdiagonal: [ x var_zm(k-1,<t+1>) ]
        lhs_ma(km1_mdiag,i,k) = - wm_zm(i,k) * invrs_dzm(i,k) * weights_zm2zt(i,k-1,m_below)
        
      end do
    end do ! k = 2, nzm-1, 1
    !$acc end parallel loop

    ! Set upper boundary array to 0
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do b = 1, ndiags3
        lhs_ma(b,i,nzm) = zero
      end do
    end do
    !$acc end parallel loop

    return

  end subroutine term_ma_zm_lhs

  !=============================================================================

end module mean_adv
