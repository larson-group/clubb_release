!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module turbulent_adv_pdf

  ! Description:
  ! Calculates the turbulent advection term in the predictive equation for a
  ! variance or covariance where turbulent advection is calculated by
  ! integrating over the PDF.  This includes the following predictive fields:
  ! <w'thl'>, <w'rt'>, <rt'^2>, <thl'^2>, and <rt'thl'>, as well as passive
  ! scalar fields <w'sclr'>, <sclr'^2>, <sclr'rt'>, and <sclr'thl'>.  CLUBB
  ! does not produce a PDF for horizontal wind components u and v.  However, the
  ! horizontal wind variances <u'^2> and <v'^2> still use this code, as well.

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: xpyp_term_ta_pdf_lhs,         &
            xpyp_term_ta_pdf_lhs_godunov, &
            xpyp_term_ta_pdf_rhs,         &
            xpyp_term_ta_pdf_rhs_godunov

  private    ! Set default scope

  integer, parameter :: &
    ndiags3 = 3

  contains

  !=============================================================================
  subroutine xpyp_term_ta_pdf_lhs( nz, ngrdcol, gr, coef_wpxpyp_implicit, & ! In
                                        rho_ds_zt, rho_ds_zm,            & ! In
                                        invrs_rho_ds_zm,                 & ! In
                                        l_upwind_xpyp_turbulent_adv,     & ! In
                                        sgn_turbulent_vel,               & ! In
                                        coef_wpxpyp_implicit_zm,         & ! In
                                        lhs_ta                           ) ! Out

    ! Description:
    ! Turbulent advection of <w'x'>, <x'^2>, and <x'y'>:  implicit portion of
    ! the code.
    !
    ! 1) <w'x'>
    !
    ! The d<w'x'>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'^2 x'> )/dz.
    !
    ! The value of <w'^2 x'> is found by integrating over the multivariate PDF
    ! of w and x, as detailed in function calc_wp2xp_pdf, which is found in
    ! module pdf_closure_module in pdf_closure_module.F90.
    !
    ! The equation obtained for <w'^2 x'> is written in terms of PDF parameters.
    ! Substitutions that are specific to the type of PDF used are made for the
    ! PDF parameters in order to write the <w'x'> turbulent advection term in
    ! the following form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit.
    !
    ! For the ADG1 PDF, coef_wp2xp_implicit is a_1 * < w'^3 > / < w'^2 >, where
    ! a_1 is 1 / ( 1 - sigma_sqd_w ).  The value of term_wp2xp_explicit is 0, as
    ! the <w'x'> turbulent advection term is entirely implicit.
    !
    ! For the new PDF, the calculations of both coef_wp2xp_implicit and
    ! term_wp2xp_explicit are detailed in function calc_coefs_wp2xp_semiimpl,
    ! which is found in module new_pdf in new_pdf.F90.
    !
    ! For the new hybrid PDF, the calculation of coef_wp2xp_implicit is
    ! detailed in function calc_coef_wp2xp_implicit, which is found in module
    ! new_hybrid_pdf in new_hybrid_pdf.F90.  The value of term_wp2xp_explicit
    ! is 0, as the <w'x'> turbulent advection term is entirely implicit.
    !
    ! For explicit turbulent advection, the value of coef_wp2xp_implicit is 0
    ! and the value of term_wp2xp_explicit is <w'^2 x'>, as calculated by
    ! retaining the equation for <w'^2 x'> that is written in terms of PDF
    ! parameters.  This is a general form that can be used with any type of PDF.
    !
    ! The <w'x'> turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * ( coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit ) )
    !     /dz.
    !
    ! The variable <w'x'> is evaluated at the (t+1) timestep, which allows the
    ! <w'x'> turbulent advection term to be expressed semi-implicitly as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wp2xp_implicit * <w'x'>(t+1) )/dz
    ! - (1/rho_ds) * d( rho_ds * term_wp2xp_explicit )/dz.
    !
    ! The implicit portion of <w'x'> turbulent advection term is:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wp2xp_implicit * <w'x'>(t+1) )/dz.
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of the implicit
    !        d[ ] / dz term is changed to a "+".
    !
    ! The timestep index (t+1) means that the value of <w'x'> being used is from
    ! the next timestep, which is being advanced to in solving the d<w'x'>/dt
    ! equation.
    !
    ! 2) <x'^2>
    !
    ! The d<x'^2>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'x'^2> )/dz;
    !
    ! The value of <w'x'^2> is found by integrating over the multivariate PDF of
    ! w and x, as detailed in function calc_wpxp2_pdf, which is found in module
    ! pdf_closure_module in pdf_closure_module.F90.
    !
    ! The equation obtained for <w'x'^2> is written in terms of PDF parameters.
    ! Substitutions that are specific to the type of PDF used are made for the
    ! PDF parameters in order to write the <x'^2> turbulent advection term in
    ! the following form:
    !
    ! <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit.
    !
    ! For the ADG1 PDF, the value of coef_wpxp2_implicit is
    ! (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxp2_explicit
    ! is (1-(1/3)*beta) * a_1^2 * < w'x' >^2 * < w'^3 > / < w'^2 >^2, where
    ! a_1 is 1 / ( 1 - sigma_sqd_w ).
    !
    ! For the new PDF, the calculation of coef_wpxp2_implicit is detailed in
    ! function calc_coef_wpxp2_implicit, which is found in module new_pdf in
    ! new_pdf.F90.  The value of term_wpxp2_explicit is 0, as the <x'^2>
    ! turbulent advection term is entirely implicit.
    !
    ! For the new hybrid PDF, the calculation of both coef_wpxp2_implicit and
    ! term_wpxp2_explicit are detailed in subroutine calc_coefs_wpxp2_semiimpl,
    ! which is found in module new_hybrid_pdf in new_hybrid_pdf.F90.
    !
    ! For explicit turbulent advection, the value of coef_wpxp2_implicit is 0
    ! and the value of term_wpxp2_explicit is <w'x'^2>, as calculated by
    ! retaining the equation for <w'x'^2> that is written in terms of PDF
    ! parameters.  This is a general form that can be used with any type of PDF.
    !
    ! The <x'^2> turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * ( coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit ) )
    !     /dz;
    !
    ! The variable <x'^2> is evaluated at the (t+1) timestep, which allows the
    ! <x'^2> turbulent advection term to be expressed semi-implicitly as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wpxp2_implicit * <x'^2>(t+1) )/dz;
    ! - (1/rho_ds) * d( rho_ds * term_wpxp2_explicit )/dz.
    !
    ! The implicit portion of <x'^2> turbulent advection term is:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wpxp2_implicit * <x'^2>(t+1) )/dz.
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of all implicit
    !        d[ ] / dz terms is changed to a "+".
    !
    ! The timestep index (t+1) means that the value of <x'^2> being used is from
    ! the next timestep, which is being advanced to in solving the d<x'^2>/dt
    ! equation.
    !
    ! 3) <x'y'>
    !
    ! The d<x'y'>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'x'y'> )/dz.
    !
    ! The value of <w'x'y'> is found by integrating over the multivariate PDF of
    ! w, x, and y, as detailed in function calc_wpxpyp_pdf, which is found in
    ! module pdf_closure_module in pdf_closure_module.F90.
    !
    ! The equation obtained for <w'x'y'> is written in terms of PDF parameters.
    ! Substitutions that are specific to the type of PDF used are made for the
    ! PDF parameters in order to write the <x'y'> turbulent advection term in
    ! the following form:
    !
    ! <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit.
    !
    ! For the ADG1 PDF, the value of coef_wpxpyp_implicit is
    ! (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxpyp_explicit
    ! is (1-(1/3)*beta) * a_1^2 * < w'x' > * < w'y' > * < w'^3 > / < w'^2 >^2,
    ! where a_1 is 1 / ( 1 - sigma_sqd_w ).
    !
    ! For the new PDF, the calculation of both coef_wpxpyp_implicit and
    ! term_wpxpyp_explicit are detailed in function calc_coefs_wpxpyp_semiimpl,
    ! which is found in module new_pdf in new_pdf.F90.
    !
    ! For the new hybrid PDF, the calculation of both coef_wpxpyp_implicit
    ! and term_wpxpyp_explicit are detailed in subroutine
    ! calc_coefs_wpxpyp_semiimpl, which is found in module new_hybrid_pdf in
    ! new_hybrid_pdf.F90.
    !
    ! For explicit turbulent advection, the value of coef_wpxpyp_implicit is 0
    ! and the value of term_wpxpyp_explicit is <w'x'y'>, as calculated by
    ! retaining the equation for <w'x'y'> that is written in terms of PDF
    ! parameters.  This is a general form that can be used with any type of PDF.
    !
    ! The <x'y'> turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * ( coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit ) )
    !     /dz.
    !
    ! The variable <x'y'> is evaluated at the (t+1) timestep, which allows the
    ! <x'y'> turbulent advection term to be expressed semi-implicitly as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz
    ! - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.
    !
    ! The implicit portion of <x'y'> turbulent advection term is:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz.
    !
    ! Note:  When the implicit term is brought over to the left-hand side, the
    !        sign is reversed and the leading "-" in front of all implicit
    !        d[ ] / dz terms is changed to a "+".
    !
    ! The timestep index (t+1) means that the value of <x'y'> being used is from
    ! the next timestep, which is being advanced to in solving the d<x'y'>/dt
    ! equation.
    !
    ! When x and y are the same variable, <x'y'> reduces to <x'^2> and <w'x'y'>
    ! reduces to <w'x'^2>.  Likewise, when y is set equal to w, <x'y'> becomes
    ! <w'x'> and <w'x'y'> reduces to <w'^2 x'>.  The discretization and the code
    ! used in this function will be written generally in terms of <x'y'> and
    ! coef_wpxpyp_implicit, but also applies to <x'^2> and coef_wpxp2_implicit,
    ! as well as to <w'x'> and coef_wp2xp_implicit.
    !
    ! The implicit discretization of this term is as follows:
    !
    ! 1) Centered Discretization
    !
    ! The values of <x'y'> are found on the momentum levels, while the values of
    ! coef_wpxpyp_implicit are found on the thermodynamic levels, which is where
    ! they were originally calculated by the PDF.  Additionally, the values of
    ! rho_ds_zt are found on the thermodynamic levels, and the values of
    ! invrs_rho_ds_zm are found on the momentum levels.  The values of <x'y'>
    ! are interpolated to the intermediate thermodynamic levels as <x'y'>|_zt.
    ! At the thermodynamic levels, the values of coef_wpxpyp_implicit are
    ! multiplied by <x'y'>|_zt, and their product is multiplied by rho_ds_zt.
    ! Then, the derivative (d/dz) of that expression is taken over the central
    ! momentum level, where it is multiplied by -invrs_rho_ds_zm.  This yields
    ! the desired result.
    !
    ! =xpyp============================================================== m(k+1)
    !
    ! -xpyp_zt(interp)-------coef_wpxpyp_implicit---------rho_ds_zt------ t(k+1)
    !
    ! =xpyp=d(rho_ds_zt*coef_wpxpyp_implicit*xpyp_zt)/dz=invrs_rho_ds_zm= m(k)
    !
    ! -xpyp_zt(interp)-------coef_wpxpyp_implicit---------rho_ds_zt------ t(k)
    !
    ! =xpyp============================================================== m(k-1)
    !
    ! The vertical indices m(k+1), t(k+1), m(k), t(k), and m(k-1) correspond
    ! with altitudes zm(k+1), zt(k+1), zm(k), zt(k), and zm(k-1), respectively.
    ! The letter "t" is used for thermodynamic levels and the letter "m" is used
    ! for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )
    !
    ! 2) "Upwind" Discretization.
    !
    ! The values of <x'y'> are found on the momentum levels.  The values of
    ! coef_wpxpyp_implicit are originally calculated by the PDF on the
    ! thermodynamic levels.  They are interpolated to the intermediate momentum
    ! levels as coef_wpxpyp_implicit_zm.  Additionally, the values of rho_ds_zm
    ! and the values of invrs_rho_ds_zm are found on the momentum levels.  The
    ! sign of the turbulent velocity is found on the central momentum level.  At
    ! the momentum levels, the values of coef_wpxpyp_implicit_zm are multiplied
    ! by <x'y'>, and their product is multiplied by rho_ds_zm.  Then, the
    ! derivative (d/dz) of that expression is taken.  When the sign of the
    ! turbulent velocity is positive, the "wind" is coming from below, and the
    ! derivative involves the central momentum level and the momentum level
    ! immediately below it.  When the sign of the turbulent velocity is
    ! negative, the "wind" is coming from above, and the derivative involves the
    ! central momentum level and the momenum level immediately above it.  After
    ! the derivative is taken, it is multiplied by -invrs_rho_ds_zm at the
    ! central momentum level.  This yields the desired result.
    !
    ! The turbulent velocity for <x'y'> is <w'x'y'> / <x'y'>, which has units of
    ! m/s.  The sign of the turbulent velocity is sgn( <w'x'y'> / <x'y'> ),
    ! where:
    !
    ! sgn(x) = | 1; when x >= 0
    !          | -1; when x < 0.
    !
    ! The sign of the turbulent velocity can also be rewritten as
    ! sgn( <w'x'y'> ) / sgn( <x'y'> ).  When a variance (<x'^2>) is being solved
    ! for, y = x, and sgn( <x'^2> ) is always 1.  The sign of the turbulent
    ! velocity reduces to simply sgn( <w'x'^2> ).
    !
    ! ---------coef_wpxpyp_implicit-------------------------------------- t(k+2)
    !
    ! =xpyp=======coef_wpxpyp_implicit_zm(interp.)=======rho_ds_zm======= m(k+1)
    !
    ! ---------coef_wpxpyp_implicit-------------------------------------- t(k+1)
    !
    ! =xpyp===coef_wpxpyp_implicit_zm(interp.)=rho_ds_zm=invrs_rho_ds_zm= m(k)
    !
    ! ---------coef_wpxpyp_implicit-------------------------------------- t(k)
    !
    ! =xpyp=======coef_wpxpyp_implicit_zm(interp.)=======rho_ds_zm======= m(k-1)
    !
    ! ---------coef_wpxpyp_implicit-------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+2), m(k+1), t(k+1), m(k), t(k), m(k-1), and
    ! t(k-1) correspond with altitudes zt(k+2), zm(k+1), zt(k+1), zm(k), zt(k),
    ! zm(k-1), and zt(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k+1) = 1 / ( zm(k+1) - zm(k) ); and
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  & ! gr%weights_zm2zt
        grid ! Type

    use constants_clubb, only: &
        zero    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Constant parameters
    integer, parameter :: &
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    integer, parameter :: &
      m_above = 1, & ! Index for upper momentum level grid weight.
      m_below = 2    ! Index for lower momentum level grid weight.

    ! ------------------------------ Input Variables ------------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      coef_wpxpyp_implicit,   & ! Coef. of <x'y'> in <w'x'y'>; t-levs  [m/s]
      rho_ds_zt,              & ! Dry, static density at t-levels      [kg/m^3]
      invrs_rho_ds_zm           ! Inv dry, static density at m-levels  [m^3/kg]

    logical, intent(in) :: &
      l_upwind_xpyp_turbulent_adv    ! Flag to use "upwind" discretization

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      sgn_turbulent_vel,       & ! Sign of the turbulent velocity      [-]
      coef_wpxpyp_implicit_zm, & ! coef_wpxpyp_implicit interp m-levs  [m/s]
      rho_ds_zm                  ! Dry, static density at m-levs       [kg/m^3]

    ! ------------------------------ Return Variable ------------------------------
    real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nz), intent(out) :: &
      lhs_ta    ! LHS coefficient of xpyp turbulent advection  [1/s]

    ! ------------------------------ Local Variables ------------------------------
    integer :: i, k, b    ! Vertical level index

    ! ------------------------------ Begin Code ------------------------------

    !$acc data copyin( gr, gr%weights_zm2zt, gr%invrs_dzm, gr%invrs_dzt, &
    !$acc              coef_wpxpyp_implicit, rho_ds_zt, invrs_rho_ds_zm, &
    !$acc              sgn_turbulent_vel, coef_wpxpyp_implicit_zm, rho_ds_zm ) &
    !$acc      copyout( lhs_ta )


    ! Set lower boundary array to 0
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do b = 1, ndiags3
        lhs_ta(b,i,1) = zero
      end do
    end do
    !$acc end parallel loop

    if ( .not. l_upwind_xpyp_turbulent_adv ) then

      ! Centered discretization.
      !$acc parallel loop gang vector collapse(2) default(present) 
      do k = 2, nz-1, 1
        do i = 1, ngrdcol
          
          ! Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
          lhs_ta(kp1_mdiag,i,k) &
          = invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
            * rho_ds_zt(i,k+1) * coef_wpxpyp_implicit(i,k+1) &
            * gr%weights_zm2zt(i,k+1,m_above)

          ! Momentum main diagonal: [ x xpyp(k,<t+1>) ]
          lhs_ta(k_mdiag,i,k) &
          = invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
            * ( rho_ds_zt(i,k+1) * coef_wpxpyp_implicit(i,k+1) &
                * gr%weights_zm2zt(i,k+1,m_below) &
                - rho_ds_zt(i,k) * coef_wpxpyp_implicit(i,k) &
                  * gr%weights_zm2zt(i,k,m_above) )

          ! Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
          lhs_ta(km1_mdiag,i,k) &
          = - invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
              * rho_ds_zt(i,k) * coef_wpxpyp_implicit(i,k) &
              * gr%weights_zm2zt(i,k,m_below)

        end do
      end do
      !$acc end parallel loop
      
    else ! l_upwind_xpyp_turbulent_adv

      ! "Upwind" discretization
      !$acc parallel loop gang vector collapse(2) default(present) 
      do k = 2, nz-1, 1
        do i = 1, ngrdcol
        
          if ( sgn_turbulent_vel(i,k) > zero ) then

            ! The "wind" is blowing upward.

            ! Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
            lhs_ta(kp1_mdiag,i,k) = zero

            ! Momentum main diagonal: [ x xpyp(k,<t+1>) ]
            lhs_ta(k_mdiag,i,k) &
             = invrs_rho_ds_zm(i,k) * gr%invrs_dzt(i,k) &
               * rho_ds_zm(i,k) * coef_wpxpyp_implicit_zm(i,k)

            ! Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
            lhs_ta(km1_mdiag,i,k) &
             = - invrs_rho_ds_zm(i,k) * gr%invrs_dzt(i,k) &
                 * rho_ds_zm(i,k-1) * coef_wpxpyp_implicit_zm(i,k-1)

          else ! sgn_turbulent_vel < 0

            ! The "wind" is blowing downward.

            ! Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
            lhs_ta(kp1_mdiag,i,k) &
             = invrs_rho_ds_zm(i,k) * gr%invrs_dzt(i,k+1) &
               * rho_ds_zm(i,k+1) * coef_wpxpyp_implicit_zm(i,k+1)

            ! Momentum main diagonal: [ x xpyp(k,<t+1>) ]
            lhs_ta(k_mdiag,i,k) &
             = - invrs_rho_ds_zm(i,k) * gr%invrs_dzt(i,k+1) &
                 * rho_ds_zm(i,k) * coef_wpxpyp_implicit_zm(i,k)

            ! Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
            lhs_ta(km1_mdiag,i,k) = zero

          end if ! sgn_turbulent_vel

        end do
      end do
      !$acc end parallel loop

    endif

    ! Set upper boundary array to 
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do b = 1, ndiags3
        lhs_ta(b,i,nz) = zero
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine xpyp_term_ta_pdf_lhs

  !=============================================================================================
  subroutine xpyp_term_ta_pdf_lhs_godunov( nz, ngrdcol, gr, & ! Intent(in)
                                                coef_wpxpyp_implicit, & ! Intent(in)
                                                invrs_rho_ds_zm, rho_ds_zm,  & ! Intent(in)
                                                lhs_ta )
  ! Intent(out)
  ! Description:
  !   This subroutine is a revised version of xpyp_term_ta_pdf_lhs_all. The
  !   revisions are maded to use the  Godunov-like upwind scheme for the
  !   vertical discretization of the turbulent advection term. This subroutine 
  !   returns an array of 3 dimensional arrays, one for every grid level not
  !   including
  !   boundary values.
  ! 
  ! Optional Arguements:
  !   The optional arguements can be used to override the default indices. 
  !   from_level - low index, default 2
  !   to level   - high index, default gr%nz-1
  ! 
  ! Notes:
  !   This subroutine exists for testing of Godunov-like upwind scheme. 
  !   THIS SUBROUTINE DOES NOT HANDLE BOUNDARY CONDITIONS AND SETS THEM TO 0
  !---------------------------------------------------------------------------------------------

    use grid_class, only:  & ! for gr%weights_zm2zt
    grid ! Type

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none
    
    integer, parameter :: &
      kp1_mdiag = 1, & ! Momentum superdiagonal index.
      k_mdiag   = 2, & ! Momentum main diagonal index.
      km1_mdiag = 3    ! Momentum subdiagonal index.

    !------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
        coef_wpxpyp_implicit,     & ! Coef. of <x'y'> in <w'x'y'>; t-lev [m/s]
        invrs_rho_ds_zm,          & ! Inv dry, static density @ m-level [m^3/kg]
        rho_ds_zm                   ! Dry, static density at m-lev [kg/m^3]

    !------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(ndiags3,ngrdcol,nz), intent(out) :: &
        lhs_ta

    !---------------- Local Variables -------------------
    integer :: &
        i, k, b          ! Loop variable for current grid level

    !---------------- Begin Code -------------------

    !$acc data copyin( gr, gr%invrs_dzm, &
    !$acc              coef_wpxpyp_implicit, invrs_rho_ds_zm, rho_ds_zm ) &
    !$acc      copyout( lhs_ta )

    ! Set lower boundary array to 0
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do b = 1, ndiags3
        lhs_ta(b,i,1) = 0.0_core_rknd
      end do
    end do
    !$acc end parallel loop

    ! Godunov-like upwind discretization
    !$acc parallel loop gang vector collapse(2) default(present) 
    do k = 2, nz-1
      do i = 1, ngrdcol
        
        ! Momentum superdiagonal: [ x xpyp(k+1,<t+1>) ]
        lhs_ta(kp1_mdiag,i,k) = invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
                                * rho_ds_zm(i,k+1) &
                                * min(0.0_core_rknd,coef_wpxpyp_implicit(i,k+1))

        ! Momentum main diagonal: [ x xpyp(k,<t+1>) ]
        lhs_ta(k_mdiag,i,k) = invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
                              * rho_ds_zm(i,k) &
                              * ( max(0.0_core_rknd,coef_wpxpyp_implicit(i,k+1)) - &
                                  min(0.0_core_rknd,coef_wpxpyp_implicit(i,k)) )

        ! Momentum subdiagonal: [ x xpyp(k-1,<t+1>) ]
        lhs_ta(km1_mdiag,i,k) = - invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
                                  * rho_ds_zm(i,k-1) &
                                  * max(0.0_core_rknd,coef_wpxpyp_implicit(i,k) )
      end do
    end do
    !$acc end parallel loop

    ! Set upper boundary array to 0
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do b = 1, ndiags3
        lhs_ta(b,i,nz) = 0.0_core_rknd
      end do
    end do
    !$acc end parallel loop

    !$acc end data

    return

  end subroutine xpyp_term_ta_pdf_lhs_godunov

  !=============================================================================
  subroutine xpyp_term_ta_pdf_rhs( nz, ngrdcol, gr, term_wpxpyp_explicit,  & ! In
                                        rho_ds_zt, rho_ds_zm,                   & ! In
                                        invrs_rho_ds_zm,                        & ! In
                                        l_upwind_xpyp_turbulent_adv,            & ! In
                                        sgn_turbulent_vel,                      & ! In
                                        term_wpxpyp_explicit_zm,                & ! In
                                        rhs_ta )                                  ! Out

    ! Description:
    ! Turbulent advection of <w'x'>, <x'^2>, and <x'y'>:  explicit portion of
    ! the code.
    !
    ! 1) <w'x'>
    !
    ! The d<w'x'>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'^2 x'> )/dz.
    !
    ! The value of <w'^2 x'> is found by integrating over the multivariate PDF
    ! of w and x, as detailed in function calc_wp2xp_pdf, which is found in
    ! module pdf_closure_module in pdf_closure_module.F90.
    !
    ! The equation obtained for <w'^2 x'> is written in terms of PDF parameters.
    ! Substitutions that are specific to the type of PDF used are made for the
    ! PDF parameters in order to write the <w'x'> turbulent advection term in
    ! the following form:
    !
    ! <w'^2 x'> = coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit.
    !
    ! For the ADG1 PDF, coef_wp2xp_implicit is a_1 * < w'^3 > / < w'^2 >, where
    ! a_1 is 1 / ( 1 - sigma_sqd_w ).  The value of term_wp2xp_explicit is 0, as
    ! the <w'x'> turbulent advection term is entirely implicit.
    !
    ! For the new PDF, the calculations of both coef_wp2xp_implicit and
    ! term_wp2xp_explicit are detailed in function calc_coefs_wp2xp_semiimpl,
    ! which is found in module new_pdf in new_pdf.F90.
    !
    ! For the new hybrid PDF, the calculation of coef_wp2xp_implicit is
    ! detailed in function calc_coef_wp2xp_implicit, which is found in module
    ! new_hybrid_pdf in new_hybrid_pdf.F90.  The value of term_wp2xp_explicit
    ! is 0, as the <w'x'> turbulent advection term is entirely implicit.
    !
    ! For explicit turbulent advection, the value of coef_wp2xp_implicit is 0
    ! and the value of term_wp2xp_explicit is <w'^2 x'>, as calculated by
    ! retaining the equation for <w'^2 x'> that is written in terms of PDF
    ! parameters.  This is a general form that can be used with any type of PDF.
    !
    ! The <w'x'> turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * ( coef_wp2xp_implicit * <w'x'> + term_wp2xp_explicit ) )
    !     /dz.
    !
    ! The variable <w'x'> is evaluated at the (t+1) timestep, which allows the
    ! <w'x'> turbulent advection term to be expressed semi-implicitly as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wp2xp_implicit * <w'x'>(t+1) )/dz
    ! - (1/rho_ds) * d( rho_ds * term_wp2xp_explicit )/dz.
    !
    ! The explicit portion of <w'x'> turbulent advection term is:
    !
    ! - (1/rho_ds) * d( rho_ds * term_wp2xp_explicit )/dz.
    !
    ! 2) <x'^2>
    !
    ! The d<x'^2>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'x'^2> )/dz;
    !
    ! The value of <w'x'^2> is found by integrating over the multivariate PDF of
    ! w and x, as detailed in function calc_wpxp2_pdf, which is found in module
    ! pdf_closure_module in pdf_closure_module.F90.
    !
    ! The equation obtained for <w'x'^2> is written in terms of PDF parameters.
    ! Substitutions that are specific to the type of PDF used are made for the
    ! PDF parameters in order to write the <x'^2> turbulent advection term in
    ! the following form:
    !
    ! <w'x'^2> = coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit.
    !
    ! For the ADG1 PDF, the value of coef_wpxp2_implicit is
    ! (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxp2_explicit
    ! is (1-(1/3)*beta) * a_1^2 * < w'x' >^2 * < w'^3 > / < w'^2 >^2, where
    ! a_1 is 1 / ( 1 - sigma_sqd_w ).
    !
    ! For the new PDF, the calculation of coef_wpxp2_implicit is detailed in
    ! function calc_coef_wpxp2_implicit, which is found in module new_pdf in
    ! new_pdf.F90.  The value of term_wpxp2_explicit is 0, as the <x'^2>
    ! turbulent advection term is entirely implicit.
    !
    ! For the new hybrid PDF, the calculation of both coef_wpxp2_implicit and
    ! term_wpxp2_explicit are detailed in subroutine calc_coefs_wpxp2_semiimpl,
    ! which is found in module new_hybrid_pdf in new_hybrid_pdf.F90.
    !
    ! For explicit turbulent advection, the value of coef_wpxp2_implicit is 0
    ! and the value of term_wpxp2_explicit is <w'x'^2>, as calculated by
    ! retaining the equation for <w'x'^2> that is written in terms of PDF
    ! parameters.  This is a general form that can be used with any type of PDF.
    !
    ! The <x'^2> turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * ( coef_wpxp2_implicit * <x'^2> + term_wpxp2_explicit ) )
    !     /dz;
    !
    ! The variable <x'^2> is evaluated at the (t+1) timestep, which allows the
    ! <x'^2> turbulent advection term to be expressed semi-implicitly as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wpxp2_implicit * <x'^2>(t+1) )/dz;
    ! - (1/rho_ds) * d( rho_ds * term_wpxp2_explicit )/dz.
    !
    ! The explicit portion of <x'^2> turbulent advection term is:
    !
    ! - (1/rho_ds) * d( rho_ds * term_wpxp2_explicit )/dz.
    !
    ! 3) <x'y'>
    !
    ! The d<x'y'>/dt equation contains a turbulent advection term:
    !
    ! - (1/rho_ds) * d( rho_ds * <w'x'y'> )/dz.
    !
    ! The value of <w'x'y'> is found by integrating over the multivariate PDF of
    ! w, x, and y, as detailed in function calc_wpxpyp_pdf, which is found in
    ! module pdf_closure_module in pdf_closure_module.F90.
    !
    ! The equation obtained for <w'x'y'> is written in terms of PDF parameters.
    ! Substitutions that are specific to the type of PDF used are made for the
    ! PDF parameters in order to write the <x'y'> turbulent advection term in
    ! the following form:
    !
    ! <w'x'y'> = coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit.
    !
    ! For the ADG1 PDF, the value of coef_wpxpyp_implicit is
    ! (1/3)*beta * a_1 * < w'^3 > / < w'^2 >.  The value of term_wpxpyp_explicit
    ! is (1-(1/3)*beta) * a_1^2 * < w'x' > * < w'y' > * < w'^3 > / < w'^2 >^2,
    ! where a_1 is 1 / ( 1 - sigma_sqd_w ).
    !
    ! For the new PDF, the calculation of both coef_wpxpyp_implicit and
    ! term_wpxpyp_explicit are detailed in function calc_coefs_wpxpyp_semiimpl,
    ! which is found in module new_pdf in new_pdf.F90.
    !
    ! For the new hybrid PDF, the calculation of both coef_wpxpyp_implicit
    ! and term_wpxpyp_explicit are detailed in subroutine
    ! calc_coefs_wpxpyp_semiimpl, which is found in module new_hybrid_pdf in
    ! new_hybrid_pdf.F90.
    !
    ! For explicit turbulent advection, the value of coef_wpxpyp_implicit is 0
    ! and the value of term_wpxpyp_explicit is <w'x'y'>, as calculated by
    ! retaining the equation for <w'x'y'> that is written in terms of PDF
    ! parameters.  This is a general form that can be used with any type of PDF.
    !
    ! The <x'y'> turbulent advection term is rewritten as:
    !
    ! - (1/rho_ds)
    !   * d( rho_ds * ( coef_wpxpyp_implicit * <x'y'> + term_wpxpyp_explicit ) )
    !     /dz.
    !
    ! The variable <x'y'> is evaluated at the (t+1) timestep, which allows the
    ! <x'y'> turbulent advection term to be expressed semi-implicitly as:
    !
    ! - (1/rho_ds) * d( rho_ds * coef_wpxpyp_implicit * <x'y'>(t+1) )/dz
    ! - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.
    !
    ! The explicit portion of <x'y'> turbulent advection term is:
    !
    ! - (1/rho_ds) * d( rho_ds * term_wpxpyp_explicit )/dz.
    !
    ! When x and y are the same variable, <x'y'> reduces to <x'^2> and <w'x'y'>
    ! reduces to <w'x'^2>.  Likewise, when y is set equal to w, <x'y'> becomes
    ! <w'x'> and <w'x'y'> reduces to <w'^2 x'>.  The discretization and the code
    ! used in this function will be written generally in terms of <x'y'> and
    ! term_wpxpyp_explicit, but also applies to <x'^2> and term_wpxp2_explicit,
    ! as well as to <w'x'> and term_wp2xp_explicit.
    !
    ! The explicit discretization of this term is as follows:
    !
    ! 1) Centered Discretization
    !
    ! The values of <x'y'> are found on the momentum levels, while the values of
    ! term_wpxpyp_explicit are found on the thermodynamic levels, which is where
    ! they were originally calculated by the PDF.  Additionally, the values of
    ! rho_ds_zt are found on the thermodynamic levels, and the values of
    ! invrs_rho_ds_zm are found on the momentum levels.  At the thermodynamic
    ! levels, the values of term_wpxpyp_explicit are multiplied by rho_ds_zt.
    ! Then, the derivative (d/dz) of that expression is taken over the central
    ! momentum level, where it is multiplied by -invrs_rho_ds_zm.  This yields
    ! the desired result.
    !
    ! ---------term_wpxpyp_explicitp1-------rho_ds_ztp1------------------ t(k+1)
    !
    ! ==xpyp==d( rho_ds_zt * term_wpxpyp_explicit )/dz==invrs_rho_ds_zm== m(k)
    !
    ! ---------term_wpxpyp_explicit---------rho_ds_zt-------------------- t(k)
    !
    ! The vertical indices t(k+1), m(k), and t(k) correspond with altitudes
    ! zt(k+1), zm(k), and zt(k), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzm(k) = 1 / ( zt(k+1) - zt(k) )
    !
    ! 2) "Upwind" Discretization.
    !
    ! The values of <x'y'> are found on the momentum levels.  The values of
    ! term_wpxpyp_explicit are originally calculated by the PDF on the
    ! thermodynamic levels.  They are interpolated to the intermediate momentum
    ! levels as term_wpxpyp_explicit_zm.  Additionally, the values of rho_ds_zm
    ! and the values of invrs_rho_ds_zm are found on the momentum levels.  The
    ! sign of the turbulent velocity is found on the central momentum level.  At
    ! the momentum levels, the values of term_wpxpyp_explicit_zm are multiplied
    ! by rho_ds_zm.  Then, the derivative (d/dz) of that expression is taken.
    ! When the sign of the turbulent velocity is positive, the "wind" is coming
    ! from below, and the derivative involves the central momentum level and the
    ! momentum level immediately below it.  When the sign of the turbulent
    ! velocity is negative, the "wind" is coming from above, and the derivative
    ! involves the central momentum level and the momenum level immediately
    ! above it.  After the derivative is taken, it is multiplied by
    ! -invrs_rho_ds_zm at the central momentum level.  This yields the desired
    ! result.
    !
    ! The turbulent velocity for <x'y'> is <w'x'y'> / <x'y'>, which has units of
    ! m/s.  The sign of the turbulent velocity is sgn( <w'x'y'> / <x'y'> ),
    ! where:
    !
    ! sgn(x) = | 1; when x >= 0
    !          | -1; when x < 0.
    !
    ! The sign of the turbulent velocity can also be rewritten as
    ! sgn( <w'x'y'> ) / sgn( <x'y'> ).  When a variance (<x'^2>) is being solved
    ! for, y = x, and sgn( <x'^2> ) is always 1.  The sign of the turbulent
    ! velocity reduces to simply sgn( <w'x'^2> ).
    !
    ! ---------term_wpxpyp_explicit-------------------------------------- t(k+2)
    !
    ! ========term_wpxpyp_explicit_zm(interp.)=======rho_ds_zm=========== m(k+1)
    !
    ! ---------term_wpxpyp_explicit-------------------------------------- t(k+1)
    !
    ! =xpyp===term_wpxpyp_explicit_zm(interp.)=rho_ds_zm=invrs_rho_ds_zm= m(k)
    !
    ! ---------term_wpxpyp_explicit-------------------------------------- t(k)
    !
    ! ========term_wpxpyp_explicit_zm(interp.)=======rho_ds_zm=========== m(k-1)
    !
    ! ---------term_wpxpyp_explicit-------------------------------------- t(k-1)
    !
    ! The vertical indices t(k+2), m(k+1), t(k+1), m(k), t(k), m(k-1), and
    ! t(k-1) correspond with altitudes zt(k+2), zm(k+1), zt(k+1), zm(k), zt(k),
    ! zm(k-1), and zt(k-1), respectively.  The letter "t" is used for
    ! thermodynamic levels and the letter "m" is used for momentum levels.
    !
    ! invrs_dzt(k+1) = 1 / ( zm(k+1) - zm(k) ); and
    ! invrs_dzt(k) = 1 / ( zm(k) - zm(k-1) ).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only:  &
        grid ! Type

    use constants_clubb, only: &
        zero    ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! ----------------------------- Input Variables -----------------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type (grid), target, intent(in) :: gr
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      term_wpxpyp_explicit, & ! RHS: <w'x'y'> eq; t-levs      [m/s(x un)(y un)]
      rho_ds_zt,            & ! Dry, static density at t-levs          [kg/m^3]
      invrs_rho_ds_zm         ! Inv dry, static density at m-levs      [m^3/kg]

    logical, intent(in) :: &
      l_upwind_xpyp_turbulent_adv    ! Flag to use "upwind" discretization

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      sgn_turbulent_vel,       & ! Sign of the turbulent velocity           [-]
      term_wpxpyp_explicit_zm, & ! term_wpxpyp_expl. zm       [m/s(x un)(y un)]
      rho_ds_zm                  ! Dry, static density at m-levs       [kg/m^3]

    ! ----------------------------- Return Variable -----------------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      rhs_ta    ! RHS portion of xpyp turbulent advection      [(x un)(y un)/s]

    ! ----------------------------- Local Variables -----------------------------
    integer :: i, k    ! Vertical level index

    ! ----------------------------- Begin Code -----------------------------

    !$acc data copyin( gr, gr%invrs_dzm, gr%invrs_dzt, &
    !$acc              term_wpxpyp_explicit, rho_ds_zt, invrs_rho_ds_zm, &
    !$acc              sgn_turbulent_vel, term_wpxpyp_explicit_zm, rho_ds_zm ) &
    !$acc      copyout( rhs_ta )


    ! Set lower boundary value to 0
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      rhs_ta(i,1) = zero
    end do
    !$acc end parallel loop

    if ( .not. l_upwind_xpyp_turbulent_adv ) then

      ! Centered discretization.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 2, nz-1, 1
        do i = 1, ngrdcol
          
          rhs_ta(i,k) &
          = - invrs_rho_ds_zm(i,k) &
              * gr%invrs_dzm(i,k) * ( rho_ds_zt(i,k+1) * term_wpxpyp_explicit(i,k+1) &
                                 - rho_ds_zt(i,k) * term_wpxpyp_explicit(i,k) )
        end do
      end do ! k = 2, nz-1, 1
      !$acc end parallel loop

    else ! l_upwind_xpyp_turbulent_adv

      ! "Upwind" discretization
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 2, nz-1, 1
        do i = 1, ngrdcol

          if ( sgn_turbulent_vel(i,k) > zero ) then

             ! The "wind" is blowing upward.
             rhs_ta(i,k) &
             = - invrs_rho_ds_zm(i,k) &
                 * gr%invrs_dzt(i,k) &
                 * ( rho_ds_zm(i,k) * term_wpxpyp_explicit_zm(i,k) &
                     - rho_ds_zm(i,k-1) * term_wpxpyp_explicit_zm(i,k-1) )

          else ! sgn_turbulent_vel < 0

             ! The "wind" is blowing downward.
             rhs_ta(i,k) &
             = - invrs_rho_ds_zm(i,k) &
                 * gr%invrs_dzt(i,k+1) &
                 * ( rho_ds_zm(i,k+1) * term_wpxpyp_explicit_zm(i,k+1) &
                     - rho_ds_zm(i,k) * term_wpxpyp_explicit_zm(i,k) )

          endif ! sgn_turbulent_vel
          
        end do
      end do ! k = 2, nz-1, 1
      !$acc end parallel loop
      
    end if

    ! Set upper boundary value to 0
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      rhs_ta(i,nz) = zero
    end do
    !$acc end parallel loop

    !$acc end data 

    return

  end subroutine xpyp_term_ta_pdf_rhs

  !=============================================================================
  subroutine xpyp_term_ta_pdf_rhs_godunov( nz, ngrdcol, gr, & ! Intent(in)
                                                term_wpxpyp_explicit_zm, & ! Intent(in)
                                                invrs_rho_ds_zm, & ! Intent(in)
                                                sgn_turbulent_vel, & ! Intent(in)
                                                rho_ds_zm, & ! Intent(in)
                                                rhs_ta ) ! Intent(out)
    ! Description:
    !   This subroutine intends to add godunov upwind difference scheme based
    !   on xpyp_term_ta_pdf_rhs.  The revisions are maded to use the Godunov-like 
    !   upwind scheme for the vertical discretization. 
    !   This subroutine returns an array of values for every grid level.
    !
    ! Optional Arguements:
    !   The optional arguements can be used to override the default indices. 
    !   from_level - low index, default 2
    !   to level - high index, default gr%nz-1
    ! 
    ! Notes:
    !   This subroutine exists for testing of Godunov-like upwind scheme. 
    !   THIS SUBROUTINE DOES NOT HANDLE BOUNDARY CONDITIONS AND SETS THEM TO 0
    !--------------------------------------------------------------------------------------------
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: &
      grid ! Type

    implicit none

    !------------------- Input Variables -------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    type (grid), target, intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      term_wpxpyp_explicit_zm,      & ! RHS: <w'x'y'> eq; m-lev(k)   [m/s(x un)(y un)]
      invrs_rho_ds_zm,              & ! Inv dry, static density at m-lev (k) [m^3/kg]
      sgn_turbulent_vel,            & ! Sign of the turbulent velocity [-]
      rho_ds_zm                       ! Dry, static density at m-lev (k) [kg/m^3]

    !------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      rhs_ta

    !---------------- Local Variables -------------------
    integer :: &
      i, k             ! Loop variable for current grid level

    !---------------- Begin Code -------------------

    !$acc data copyin( gr, gr%invrs_dzm, &
    !$acc              term_wpxpyp_explicit_zm, invrs_rho_ds_zm, &
    !$acc              sgn_turbulent_vel, rho_ds_zm ) &
    !$acc      copyout( rhs_ta )

    ! Set lower boundary value to 0
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      rhs_ta(i,1) = 0.0_core_rknd
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 2, nz-1
      do i = 1, ngrdcol 
        rhs_ta(i,k) = - invrs_rho_ds_zm(i,k) * gr%invrs_dzm(i,k) &
                    * ( min(0.0_core_rknd, sgn_turbulent_vel(i,k+1)) &
                          * rho_ds_zm(i,k+1) * term_wpxpyp_explicit_zm(i,k+1) &
                      + max(0.0_core_rknd, sgn_turbulent_vel(i,k+1)) &   
                          * rho_ds_zm(i,k)   * term_wpxpyp_explicit_zm(i,k) &
                      - min(0.0_core_rknd, sgn_turbulent_vel(i,k))&
                          * rho_ds_zm(i,k)   * term_wpxpyp_explicit_zm(i,k) &
                      - max(0.0_core_rknd, sgn_turbulent_vel(i,k)) &
                          * rho_ds_zm(i,k-1) * term_wpxpyp_explicit_zm(i,k-1) &
                      )
      end do
    end do
    !$acc end parallel loop

    ! Set upper boundary value to 0
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      rhs_ta(i,nz) = 0.0_core_rknd
    end do
    !$acc end parallel loop

    !$acc end data 

    return

  end subroutine xpyp_term_ta_pdf_rhs_godunov

!===============================================================================

end module turbulent_adv_pdf
