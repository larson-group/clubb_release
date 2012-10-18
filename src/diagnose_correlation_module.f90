! $Id$
module diagnose_correlations_module 

  use clubb_precision, only: &
      core_rknd
  
  use grid_class, only: &
      gr

  implicit none 

  public :: diagnose_KK_corr, diagnose_LH_corr, &
            calc_mean, calc_varnce
            

  private :: diagnose_corr, calc_w_corr

  real( kind = core_rknd ) :: &
      spwp,  &  ! Covariance of s and w on the zt-grid
      rrpwp, &  ! Covariance of rrain and w on the zt-grid
      Nrpwp, &  ! Covariance of Nr and w on the zt-grid
      Ncpwp, &  ! Covariance of Nc and w on the zt-grid
      stdev_w

  contains 



!-----------------------------------------------------------------------
  subroutine diagnose_KK_corr( Ncm, rrainm, Nrm, & ! intent(in)
                               Ncp2_on_Ncm2, rrp2_on_rrm2, Nrp2_on_Nrm2, &
                               corr_sw, corr_rrw, corr_Nrw, corr_Ncw, &
                               pdf_params, &
                               corr_rrNr, corr_srr, corr_sNr, corr_sNc ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it 
    !   into KK microphysics.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB-Trac:ticket:514)
    !-----------------------------------------------------------------------


    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use grid_class, only: &
        gr ! Variable(s)

    use variables_diagnostic_module, only: &
        Kh_zt

    use model_flags, only: &
        l_calc_w_corr

    use pdf_parameter_module, only: &
        pdf_parameter  ! Type

    use constants_clubb, only: &
        w_tol,         & ! [m/s]
        s_mellor_tol,  & ! [kg/kg]
        Nc_tol,        & ! [#/kg]
        rr_tol,        & ! [kg/kg] 
        Nr_tol           ! [#/kg]

    implicit none

    intrinsic :: sqrt

    ! Local Constants
    integer, parameter :: &
      n_variables = 5

    ! Input Variables

    real( kind = core_rknd ), intent(in) :: &
      Ncm,            &  
      rrainm,         &
      Nrm,            &
      Ncp2_on_Ncm2,   &
      rrp2_on_rrm2,   &
      Nrp2_on_Nrm2,   &
      corr_sw,        &  ! Correlation between s_mellor and w
      corr_rrw,       &  ! Correlation between rrain and w
      corr_Nrw,       &  ! Correlation between Nr and w
      corr_Ncw           ! Correlation between Nc and w
      
    type(pdf_parameter), target, intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout) :: &
      corr_rrNr, &
      corr_srr,  &
      corr_sNr,  &
      corr_sNc


    ! Local Variables
    real( kind = core_rknd ), dimension(n_variables, n_variables) :: &
      corr_matrix_approx

    real( kind = core_rknd ), dimension(n_variables) :: &
      sqrt_xp2_on_xm2, & ! sqrt of x_variance / x_mean^2
      xm,              & ! means of the hydrometeors
      stdev_x,         & ! ratios of x_variance / x_mean^2
      wpxp,            & ! Covariances of w with the hydrometeors
      xp2_on_xm2,      & ! ratios of x_variance over x_mean^2
      x_tol              ! tolerances for the x variables

    real( kind = core_rknd ) :: &
      s_mellorp2_on_s_mellorm2

    ! Indices of the hydrometeors
    integer :: &
      ii_w = 1, &
      ii_s = 2, &
      ii_rrain = 3, &
      ii_Nr = 4, &
      ii_Nc = 5

    integer :: i, j ! Loop Iterators


    !-------------------- Begin code --------------------

    ! S_i is set to 1 for s_mellor, because s_mellorm could be 0
    s_mellorp2_on_s_mellorm2 = 1

    ! set up xp2_on_xm2
    sqrt_xp2_on_xm2(ii_w) = 1
    sqrt_xp2_on_xm2(ii_s) = 1
    sqrt_xp2_on_xm2(ii_rrain) = sqrt(rrp2_on_rrm2)
    sqrt_xp2_on_xm2(ii_Nr) = sqrt(Nrp2_on_Nrm2)
    sqrt_xp2_on_xm2(ii_Nc) = sqrt(Ncp2_on_Ncm2)

    ! initialize the correlation matrix with 0
    do i=2, n_variables
       do j=1, n_variables
          corr_matrix_approx(i,j) = 0
       end do
    end do

    ! set diagonal of the corraltion matrix to 1
    do i = 1, n_variables
       corr_matrix_approx(i,i) = 1
    end do


    if ( l_calc_w_corr ) then

      ! set up xm
      xm(ii_w) = 1
      xm(ii_s) = calc_mean( pdf_params%mixt_frac, pdf_params%s1, pdf_params%s2)
      xm(ii_rrain) = rrainm
      xm(ii_Nr) = Nrm
      xm(ii_Nc) = Ncm

      ! calculate the standard deviations
      stdev_x(ii_w) = stdev_w

      ! same formula as wpsp, but spsp=sp2 and therefore sqrt(sp2)=stdev_s
      stdev_x(ii_s) = sqrt( pdf_params%mixt_frac * ( 1 - pdf_params%mixt_frac ) ) &
                      * ( pdf_params%s1 - pdf_params%s2 ) 

      stdev_x(ii_rrain) = sqrt_xp2_on_xm2(ii_rrain) * xm(ii_rrain)
      stdev_x(ii_Nr) = sqrt_xp2_on_xm2(ii_Nr) * xm(ii_Nr)
      stdev_x(ii_Nc) = sqrt_xp2_on_xm2(ii_Nc) * xm(ii_Nc)

      ! set up wpxp
      wpxp(ii_w) = 1
      wpxp(ii_s) = spwp
      wpxp(ii_rrain) = rrpwp
      wpxp(ii_Nr) = Nrpwp
      wpxp(ii_Nc) = Ncpwp

      ! set up x_tol
      x_tol(ii_w) = w_tol
      x_tol(ii_s) = s_mellor_tol
      x_tol(ii_rrain) = rr_tol
      x_tol(ii_Nr) = Nr_tol
      x_tol(ii_Nc) = Nc_tol

      call calc_w_corr( n_variables, wpxp, stdev_x, x_tol, & ! intent(in) 
                          corr_matrix_approx )  ! intent(inout)


    else

      ! set the first row to the corresponding prescribed correlations
      corr_matrix_approx(1,ii_s) = corr_sw
      corr_matrix_approx(1,ii_rrain) = corr_rrw
      corr_matrix_approx(1,ii_Nr) = corr_Nrw
      corr_matrix_approx(1,ii_Nc) = corr_Ncw

    end if ! l_calc_w_corr

    call diagnose_corr( n_variables, sqrt_xp2_on_xm2, & !intent(in)
                        corr_matrix_approx ) ! intent(inout)    

    corr_rrNr = corr_matrix_approx(ii_rrain, ii_Nr)
    corr_srr = corr_matrix_approx(ii_s, ii_rrain)
    corr_sNr = corr_matrix_approx(ii_s, ii_Nr)
    corr_sNc = corr_matrix_approx(ii_s, ii_Nc)

  end subroutine diagnose_KK_corr

!-----------------------------------------------------------------------
  subroutine diagnose_LH_corr( xp2_on_xm2, & !intent(in)
                               corr_matrix_approx ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it 
    !   into SILHS microphysics.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    ! Local Constants
    integer, parameter :: &
      n_variables = 5

    ! Input/Output variables
    real( kind = core_rknd ), dimension(n_variables, n_variables), intent(inout) :: &
      corr_matrix_approx

    ! Local Variables
    real( kind = core_rknd ), dimension(n_variables) :: & 
      xp2_on_xm2 ! ratios of x_variance over x_mean^2

  
  end subroutine diagnose_LH_corr

!-----------------------------------------------------------------------
  subroutine diagnose_corr( n_variables, sqrt_xp2_on_xm2, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix for each timestep.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use error_code, only:  & 
        clubb_at_least_debug_level ! Function

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use parameters_tunable, only:  & 
        alpha_corr ! Constant(s)

    use constants_clubb, only: &
      max_mag_correlation

    implicit none

    intrinsic :: &
      sqrt, abs, sign

    ! Input Variables
    integer :: &
      n_variables  ! number of variables in the correlation matrix
    
    real( kind = core_rknd ), dimension(n_variables), intent(in) :: & 
      sqrt_xp2_on_xm2    ! sqrt of x_variance / x_mean^2

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(inout) :: &
      corr_matrix_approx ! correlation matrix


    ! Local Variables
    integer :: i, j ! Loop iterator
    real( kind = core_rknd ) :: f_ij
    real( kind = core_rknd ), dimension(n_variables) :: &
      s_1j ! s_1j = sqrt(1-c_1j^2)


    !-------------------- Begin code --------------------

    ! calculate all square roots
    do i = 1, n_variables

       s_1j(i) = sqrt(1._core_rknd-corr_matrix_approx(1,i)**2)

    end do


    ! Diagnose the missing correlations (upper triangle)
    do i = 2, (n_variables-1)
      do j = (i+1), n_variables

        ! formula (16) in the ref. paper 
        f_ij = alpha_corr * sqrt_xp2_on_xm2(i) * sqrt_xp2_on_xm2(j) &
        * sign(1.0_core_rknd,corr_matrix_approx(1,i)*corr_matrix_approx(1,j)) 

        ! make sure -1 < f_ij < 1
        if ( clubb_at_least_debug_level( 2 ) ) then

           if ( f_ij < -max_mag_correlation ) then

              f_ij = -max_mag_correlation

           else if ( f_ij > max_mag_correlation ) then

              f_ij = max_mag_correlation

           end if

        end if

        ! formula (15) in the ref. paper
        corr_matrix_approx(i,j) = corr_matrix_approx(1,i) * corr_matrix_approx(1,j) &
        + f_ij * s_1j(i) * s_1j(j)

      end do ! do j
    end do ! do i
    
  end subroutine diagnose_corr 

  !------------------------------------------------------------------------
  subroutine setup_w_covars( spwp_zt, rrpwp_zt, Nrpwp_zt, Ncpwp_zt, std_dev_w )  ! intent(in)
    ! Description:
    ! Setup the covariances of w with the hydrometeors. The covariances at 
    ! the current height level are stored in global variales. These will be 
    ! needed to calculate the corresponding correlations. 

    ! References:
    ! clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    implicit none

    ! Input Variables    
    real( kind = core_rknd ), intent(in) ::  &
      spwp_zt,  &  ! Covariance of s and w on the zt-grid
      rrpwp_zt, &  ! Covariance of rrain and w on the zt-grid
      Nrpwp_zt, &  ! Covariance of Nr and w on the zt-grid
      Ncpwp_zt, &  ! Covariance of Nc and w on the zt-grid
      std_dev_w

    ! --- Begin Code ---

    spwp = spwp_zt
    rrpwp = rrpwp_zt
    Nrpwp = Nrpwp_zt
    Ncpwp = Ncpwp_zt
    stdev_w = std_dev_w

  end subroutine setup_w_covars

  !-----------------------------------------------------------------------
  subroutine calc_w_corr( n_variables, wpxp, stdev_x, x_tol, & ! intent(in)
                          corr_matrix_approx )  ! intent(inout)
    ! Description:
    ! Compute the correlations of w with the hydrometeors.

    ! References:
    ! clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: &
      max_mag_correlation

    implicit none

    intrinsic :: max

    ! Input Variables
    integer, intent(in) :: &
      n_variables  ! Number of variables involved (size of xm et al.)

    real( kind = core_rknd ), dimension(n_variables), intent(in) :: &
      stdev_x,  & ! standard deviation of x ( x(1)=sqrt(wp2) )
      wpxp,     & ! Covariances of w with the hydrometeors ( wpxp(1)=1 to keep the vector length equal )
      x_tol       ! tolerances for the x variables
    
    ! Input/Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(inout) :: &
      corr_matrix_approx ! correlation matrix

    ! Local Variables
    integer :: i ! Loop iterator

    ! --- Begin Code ---
  
    do i = 2, n_variables

      corr_matrix_approx(1, i) = wpxp(i) / ( max(stdev_x(i), x_tol(i)) * max(stdev_x(1), x_tol(1)) )

      ! Make sure the correlation is in [-1,1]
      if ( corr_matrix_approx(1, i) < -max_mag_correlation ) then

        corr_matrix_approx(1, i) = -max_mag_correlation

      else if ( corr_matrix_approx(1, i) > max_mag_correlation ) then

        corr_matrix_approx(1, i) = max_mag_correlation

      end if

    end do
    
  end subroutine calc_w_corr

  !-----------------------------------------------------------------------
  function calc_varnce( mixt_frac, x1, x2, xm, x1p2, x2p2 )

    ! Description:
    ! Calculate the variance xp2 from the components x1, x2.

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
    !   page 3535
    !-----------------------------------------------------------------------
    
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mixt_frac, &
      x1, &         ! first component of the double gaussian
      x2, &         ! second component of the double gaussian
      xm, &         ! mean of x
      x1p2, &       ! variance of the first component
      x2p2          ! variance of the second component

    ! Return Variable
    real( kind = core_rknd ) :: &
      calc_varnce

    ! --- Begin Code ---

    calc_varnce = mixt_frac * ((x1 - xm)**2 + x1p2) + (1 - mixt_frac) * ((x2 - xm)**2 + x2p2)

    return
  end function calc_varnce

  !-----------------------------------------------------------------------
  function calc_mean( mixt_frac, x1, x2 )

    ! Description:
    ! Calculate the mean xm from the components x1, x2.

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02,
    !   page 3535
    !-----------------------------------------------------------------------
    
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ) :: &
      mixt_frac, &
      x1, &         ! first component of the double gaussian
      x2            ! second component of the double gaussian

    ! Return Variable
    real( kind = core_rknd ) :: &
      calc_mean

    ! --- Begin Code ---

    calc_mean = mixt_frac * x1 + (1 - mixt_frac) * x2

    return
  end function calc_mean

end module diagnose_correlations_module 
