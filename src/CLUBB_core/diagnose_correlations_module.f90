! $Id$
module diagnose_correlations_module 

  use clubb_precision, only: &
      core_rknd

  implicit none 

  public :: diagnose_KK_corr, diagnose_LH_corr, &
            calc_mean, calc_varnce, calc_w_corr, &
            corr_stat_output
            

  private :: diagnose_corr 


  contains 

!-----------------------------------------------------------------------
  subroutine diagnose_KK_corr( Ncnm, rrainm, Nrm, & ! intent(in)
                               Ncnp2_on_Ncnm2, rrp2_on_rrm2, Nrp2_on_Nrm2, &
                               corr_ws, corr_wrr, corr_wNr, corr_wNcn, &
                               pdf_params, &
                               corr_rrNr_p, corr_srr_p, corr_sNr_p, corr_sNcn_p, &
                               corr_rrNr, corr_srr, corr_sNr, corr_sNcn ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it 
    !   into KK microphysics.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB-Trac:ticket:514)
    !-----------------------------------------------------------------------


    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Type

    use constants_clubb, only: &
        w_tol,         & ! [m/s]
        s_mellor_tol,  & ! [kg/kg]
        Ncn_tol,       & ! [num/kg]
        rr_tol,        & ! [kg/kg] 
        Nr_tol           ! [num/kg]

    use stats_type, only: & 
        stat_update_var_pt  ! Procedure(s)

    implicit none

    intrinsic :: sqrt

    ! Local Constants
    integer, parameter :: &
      n_variables = 5

    ! Input Variables

    real( kind = core_rknd ), intent(in) :: &
      Ncnm,           &  ! Mean cloud nuclei concentration             [num/kg]
      rrainm,         &  ! Mean rain water mixing ratio (overall)      [kg/kg]
      Nrm,            &  ! Mean rain drop concentration (overall)      [num/kg]
      Ncnp2_on_Ncnm2, &  ! Variance of Ncn divided by Ncnm^2           [-]
      rrp2_on_rrm2,   &  ! Variance of rr divided by rrm^2             [-]
      Nrp2_on_Nrm2,   &  ! Variance of Nr divided by Nrm^2             [-]
      corr_ws,        &  ! Correlation between s_mellor and w          [-]
      corr_wrr,       &  ! Correlation between rrain and w             [-]
      corr_wNr,       &  ! Correlation between Nr and w                [-]
      corr_wNcn,      &  ! Correlation between Ncn and w               [-]
      corr_rrNr_p,    &  ! Prescribed correlation between rr and Nr    [-]
      corr_srr_p,     &  ! Prescribed correlation between s and rr     [-]
      corr_sNr_p,     &  ! Prescribed correlation between s and Nr     [-]
      corr_sNcn_p        ! Prescribed correlation between s and Ncn    [-]
      
    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters  [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout) :: &
      corr_rrNr,   &  ! Correlation between rr and Nr    [-]
      corr_srr,    &  ! Correlation between s and rr     [-]
      corr_sNr,    &  ! Correlation between s and Nr     [-]
      corr_sNcn       ! Correlation between s and Ncn    [-]



    ! Local Variables
    real( kind = core_rknd ), dimension(n_variables, n_variables) :: &
      corr_matrix_approx, & ! [-]
      corr_matrix_prescribed ! [-]

    real( kind = core_rknd ), dimension(n_variables) :: &
      sqrt_xp2_on_xm2, & ! sqrt of x_variance / x_mean^2          [units vary]
      xm                 ! means of the hydrometeors              [units vary]

    ! Indices of the hydrometeors
    integer :: &
      ii_w     = 1, &
      ii_s     = 2, &
      ii_rrain = 3, &
      ii_Nr    = 4, &
      ii_Ncn   = 5

    integer :: i, j ! Loop Iterators


    !-------------------- Begin code --------------------

    ! Remove compiler warnings about unused variables.
    if ( .false. ) then
       xm(ii_rrain) = rrainm
       xm(ii_Nr)    = Nrm
       xm(ii_Ncn)   = Ncnm
       print *, "pdf_params = ", pdf_params
    endif

    ! set up xp2_on_xm2

    ! TODO Why is wp2_on_wm2=1
    ! S_i is set to 1 for s_mellor and w, because s_mellorm could be 0
    sqrt_xp2_on_xm2(ii_w) = 1._core_rknd
    sqrt_xp2_on_xm2(ii_s) = 1._core_rknd

    sqrt_xp2_on_xm2(ii_rrain) = sqrt(rrp2_on_rrm2)
    sqrt_xp2_on_xm2(ii_Nr)    = sqrt(Nrp2_on_Nrm2)
    sqrt_xp2_on_xm2(ii_Ncn)   = sqrt(Ncnp2_on_Ncnm2)

    ! initialize the correlation matrix with 0
    do i=1, n_variables
       do j=1, n_variables
          corr_matrix_approx(i,j) = 0._core_rknd
          corr_matrix_prescribed(i,j) = 0._core_rknd
       end do
    end do

    ! set diagonal of the correlation matrix to 1
    do i = 1, n_variables
       corr_matrix_approx(i,i) = 1._core_rknd
       corr_matrix_prescribed(i,i) = 1._core_rknd
    end do


    ! set the first row to the corresponding prescribed correlations
    corr_matrix_approx(ii_s,1)     = corr_ws
    corr_matrix_approx(ii_rrain,1) = corr_wrr
    corr_matrix_approx(ii_Nr,1)    = corr_wNr
    corr_matrix_approx(ii_Ncn,1)   = corr_wNcn

    !corr_matrix_prescribed = corr_matrix_approx

    ! set up the prescribed correlation matrix
    if( ii_rrain > ii_Nr ) then
      corr_matrix_prescribed(ii_rrain, ii_Nr) = corr_rrNr_p
    else
      corr_matrix_prescribed(ii_Nr, ii_rrain) = corr_rrNr_p
    end if

    if ( ii_s > ii_rrain ) then
      corr_matrix_prescribed(ii_s, ii_rrain) = corr_srr_p
    else
      corr_matrix_prescribed(ii_rrain, ii_s) = corr_srr_p
    end if

    if ( ii_s > ii_Nr ) then
      corr_matrix_prescribed(ii_s, ii_Nr) = corr_sNr_p
    else
      corr_matrix_prescribed(ii_Nr, ii_s) = corr_sNr_p
    end if

    if ( ii_s > ii_Ncn ) then
      corr_matrix_prescribed(ii_s, ii_Ncn) = corr_sNcn_p
    else
      corr_matrix_prescribed(ii_Ncn, ii_s) = corr_sNcn_p
    end if
    
    call diagnose_corr( n_variables, sqrt_xp2_on_xm2, corr_matrix_prescribed, & !intent(in)
                        corr_matrix_approx ) ! intent(inout)    

    if( ii_rrain > ii_Nr ) then
      corr_rrNr = corr_matrix_approx(ii_rrain, ii_Nr)
    else
      corr_rrNr = corr_matrix_approx(ii_Nr, ii_rrain)
    end if

    if ( ii_s > ii_rrain ) then
      corr_srr = corr_matrix_approx(ii_s, ii_rrain)
    else
      corr_srr = corr_matrix_approx(ii_rrain, ii_s)
    end if

    if ( ii_s > ii_Nr ) then
      corr_sNr = corr_matrix_approx(ii_s, ii_Nr)
    else
      corr_sNr = corr_matrix_approx(ii_Nr, ii_s)
    end if

    if ( ii_s > ii_Ncn ) then
      corr_sNcn = corr_matrix_approx(ii_s, ii_Ncn)
    else
      corr_sNcn = corr_matrix_approx(ii_Ncn, ii_s)
    end if

  end subroutine diagnose_KK_corr

!-----------------------------------------------------------------------
  subroutine diagnose_LH_corr( xp2_on_xm2, d_variables, corr_matrix_prescribed, & !intent(in)
                               corr_array ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it 
    !   into SILHS microphysics.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_matrix_module, only: &
      iiLH_w ! Variable(s)

    implicit none

    intrinsic :: max, sqrt, transpose

    ! Input Variables
    integer, intent(in) :: d_variables

    real( kind = core_rknd ), dimension(d_variables, d_variables), intent(in) :: &
      corr_matrix_prescribed

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      xp2_on_xm2 ! ratios of x_variance over x_mean^2

    ! Input/Output variables
    real( kind = core_rknd ), dimension(d_variables, d_variables), intent(inout) :: &
      corr_array

    ! Local Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables) :: &
      corr_matrix_pre_swapped

    real( kind = core_rknd ), dimension(d_variables) :: &
      swap_array

    !-------------------- Begin code --------------------

    ! Swap the w-correlations to the first row
    swap_array = corr_array(:, 1)
    corr_array(1:iiLH_w, 1) = corr_array(iiLH_w, iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, 1) = corr_array((iiLH_w+1):d_variables, iiLH_w)
    corr_array(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

    corr_matrix_pre_swapped = corr_matrix_prescribed
    swap_array = corr_matrix_pre_swapped (:,1)
    corr_matrix_pre_swapped(1:iiLH_w, 1) = corr_matrix_pre_swapped(iiLH_w, iiLH_w:1:-1)
    corr_matrix_pre_swapped((iiLH_w+1):d_variables, 1) = corr_matrix_pre_swapped( &
                                                         (iiLH_w+1):d_variables, iiLH_w)
    corr_matrix_pre_swapped(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_matrix_pre_swapped((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

    ! diagnose correlations
    call diagnose_corr( d_variables, sqrt(xp2_on_xm2), corr_matrix_pre_swapped, &
                        corr_array)

    ! Swap rows back
    swap_array = corr_array(:, 1)
    corr_array(1:iiLH_w, 1) = corr_array(iiLH_w, iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, 1) = corr_array((iiLH_w+1):d_variables, iiLH_w)
    corr_array(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

  end subroutine diagnose_LH_corr

!-----------------------------------------------------------------------
  subroutine diagnose_correlations( d_variables, corr_array_pre, & ! Intent(in)
                                    corr_ws, corr_wrr, &           ! Intent(out)
                                    corr_wNr, corr_wNcn, &         ! Intent(out)
                                    corr_st, corr_srr, &           ! Intent(out)
                                    corr_sNr, corr_sNcn, &         ! Intent(out)
                                    corr_trr, corr_tNr, &          ! Intent(out)
                                    corr_tNcn, corr_rrNr )         ! Intent(out)
    ! Description:
    !   This subroutine diagnoses the correlation matrix in order to feed it
    !   into SILHS microphysics.

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_matrix_module, only: &
        iiLH_w ! Variable(s)

    use constants_clubb, only: &
        zero

    use model_flags, only: &
        l_calc_w_corr ! Flag(s)

    implicit none

    intrinsic :: max, sqrt, transpose

    ! Input Variables
    integer, intent(in) :: &
      d_variables   ! number of diagnosed correlations

    real( kind = core_rknd ), dimension(d_variables, d_variables), intent(in) :: &
      corr_array_pre   ! Prescribed correlations

    ! Output variables
    real( kind = core_rknd ), intent(out) :: &
      corr_ws,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_wrr,    & ! Correlation between w and rr (1st PDF component) ip [-]
      corr_wNr,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      corr_wNcn,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      corr_st,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_srr,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_sNr,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNcn,   & ! Correlation between s and Ncn (1st PDF component)   [-]
      corr_trr,    & ! Correlation between t and rr (1st PDF component) ip [-]
      corr_tNr,    & ! Correlation between t and Nr (1st PDF component) ip [-]
      corr_tNcn,   & ! Correlation between t and Ncn (1st PDF component)   [-]
      corr_rrNr      ! Correlation between rr & Nr (1st PDF component) ip  [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables) :: &
      corr_array_pre_swapped
      
    real( kind = core_rknd ), dimension(d_variables) :: &
      swap_array

    real( kind = core_rknd ), dimension(d_variables, d_variables) :: &
      corr_array

    ! We actually don't need this right now
    real( kind = core_rknd ), dimension(d_variables) :: &
      xp2_on_xm2_array   ! ratios of x_variance over x_mean^2

    integer :: i ! Loop iterator

    !-------------------- Begin code --------------------

    ! Initialize xp2_on_xm2_array
    do i = 1, d_variables
       xp2_on_xm2_array(i) = zero
    end do 

    ! Swap the w-correlations to the first row for the prescribed correlations
    corr_array_pre_swapped = corr_array_pre
    swap_array = corr_array_pre_swapped (:,1)
    corr_array_pre_swapped(1:iiLH_w, 1) = corr_array_pre_swapped(iiLH_w, iiLH_w:1:-1)
    corr_array_pre_swapped((iiLH_w+1):d_variables, 1) = corr_array_pre_swapped( &
                                                         (iiLH_w+1):d_variables, iiLH_w)
    corr_array_pre_swapped(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_array_pre_swapped((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

    ! diagnose correlations
    
    if ( .not. l_calc_w_corr ) then
       corr_array = corr_array_pre_swapped
    endif

    call diagnose_corr( d_variables, sqrt(xp2_on_xm2_array), corr_array_pre_swapped, &
                        corr_array )

    ! Swap rows back
    swap_array = corr_array(:, 1)
    corr_array(1:iiLH_w, 1) = corr_array(iiLH_w, iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, 1) = corr_array((iiLH_w+1):d_variables, iiLH_w)
    corr_array(iiLH_w, 1:iiLH_w) = swap_array(iiLH_w:1:-1)
    corr_array((iiLH_w+1):d_variables, iiLH_w) = swap_array((iiLH_w+1):d_variables)

    ! Unpack the corr_array to the corresponding variables
    call unpack_correlations( d_variables, corr_array, & ! Intent(in)
                              corr_ws, corr_wrr, corr_wNr, corr_wNcn, & ! Intent(out)
                              corr_st, corr_srr, corr_sNr, corr_sNcn, &
                              corr_trr, corr_tNr, corr_tNcn, corr_rrNr )  


  end subroutine diagnose_correlations

!-----------------------------------------------------------------------
  subroutine diagnose_corr( n_variables, sqrt_xp2_on_xm2, corr_matrix_prescribed, & !intent(in)
                            corr_matrix_approx ) ! intent(inout)

    ! Description:
    !   This subroutine diagnoses the correlation matrix for each timestep.   

    ! References:
    !   Larson et al. (2011), J. of Geophysical Research, Vol. 116, D00T02
    !   (see CLUBB Trac ticket#514)
    !-----------------------------------------------------------------------

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
    integer, intent(in) :: &
      n_variables  ! number of variables in the correlation matrix [-]
    
    real( kind = core_rknd ), dimension(n_variables), intent(in) :: & 
      sqrt_xp2_on_xm2    ! sqrt of x_variance / x_mean^2 [units vary]

    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(in) :: &
      corr_matrix_prescribed ! correlation matrix [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(n_variables,n_variables), intent(inout) :: &
      corr_matrix_approx ! correlation matrix [-]


    ! Local Variables
    integer :: i, j ! Loop iterator

    real( kind = core_rknd ) :: &
      f_ij, &
      f_ij_o

    real( kind = core_rknd ), dimension(n_variables) :: &
      s_1j ! s_1j = sqrt(1-c_1j^2)


    !-------------------- Begin code --------------------

    ! Remove compiler warnings about unused variables.
    if ( .false. ) then
       print *, "sqrt_xp2_on_xm2 = ", sqrt_xp2_on_xm2
       print *, "f_ij_o = ", f_ij_o
    endif

    ! calculate all square roots
    do i = 1, n_variables

       s_1j(i) = sqrt(1._core_rknd-corr_matrix_approx(i,1)**2)

    end do


    ! Diagnose the missing correlations (upper triangle)
    do j = 2, (n_variables-1)
      do i = (j+1), n_variables

        ! formula (16) in the ref. paper (Larson et al. (2011))
        !f_ij = alpha_corr * sqrt_xp2_on_xm2(i) * sqrt_xp2_on_xm2(j) &
        !        * sign(1.0_core_rknd,corr_matrix_approx(1,i)*corr_matrix_approx(1,j))

        ! If the predicting c1i's are small then cij will be closer to the prescribed value. If
        ! the c1i's are bigger, then cij will be closer to formular (15) from the ref. paper. See
        ! clubb:ticket:514:comment:61 for details.
        !f_ij = (1-abs(corr_matrix_approx(1,i)*corr_matrix_approx(1,j)))*corr_matrix_prescribed(i,j) &
        !       + abs(corr_matrix_approx(1,i)*corr_matrix_approx(1,j))*f_ij_o

        f_ij = corr_matrix_prescribed(i,j)

        ! make sure -1 < f_ij < 1
        if ( f_ij < -max_mag_correlation ) then

           f_ij = -max_mag_correlation

        else if ( f_ij > max_mag_correlation ) then

           f_ij = max_mag_correlation

        end if


        ! formula (15) in the ref. paper (Larson et al. (2011))
        corr_matrix_approx(i,j) = corr_matrix_approx(i,1) * corr_matrix_approx(j,1) &
        + f_ij * s_1j(i) * s_1j(j)

      end do ! do j
    end do ! do i
    
  end subroutine diagnose_corr 

  !-----------------------------------------------------------------------
  subroutine approx_w_corr( nz, d_variables, pdf_params, & ! Intent(in)
                            rrainm, Nrm, Ncnm, &
                            stdev_w, sigma_rr_1, &
                            sigma_Nr_1, sigma_Ncn_1, &
                            corr_array) ! Intent(out)
    ! Description:
    ! Approximate the correlations of w with the hydrometeors.

    ! References:
    ! clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use pdf_parameter_module, only:  &
        pdf_parameter  ! Type

    use constants_clubb, only:  &
        one,          & ! Constant(s)
        rr_tol,       &
        Nr_tol,       &
        Ncn_tol,      &
        w_tol,        & ! [m/s]
        s_mellor_tol    ! [kg/kg]

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Number of diagnosed correlations
      nz             ! Number of model vertical grid levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >    [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >    [num/kg]
      Ncnm,            & ! Mean cloud nuclei conc., < N_cn >        [num/kg]
      stdev_w            ! Standard deviation of w                              [m/s]

    real( kind = core_rknd ), intent(in) :: &
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (2nd PDF component)   [num/kg]
      sigma_rr_1       ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables, nz), intent(out) :: &
      corr_array

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      corr_sw,       & ! Correlation between w and s (both components)       [-]
      corr_wrr,      & ! Correlation between w and rr (both components)      [-]
      corr_wNr,      & ! Correlation between w and Nr (both components)      [-]
      corr_wNcn        ! Correlation between w and Ncn (both components)     [-]

    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zt,   & ! Covariance of s and w on the zt-grid    [(m/s)(kg/kg)]
      wprrp_zt,  & ! Covariance of r_r and w on the zt-grid  [(m/s)(kg/kg)]
      wpNrp_zt,  & ! Covariance of N_r and w on the zt-grid  [(m/s)(#/kg)]
      wpNcnp_zt    ! Covariance of N_cn and w on the zt-grid  [(m/s)(#/kg)]

    real( kind = core_rknd ) :: &
      s_mellor_m,      & ! Mean of s_mellor                              [kg/kg]
      stdev_s_mellor     ! Standard deviation of s_mellor                [kg/kg]

    integer :: k ! vertical loop iterator

    ! ----- Begin Code -----

    call approx_w_covar( nz, pdf_params, rrainm, Nrm, Ncnm, & ! Intent(in)
                         wpsp_zt, wprrp_zt, wpNrp_zt, wpNcnp_zt ) ! Intent(out)

    do k = 1, nz

       s_mellor_m &
       = calc_mean( pdf_params(k)%mixt_frac, pdf_params(k)%s1, &
                    pdf_params(k)%s2 )

       stdev_s_mellor &
       = sqrt( pdf_params(k)%mixt_frac &
               * ( ( pdf_params(k)%s1 - s_mellor_m )**2 &
                     + pdf_params(k)%stdev_s1**2 ) &
             + ( one - pdf_params(k)%mixt_frac ) &
               * ( ( pdf_params(k)%s2 - s_mellor_m )**2 &
                     + pdf_params(k)%stdev_s2**2 ) &
             )

       corr_sw(k) &
       = calc_w_corr( wpsp_zt(k), stdev_w(k), stdev_s_mellor, &
                      w_tol, s_mellor_tol )

       corr_wrr(k) &
       = calc_w_corr( wprrp_zt(k), stdev_w(k), sigma_rr_1, w_tol, rr_tol )

       corr_wNr(k) &
       = calc_w_corr( wpNrp_zt(k), stdev_w(k), sigma_Nr_1, w_tol, Nr_tol )

       corr_wNcn(k) &
       = calc_w_corr( wpNcnp_zt(k), stdev_w(k), sigma_Ncn_1, w_tol, Ncn_tol )

    enddo

    call set_w_corr( nz, d_variables, & ! Intent(in)
                         corr_sw, corr_wrr, corr_wNr, corr_wNcn, &
                         corr_array ) ! Intent(inout)

  end subroutine approx_w_corr


  !-----------------------------------------------------------------------
  subroutine approx_w_covar( nz, pdf_params, rrainm, Nrm, Ncnm, & ! Intent(in)
                             wpsp_zt, wprrp_zt, wpNrp_zt, wpNcnp_zt ) ! Intent(out)
    ! Description:
    ! Approximate the covariances of w with the hydrometeors using Eddy
    ! diffusivity.

    ! References:
    ! clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: &
        gr,  & ! Variable(s)
        zm2zt,  & ! Procedure(s)
        zt2zm

    use pdf_parameter_module, only:  &
        pdf_parameter  ! Type

    use parameters_tunable, only: &
        c_Krrainm ! Variable(s)

    use constants_clubb, only: &
        one ! Constant(s)

    use advance_windm_edsclrm_module, only: &
        xpwp_fnc ! Procedure(s)

    use variables_diagnostic_module, only: &
        Kh_zm ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                         [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) ::  &
      rrainm,          & ! Mean rain water mixing ratio, < r_r >      [kg/kg]
      Nrm,             & ! Mean rain drop concentration, < N_r >      [num/kg]
      Ncnm               ! Mean cloud nuclei concentration, < N_cn >  [num/kg]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) ::  &
      wpsp_zt,   & ! Covariance of s and w on the zt-grid     [(m/s)(kg/kg)]
      wprrp_zt,  & ! Covariance of r_r and w on the zt-grid   [(m/s)(kg/kg)]
      wpNrp_zt,  & ! Covariance of N_r and w on the zt-grid   [(m/s)(#/kg)]
      wpNcnp_zt    ! Covariance of N_cn and w on the zt-grid  [(m/s)(#/kg)]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zm,   & ! Covariance of s and w on the zm-grid     [(m/s)(kg/kg)]
      wprrp_zm,  & ! Covariance of r_r and w on the zm-grid   [(m/s)(kg/kg)]
      wpNrp_zm,  & ! Covariance of N_r and w on the zm-grid   [(m/s)(#/kg)]
      wpNcnp_zm    ! Covariance of N_cn and w on the zm-grid  [(m/s)(#/kg)]

    integer :: k ! vertical loop iterator

    ! ----- Begin Code -----

    ! calculate the covariances of w with the hydrometeors
    do k = 1, nz
      wpsp_zm(k) = pdf_params(k)%mixt_frac &
                   * ( one - pdf_params(k)%mixt_frac ) &
                   * ( pdf_params(k)%s1 - pdf_params(k)%s2 ) &
                   * ( pdf_params(k)%w1 - pdf_params(k)%w2 )
    enddo

! same for wpNrp
!    wprrp_zm(1:nz-1) &
!    = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), &
!                rrainm(1:nz-1) / max( precip_frac(1:nz-1), eps ), &
!                rrainm(2:nz) / max( precip_frac(2:nz), eps ), &
!                gr%invrs_dzm(1:nz-1) )

    wprrp_zm(1:nz-1) &
    = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), &
                rrainm(1:nz-1), rrainm(2:nz), &
                gr%invrs_dzm(1:nz-1) )

    wpNrp_zm(1:nz-1) &
    = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), &
                Nrm(1:nz-1), Nrm(2:nz), &
                gr%invrs_dzm(1:nz-1) )

    wpNcnp_zm(1:nz-1) = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), Ncnm(1:nz-1), &
                                  Ncnm(2:nz), gr%invrs_dzm(1:nz-1) )

    ! Boundary conditions; We are assuming constant flux at the top.
    wprrp_zm(nz)  = wprrp_zm(nz-1)
    wpNrp_zm(nz)  = wpNrp_zm(nz-1)
    wpNcnp_zm(nz) = wpNcnp_zm(nz-1)

    ! interpolate back to zt-grid
    wpsp_zt   = zm2zt(wpsp_zm)
    wprrp_zt  = zm2zt(wprrp_zm)
    wpNrp_zt  = zm2zt(wpNrp_zm)
    wpNcnp_zt = zm2zt(wpNcnp_zm)

  end subroutine approx_w_covar

  !-----------------------------------------------------------------------
  function calc_w_corr( wpxp, stdev_w, stdev_x, w_tol, x_tol )
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
    real( kind = core_rknd ), intent(in) :: &
      stdev_w,  & ! standard deviation of w [m/s]
      stdev_x,  & ! standard deviation of x [units vary]
      wpxp,     & ! Covariances of w with the hydrometeors [units vary]
      w_tol,    & ! tolerance for w [m/s]
      x_tol       ! tolerance for x [units vary]

    real( kind = core_rknd ) :: &
      calc_w_corr

    ! --- Begin Code ---

    calc_w_corr = wpxp / ( max(stdev_x, x_tol) * max(stdev_w, w_tol) )

    ! Make sure the correlation is in [-1,1]
    if ( calc_w_corr < -max_mag_correlation ) then

      calc_w_corr = -max_mag_correlation

    else if ( calc_w_corr > max_mag_correlation ) then

      calc_w_corr = max_mag_correlation

    end if
   
  end function calc_w_corr


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
      mixt_frac, &  ! mixing ratio [-]
      x1, &         ! first component of the double gaussian [units vary]
      x2, &         ! second component of the double gaussian [units vary]
      xm, &         ! mean of x [units vary]
      x1p2, &       ! variance of the first component [units vary]
      x2p2          ! variance of the second component [units vary]

    ! Return Variable
    real( kind = core_rknd ) :: &
      calc_varnce ! variance of x (both components) [units vary]

    ! --- Begin Code ---

    calc_varnce &
    = mixt_frac * ( ( x1 - xm )**2 + x1p2 ) &
      + ( 1.0_core_rknd - mixt_frac ) * ( ( x2 - xm )**2 + x2p2 )

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
    real( kind = core_rknd ), intent(in) :: &
      mixt_frac, &  ! mixing ratio [-]
      x1, &         ! first component of the double gaussian [units vary]
      x2            ! second component of the double gaussian [units vary]

    ! Return Variable
    real( kind = core_rknd ) :: &
      calc_mean  ! mean of x (both components) [units vary]

    ! --- Begin Code ---

    calc_mean = mixt_frac * x1 + (1.0_core_rknd - mixt_frac) * x2

    return
  end function calc_mean

  !-----------------------------------------------------------------------
  subroutine set_w_corr( nz, d_variables, & ! Intent(in)
                         corr_sw, corr_wrr, corr_wNr, corr_wNcn, &
                         corr_array ) ! Intent(inout)

    ! Description:
    ! Set the first row of corr_array to the according w-correlations.

    ! References:
    ! clubb:ticket:514
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use corr_matrix_module, only: &
      iiLH_w,           & ! Variable(s)
      iiLH_s_mellor,    &
      iiLH_rrain,       &
      iiLH_Nr,          &
      iiLH_Ncn

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of Variables to be diagnosed

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      corr_sw,       & ! Correlation between s & w (both components)         [-]
      corr_wrr,      & ! Correlation between rr & w (both components)        [-]
      corr_wNr,      & ! Correlation between Nr & w (both components)        [-]
      corr_wNcn        ! Correlation between Ncn & w (both components)       [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables, d_variables, nz), &
    intent(inout) :: &
      corr_array

    ! ----- Begin Code -----

      corr_array(iiLH_w, iiLH_s_mellor, :) = corr_sw
      corr_array(iiLH_w, iiLH_rrain, :) = corr_wrr
      corr_array(iiLH_w, iiLH_Nr, :) = corr_wNr
      corr_array(iiLH_w, iiLH_Ncn, :) = corr_wNcn

  end subroutine set_w_corr

  !=============================================================================
  subroutine corr_stat_output( d_variables, nz, corr_array )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type, only: &
        stat_update_var  ! Procedure(s)

    use stats_variables, only : &
        icorr_srr,    & ! Variable(s)
        icorr_sNr,    &
        icorr_sNcn,   &
        icorr_rrNr,   &
        icorr_sw,     &
        icorr_wrr,    &
        icorr_wNr,    &
        icorr_wNcn,   &
        zt,           &
        l_stats_samp

    use corr_matrix_module, only: &
        iiLH_w,        & ! Variable(s)
        iiLH_s_mellor, &
        iiLH_t_mellor, &
        iiLH_Ncn,      &
        iiLH_rrain,    &
        iiLH_Nr

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of vertical levels
      d_variables    ! Number of variables in correlation array

    real( kind = core_rknd ), dimension( d_variables, d_variables, nz ), &
    intent(in) :: &
      corr_array    ! Correlation array


    !!! Output the correlations

    ! Statistics
    if ( l_stats_samp ) then

       ! Correlation between w and s.
       if ( iiLH_s_mellor > iiLH_w ) then
          call stat_update_var( icorr_sw, &
                                corr_array( iiLH_s_mellor, iiLH_w, : ), zt )
       else
          call stat_update_var( icorr_sw, &
                                corr_array( iiLH_w, iiLH_s_mellor, : ), zt )
       endif

       ! Correlation (in-precip) between s and r_r.
       if ( iiLH_s_mellor > iiLH_rrain ) then
          call stat_update_var( icorr_srr, &
                                corr_array( iiLH_s_mellor, iiLH_rrain, : ), zt )
       else
          call stat_update_var( icorr_srr, &
                                corr_array( iiLH_rrain, iiLH_s_mellor, : ), zt )
       endif

       ! Correlation (in-precip) between s and N_r.
       if ( iiLH_s_mellor > iiLH_Nr ) then
          call stat_update_var( icorr_sNr, &
                                corr_array( iiLH_s_mellor, iiLH_Nr, : ), zt )
       else
          call stat_update_var( icorr_sNr, &
                                corr_array( iiLH_Nr, iiLH_s_mellor, : ), zt )
       endif

       ! Correlation between s and N_cn.
       if ( iiLH_s_mellor > iiLH_Ncn ) then
          call stat_update_var( icorr_sNcn, &
                                corr_array( iiLH_s_mellor, iiLH_Ncn, : ), zt )
       else
          call stat_update_var( icorr_sNcn, &
                                corr_array( iiLH_Ncn, iiLH_s_mellor, : ), zt )
       endif

       ! Correlation (in-precip) between r_r and N_r.
       if ( iiLH_rrain > iiLH_Nr ) then
          call stat_update_var( icorr_rrNr, &
                                corr_array( iiLH_rrain, iiLH_Nr, : ), zt )
       else
          call stat_update_var( icorr_rrNr, &
                                corr_array( iiLH_Nr, iiLH_rrain, : ), zt )
       endif

       ! Correlation between w and r_r.
       if ( iiLH_w > iiLH_rrain ) then
          call stat_update_var( icorr_wrr, &
                                corr_array( iiLH_w, iiLH_rrain, : ), zt )
       else
          call stat_update_var( icorr_wrr, &
                                corr_array( iiLH_rrain, iiLH_w, : ), zt )
       endif

       ! Correlation between w and N_r.
       if ( iiLH_w > iiLH_Nr ) then
          call stat_update_var( icorr_wNr, &
                                corr_array( iiLH_w, iiLH_Nr, : ), zt )
       else
          call stat_update_var( icorr_wNr, &
                                corr_array( iiLH_Nr, iiLH_w, : ), zt )
       endif

       ! Correlation between w and N_cn.
       if ( iiLH_w > iiLH_Ncn ) then
          call stat_update_var( icorr_wNcn, &
                                corr_array( iiLH_w, iiLH_Ncn, : ), zt )
       else
          call stat_update_var( icorr_wNcn, &
                                corr_array( iiLH_Ncn, iiLH_w, : ), zt )
       endif

    endif ! l_stats_samp


    return

  end subroutine corr_stat_output

!-----------------------------------------------------------------------
  subroutine unpack_correlations( d_variables, corr_array, & ! Intent(in)
                                  corr_ws, corr_wrr, corr_wNr, corr_wNcn, &
                                  corr_st, corr_srr, corr_sNr, corr_sNcn, &
                                  corr_trr, corr_tNr, corr_tNcn, corr_rrNr )  

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_matrix_module, only: &
        iiLH_w,        & ! Variable(s)
        iiLH_s_mellor, &
        iiLH_t_mellor, &
        iiLH_rrain,    &
        iiLH_Nr,       &
        iiLH_Ncn

    implicit none

    intrinsic :: max, sqrt, transpose

    ! Input Variables
    integer, intent(in) :: &
      d_variables   ! number of diagnosed correlations

    real( kind = core_rknd ), dimension(d_variables, d_variables), intent(in) :: &
      corr_array   ! Prescribed correlations

    ! Output variables
    real( kind = core_rknd ), intent(out) :: &
      corr_ws,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_wrr,    & ! Correlation between w and rr (1st PDF component) ip [-]
      corr_wNr,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      corr_wNcn,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      corr_st,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_srr,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_sNr,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNcn,   & ! Correlation between s and Ncn (1st PDF component)   [-]
      corr_trr,    & ! Correlation between t and rr (1st PDF component) ip [-]
      corr_tNr,    & ! Correlation between t and Nr (1st PDF component) ip [-]
      corr_tNcn,   & ! Correlation between t and Ncn (1st PDF component)   [-]
      corr_rrNr      ! Correlation between rr & Nr (1st PDF component) ip  [-]

    ! ---- Begin Code ----

!    corr_ws   = corr_array(iiLH_w, iiLH_s_mellor)
!    corr_wrr  = corr_array(iiLH_w, iiLH_rrain)
!    corr_wNr  = corr_array(iiLH_w, iiLH_Nr)
!    corr_wNcn = corr_array(iiLH_w, iiLH_Ncn)
!    corr_st   = corr_array(iiLH_s_mellor, iiLH_t_mellor)
!    corr_srr  = corr_array(iiLH_s_mellor, iiLH_rrain)
!    corr_sNr  = corr_array(iiLH_s_mellor, iiLH_Nr)
!    corr_sNcn = corr_array(iiLH_s_mellor, iiLH_Ncn)
!    corr_trr  = corr_array(iiLH_t_mellor, iiLH_rrain)
!    corr_tNr  = corr_array(iiLH_t_mellor, iiLH_Nr)
!    corr_tNcn = corr_array(iiLH_t_mellor, iiLH_Ncn)
!    corr_rrNr = corr_array(iiLH_rrain, iiLH_Nr)

    corr_ws   = corr_array(iiLH_s_mellor, iiLH_w)
    corr_wrr  = corr_array(iiLH_rrain, iiLH_w)
    corr_wNr  = corr_array(iiLH_Nr, iiLH_w)
    corr_wNcn = corr_array(iiLH_Ncn, iiLH_w)
    corr_st   = corr_array(iiLH_t_mellor, iiLH_s_mellor)
    corr_srr  = corr_array(iiLH_rrain, iiLH_s_mellor)
    corr_sNr  = corr_array(iiLH_Nr, iiLH_s_mellor)
    corr_sNcn = corr_array(iiLH_Ncn, iiLH_s_mellor)
    corr_trr  = corr_array(iiLH_rrain, iiLH_s_mellor)
    corr_tNr  = corr_array(iiLH_Nr, iiLH_s_mellor)
    corr_tNcn = corr_array(iiLH_Ncn, iiLH_s_mellor)
    corr_rrNr = corr_array(iiLH_rrain, iiLH_Nr)

  end subroutine unpack_correlations

!===============================================================================

end module diagnose_correlations_module
