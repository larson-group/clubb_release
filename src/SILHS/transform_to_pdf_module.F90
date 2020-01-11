!-------------------------------------------------------------------------------
!$Id$
!===============================================================================
module transform_to_pdf_module

  implicit none

  public :: ltqnorm, multiply_Cholesky, transform_uniform_samples_to_pdf, chi_eta_2_rtthl

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine transform_uniform_samples_to_pdf &
             ( nz, num_samples, pdf_dim, d_uniform_extra, & ! In
               Sigma_Cholesky1, Sigma_Cholesky2, & ! In
               mu1, mu2, X_mixt_comp_all_levs, & ! In
               X_u_all_levs, cloud_frac, & ! In
               l_in_precip_all_levs, & ! In
               X_nl_all_levs ) ! Out
! Description:
!   This subroutine transforms a uniform sample to a sample from CLUBB's PDF.

! References:
!   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:uniform2pdf
!
!   ``Supplying Local Microphysical Parameterizations with Information about
!     Subgrid Variability: Latin Hypercube Sampling'', JAS Vol. 62,
!     p. 4010--4026, Larson, et al. 2005
!-------------------------------------------------------------------------------

    use array_index, only: &
      iiPDF_chi, & ! Variable(s)
      iiPDF_eta, &
      iiPDF_w

    use constants_clubb, only:  &
      one, &
      zero

    use clubb_precision, only: &
      core_rknd
      
    use array_index, only: &
      iiPDF_Ncn      ! Variable

    implicit none

    ! External
    intrinsic :: max

    ! Input Variables
    integer, intent(in) :: &
      nz, & ! Number of vertical grid levels
      num_samples,  & ! Number of subcolumn samples
      pdf_dim, &  ! `d' Number of variates (normally 3 + microphysics specific variables)
      d_uniform_extra ! Number of variates included in uniform sample only (often 2)
      
    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim,nz), intent(in) :: &
      Sigma_Cholesky1, & ! Correlations Cholesky matrix, 1st component [-]
      Sigma_Cholesky2    ! Correlations Cholesky matrix, 2nd component [-]

    real( kind = core_rknd ), dimension(pdf_dim,nz), intent(in) :: &
      mu1, & ! Means of the hydrometeors,(chi, eta, w, <hydrometeors>), 1st component [units vary]
      mu2    ! Means of the hydrometeors,(chi, eta, w, <hydrometeors>), 2nd component [units vary]

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,pdf_dim+d_uniform_extra) :: &
      X_u_all_levs ! Sample drawn from uniform distribution from a particular grid level

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples) :: &
      cloud_frac   ! Cloud fraction [-]
      
    logical, intent(in), dimension(nz,num_samples) :: &
      l_in_precip_all_levs ! Whether we are in precipitation (T/F)
      
    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in the first or second mixture component
      
    ! Output Variable

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    ! Local Variables

    logical, dimension(pdf_dim) :: &
      l_d_variable_lognormal ! Whether a given variable in X_nl has a lognormal dist.

    integer :: i, k, sample
    
    real( kind = core_rknd ), dimension(nz,num_samples,pdf_dim) :: &
      std_normal ! vector of d-variate standard normal distribution [-]

    ! Flag to clip sample point values of chi in extreme situations.
    logical, parameter :: &
      l_clip_extreme_chi_sample_pts = .true.

    ! ---- Begin Code ----

    ! Determine which variables are a lognormal distribution
    i = max( iiPDF_chi, iiPDF_eta, iiPDF_w )
    l_d_variable_lognormal(1:i) = .false. ! The 1st 3 variates
    l_d_variable_lognormal(i+1:pdf_dim) = .true.  ! Hydrometeors

    !---------------------------------------------------------------------------
    ! Generate a set of sample points for a microphysics/radiation scheme
    !---------------------------------------------------------------------------

    !$acc data create(std_normal) async(1)
              
    ! From Latin hypercube sample, generate standard normal sample
    call cdfnorminv( pdf_dim, nz, num_samples, X_u_all_levs, &  ! In
                     std_normal )                               ! Out   
    
    ! Computes the nonstd_normal from the Cholesky factorization of Sigma, std_normal, and mu.
    call multiply_Cholesky( nz, num_samples, pdf_dim, std_normal, & ! In
                            Sigma_Cholesky1, Sigma_Cholesky2,     & ! In
                            mu1, mu2, X_mixt_comp_all_levs,       & ! In
                            X_nl_all_levs )                         ! Out
                            
    !$acc kernels default(present) async(1)
    
    ! Convert lognormal variates (e.g. Ncn and rr) to lognormal
    do i = max( iiPDF_chi, iiPDF_eta, iiPDF_w )+1, pdf_dim
      do sample = 1, num_samples
        do k = 1, nz
            ! Convert lognormal variates (e.g. Ncn and rr) to lognormal
            X_nl_all_levs(k,sample,i) = exp( X_nl_all_levs(k,sample,i) )
        end do
      end do
    end do
    
    ! Zero precipitation hydrometeors if not in precipitation
    do sample = 1, num_samples
      do k = 1, nz 
        
        ! Zero precipitation hydrometeors if not in precipitation
        if ( .not. l_in_precip_all_levs(k,sample) ) then
          
          X_nl_all_levs(k,sample,iiPDF_Ncn+1:pdf_dim) = zero
                  
        end if            
              
      end do
    end do

    ! Clip extreme sample point values of chi, when necessary.
    ! The values of PDF component cloud fraction have been clipped within PDF
    ! closure under extreme conditions.  This code forces the sample point
    ! values of chi to be saturated or unsaturated to match the condition
    ! enforced by the clipping of PDF component cloud fraction.
    if ( l_clip_extreme_chi_sample_pts ) then
      
      do sample = 1, num_samples
        do k = 1, nz 

          if ( cloud_frac(k,sample) < epsilon( cloud_frac(k,sample) ) ) then

             ! Cloud fraction in the 1st PDF component is 0.
             ! All sample point values of chi must be <= 0.
             ! Clip the sample point value of chi back to 0.
             X_nl_all_levs(k,sample,iiPDF_chi) = min( X_nl_all_levs(k,sample,iiPDF_chi), zero )

          elseif ( cloud_frac(k,sample) > ( one - epsilon( cloud_frac(k,sample) ) ) ) then

             ! Cloud fraction in the 1st PDF component is 1.
             ! All sample point values of chi must be > 0.
             ! Clip the sample point value of chi to epsilon.
             X_nl_all_levs(k,sample,iiPDF_chi) = max( X_nl_all_levs(k,sample,iiPDF_chi), &
                                                      epsilon( zero ) )

          endif ! cloud_frac_1

       end do
     end do

   endif ! l_clip_extreme_chi_sample_pts

   !$acc end kernels
   !$acc end data

    return
  end subroutine transform_uniform_samples_to_pdf

!-----------------------------------------------------------------------
  subroutine cdfnorminv( pdf_dim, nz, num_samples, X_u_all_levs, &
                         std_normal )
! Description:
!     This function computes the inverse of the cumulative normal distribution function.
!     The return value is the lower tail quantile for the standard normal distribution. 
!     This is equivalent to SQRT(2) * ERFINV(2*P-1), but is designed for computational 
!     efficiency on GPUs, however it also has a signficant performance boost when run 
!     on CPUs compared to the previously used ltqnorm. The GPU based performance mainly
!     comes from the reduction of the chance for warp divergence. 
!   
!     THIS FUNCTION ONLY HAS SINGLE PRECISION ACCURACY, BUT ACCEPTS DOUBLE PRECISION ARGUMENTS
!
! References:
!   This algorithm was designed based on the source code provided in 
!     M.B. Giles (2010) "Approximating the erfinv function"
!     https://people.maths.ox.ac.uk/gilesm/files/gems_erfinv.pdf
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one,  & ! Constants
      two,  &
      sqrt_2   
      
    implicit none

    ! ---------------- Input Variable(s) ----------------

    integer, intent(in) :: &
      nz,           & ! Number of vertical grid levels
      num_samples,  & ! Number of subcolumn samples
      pdf_dim         ! `d' Number of variates (normally 3 + microphysics specific variables)

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,pdf_dim) :: X_u_all_levs

    ! ---------------- Return Variable ----------------

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples,pdf_dim) :: std_normal

    ! ---------------- Local Variable(s) ----------------
    
    ! Polynomial coefficients
    real( kind = core_rknd ), dimension(9), parameter :: &
      a = (/ 2.81022636e-8, 3.43273939e-7, -3.5233877e-6, -4.39150654e-6, 0.00021858087, &
             -0.00125372503, -0.00417768164, 0.246640727, 1.50140941 /)
             
    ! Polynomial coefficients
    real( kind = core_rknd ), dimension(9), parameter :: &
      b = (/ -0.000200214257,0.000100950558,0.00134934322,-0.00367342844,0.00573950773,&
             -0.0076224613,0.00943887047,1.00167406,2.83297682 /)
             
    real( kind = core_rknd ) :: w, x
    
    integer :: &
      i, sample, k  ! Loop variables
    
    ! ---------------- Begin Code ----------------
    
    !$acc parallel loop collapse(3) async(1)
    do i = 1, pdf_dim
      do sample = 1, num_samples
        do k = 1, nz
    
          x = two * X_u_all_levs(k,sample,i) - one
          
          w = -log( ( one - x ) * ( one + x ) ) 
          
          if ( w < 5.0 ) then 
            w = w - 2.5_core_rknd
            std_normal(k,sample,i) = sqrt_2 * x &
                                   * (((((((( a(1) * w + a(2) ) * w + a(3) ) * w + a(4) ) * w &
                                        + a(5) ) * w + a(6) ) * w + a(7) ) * w + a(8) ) * w + a(9) )
          else
            w = sqrt(w) - 3._core_rknd
            std_normal(k,sample,i) = sqrt_2 * x &
                                   * (((((((( b(1) * w + b(2) ) * w + b(3) ) * w + b(4) ) * w &
                                        + b(5) ) * w + b(6) ) * w + b(7) ) * w + b(8) ) * w + b(9) )
          end if                 
          
        end do
      end do
    end do
                
  end subroutine cdfnorminv

!-----------------------------------------------------------------------
  function ltqnorm( p_core_rknd )
! Description:
!   This function is ported to Fortran from the same function written in Matlab,
!    see the following description of this function.  Hongli Jiang, 2/17/2004
!   Converted to double precision by Vince Larson 2/22/2004;
!    this improves results for input values of p near 1.

! LTQNORM Lower tail quantile for standard normal distribution.
!
!   Z = LTQNORM(P) returns the lower tail quantile for the standard normal
!   distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
!   where X has a standard normal distribution.
!
!   LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
!   more accurate value when P is close to zero.

!   The algorithm uses a minimax approximation by rational functions and the
!   result has a relative error less than 1.15e-9. A last refinement by
!   Halley's rational method is applied to achieve full machine precision.

!   Author:      Peter J. Acklam
!   Time-stamp:  2003-04-23 08:26:51 +0200
!   E-mail:      pjacklam@online.no
!   URL:         http://home.online.no/~pjacklam
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd, &    ! Constant(s)
      dp ! double precision

    use constants_clubb, only: &    
      sqrt_2_dp,   &        ! Constant(s)
      sqrt_2pi_dp, &
      two_dp,      &
      one_dp,      &
      one_half_dp

#ifdef CLUBB_CAM
    ! Some compilers cannot handle 1.0/0.0, so in CAM we import their
    ! +Inf and -Inf constants. We REALLY should find a better way to
    ! do this.
    ! Eric Raut, 24 Feb 2016
    use shr_infnan_mod, only: &
      nan => shr_infnan_nan, &
      infp => shr_infnan_posinf, &
      infn => shr_infnan_neginf, &
      assignment(=)
#endif

    implicit none

    ! External

    intrinsic :: log, sqrt, exp

    ! Constant Parameter

    ! Apply Halley's method to answer to achieve more accurate result
    logical, parameter :: &
      l_apply_halley_method = .true.

    ! Input Variable(s)

    real( kind = core_rknd ), intent(in) :: p_core_rknd

    ! Return Variable

    real( kind = core_rknd ) :: ltqnorm

    ! Local Variable(s)
    real( kind = dp ) :: p

    real( kind = dp ) a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, &
                     c1, c2, c3, c4, c5, c6, d1, d2, d3, d4

    real( kind = dp ) q, r, z, z1, plow, phigh

    real( kind = dp ) ::  e, u

! Coefficients in rational approximations.
! equivalent: a(1)=a1, a(2)=a2, and etc, when a(1) is in Matlab.
! Similarly for b, c, and d's
    parameter (a1 = -3.969683028665376E+01_dp,  &
               a2 = 2.209460984245205E+02_dp, &
               a3 = -2.759285104469687E+02_dp,  &
               a4 = 1.383577518672690E+02_dp, &
               a5 = -3.066479806614716E+01_dp,  &
               a6 = 2.506628277459239E+00_dp)
    parameter (b1 = -5.447609879822406E+01_dp,  &
               b2 = 1.615858368580409E+02_dp, &
               b3 = -1.556989798598866E+02_dp,  &
               b4 = 6.680131188771972E+01_dp, &
               b5 = -1.328068155288572E+01_dp)
    parameter (c1 = -7.784894002430293E-03_dp,  &
               c2 = -3.223964580411365E-01_dp, &
               c3 = -2.400758277161838E+00_dp,  &
               c4 = -2.549732539343734E+00_dp, &
               c5 =  4.374664141464968E+00_dp,  &
               c6 =  2.938163982698783E+00_dp)
    parameter (d1 =  7.784695709041462E-03_dp,  &
               d2 =  3.224671290700398E-01_dp, &
               d3 =  2.445134137142996E+00_dp,  &
               d4 =  3.754408661907416E+00_dp)

    p = real( p_core_rknd, kind=dp )

    ! Default initialization
    z = 0.0_dp

!  Define break-points.
    plow  = 0.02425_dp
    phigh = 1._dp - plow

!  Initialize output array. Don't need this in Fortran
!   z = zeros(size(p));

!  Rational approximation for lower region:
    if (p > 0._dp .and. p < plow) then
      q = sqrt( -2._dp * log( p ) )
      z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ &
                ((((d1*q+d2)*q+d3)*q+d4)*q+1._dp)
!  Rational approximation for central region:
    else if (p >= plow .and. p <= phigh) then
      q = p - 0.5_dp
      r = q * q
      z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q &
                 /(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1._dp)
! Rational approximation for upper region:
    else if (p > phigh .and. p < 1._dp) then
      q  = sqrt( -2._dp * log(1._dp - p) )
      z  = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) &
                  /((((d1*q+d2)*q+d3)*q+d4)*q+1._dp)
    end if

    ! Eric Raut note: In CAM, we use CAM's predefined infinity and nan
    ! constants to avoid dividing by zero. We don't have similar constants
    ! in CLUBB or SILHS "cores", so we have to divide by zero. We should
    ! fix this. --24 Feb 2016
#ifdef CLUBB_CAM
! Case when P = 1:, z=+inf
    if(p == 1._dp)then
       z = infp
    end if

!  Case when P = 0: z = -inf
    if (p == 0._dp) then
       z = infn
    end if

!  Cases when output will be NaN:
!   k = p < 0 | p > 1 | isnan(p);
    if (p < 0._dp .or. p > 1._dp) then
       z = nan
    end if
#else
!  Case when P = 0: z = -inf, to create inf z =-1.0.
!     to create NaN's inf*inf.
    z1 = 0._dp
    if (p == 0._dp) then
      z = (-1._dp)/z1
    end if

! Case when P = 1:, z=inf
    if(p == 1._dp)then
      z = 1._dp/z1
    end if

!  Cases when output will be NaN:
!   k = p < 0 | p > 1 | isnan(p);
! usually inf*inf --> NaN's.
    if (p < 0._dp .or. p > 1._dp) then
      z = (1._dp/z1)**2
    end if
#endif

!  The relative error of the approximation has absolute value less
!  than 1.15e-9. One iteration of Halley's rational method (third
!  order) gives full machine precision.
! V. Larson 20Feb04: Don't use the following if-end if loop.
!   The value of e is very different than what MATLAB produces,
!   possibly because of
!   poor values of erf from Numerical Recipes.
!   The value is close to MATLAB's
!   if I omit the following if-end if loop.
! End V. Larson comment

    ! Halley's rational method is applied to achieve a more accurate result if
    ! the flag below is true. In tests, this did increase the runtime of SILHS
    ! slightly but did improve results.
    ! Eric Raut 23Aug14
    if ( l_apply_halley_method ) then
      e = one_half_dp * erfc(-z/sqrt_2_dp) - p
      u = e * sqrt_2pi_dp * exp( (z**2) / two_dp )
      z = z - u / ( one_dp + z*u/two_dp )
    end if

! return z as double precision:
    ltqnorm = real( z, kind=core_rknd )

    return
  end function ltqnorm

!-------------------------------------------------------------------------------
  subroutine multiply_Cholesky( nz, num_samples, pdf_dim, std_normal, &
                                Sigma_Cholesky1, Sigma_Cholesky2, &
                                mu1, mu2, X_mixt_comp_all_levs, &
                                X_nl_all_levs )
! Description:
!   Computes X_nl_all_levs from the Cholesky factorization of Sigma,
!   std_normal, and mu.
!   X_nl_all_levs = Sigma_Cholesky * std_normal + mu.

! References:
!   M. E. Johnson (1987), ``Multivariate Normal and Related Distributions'' p50-55
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    implicit none
    
    ! Input Variables
    integer, intent(in) :: &
      nz,           & ! Number of vertical grid levels
      num_samples,  & ! Number of samples
      pdf_dim         ! Number of variates (normally=5)

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,pdf_dim) :: &
      std_normal ! vector of d-variate standard normal distribution [-]

    real( kind = core_rknd ), intent(in), dimension(pdf_dim,nz) :: &
      mu1, & ! d-dimensional column vector of means of Gaussian, 1st component [units vary]
      mu2    ! d-dimensional column vector of means of Gaussian, 2nd component [units vary]

    real( kind = core_rknd ), intent(in), dimension(pdf_dim,pdf_dim,nz) :: &
      Sigma_Cholesky1, & ! Cholesky factorization of the Sigma matrix, 1st component [units vary]
      Sigma_Cholesky2    ! Cholesky factorization of the Sigma matrix, 2nd component [units vary]
      
    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in the first or second mixture component

    ! Output Variables

    ! nxd matrix of n samples from d-variate normal distribution
    !   with mean mu and covariance structure Sigma
    real( kind = core_rknd ), intent(out), dimension(nz,num_samples,pdf_dim) :: &
      X_nl_all_levs
      
    ! Local Variables

    ! Loop iterators
    integer :: i, j, k, sample

    ! --- Begin Code ---
    
    !$acc kernels default(present) async(1)

    do i = 1, pdf_dim
      do sample = 1, num_samples
        do k = 1, nz
          X_nl_all_levs(k,sample,i) = 0.0_core_rknd ! Copy std_normal into 'x'                      
        end do
      end do
    end do
    
    ! Compute Sigma_Cholesky * std_normal
    do sample = 1, num_samples
      do k = 1, nz
        do  i = 1, pdf_dim
          do j = 1, i
            
            if ( X_mixt_comp_all_levs(k,sample) == 1 ) then
              X_nl_all_levs(k,sample,i) = X_nl_all_levs(k,sample,i) &
                                        + Sigma_Cholesky1(i,j,k) * std_normal(k,sample,j)
            else
              X_nl_all_levs(k,sample,i) = X_nl_all_levs(k,sample,i) &
                                        + Sigma_Cholesky2(i,j,k) * std_normal(k,sample,j)
            end if
             
          end do
        end do
      end do
    end do
     
    do i = 1, pdf_dim
      do sample = 1, num_samples
        do k = 1, nz
          ! Add mu to Sigma * std_normal
          if ( X_mixt_comp_all_levs(k,sample) == 1 ) then
            X_nl_all_levs(k,sample,i) = X_nl_all_levs(k,sample,i) + mu1(i,k)
          else
            X_nl_all_levs(k,sample,i) = X_nl_all_levs(k,sample,i) + mu2(i,k)
          end if
        end do
      end do
    end do
    
    !$acc end kernels

    return
  end subroutine multiply_Cholesky
!-----------------------------------------------------------------------
  subroutine chi_eta_2_rtthl( nz, num_samples, &
                              rt_1, thl_1, &
                              rt_2, thl_2, &
                              crt_1, cthl_1, &
                              crt_2, cthl_2, &
                              mu_chi_1, mu_chi_2, &
                              chi, eta, &
                              X_mixt_comp_all_levs, &
                              lh_rt, lh_thl )
! Description:
!   Converts from chi(s), eta(t) variables to rt, thl.  Also sets a limit on the value
!   of cthl_1 and cthl_2 to prevent extreme values of temperature.
!
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! double precision

    implicit none

    ! External

    intrinsic :: max, real

    ! Constant Parameters

    real(kind = core_rknd), parameter :: &
      thl_dev_lim = 5.0_core_rknd ! Max deviation from mean thetal [K]

    ! ------------------- Input Variables -------------------
    
    integer, intent(in) :: &
      nz, &         ! Vertical grid levels
      num_samples   ! Number of subcolumn samples

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rt_1, rt_2,    & ! n dimensional column vector of rt         [kg/kg]
      thl_1, thl_2,  & ! n dimensional column vector of thetal     [K]
      crt_1, crt_2,  & ! Constants from plumes 1 & 2 of rt
      cthl_1, cthl_2   ! Constants from plumes 1 & 2 of thetal

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      mu_chi_1, mu_chi_2 ! Mean for chi_1 and chi_2         [kg/kg]

    ! n-dimensional column vector of Mellor's chi(s) and eta(t), including mean and perturbation
    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      chi, &  ! [kg/kg]
      eta     ! [-]

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in the first or second mixture component

    ! ------------------- Output variables -------------------

    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      lh_rt, lh_thl ! n-dimensional column vectors of rt and thl, including mean and perturbation

    ! ------------------- Local Variables -------------------

    real( kind= core_rknd ) :: lh_dev_thl_lim ! Limited value of the deviation on thetal [K]
    
    integer :: k, sample  ! Loop indices

    ! ---- Begin Code ----
    
    do sample = 1, num_samples
      do k = 2, nz

        if ( X_mixt_comp_all_levs(k,sample) == 1 ) then
          
          lh_rt(k,sample)  = rt_1(k) + (0.5_core_rknd/crt_1(k))*(chi(k,sample)-mu_chi_1(k)) +  &
                                 (0.5_core_rknd/crt_1(k))*eta(k,sample)

          ! Limit the quantity that temperature can vary by (in K)
          lh_dev_thl_lim = (-0.5_core_rknd/cthl_1(k))*(chi(k,sample)-mu_chi_1(k)) &
                         + (0.5_core_rknd/cthl_1(k))*eta(k,sample)

          lh_dev_thl_lim = max( min( lh_dev_thl_lim, thl_dev_lim ), -thl_dev_lim )

          lh_thl(k,sample) = thl_1(k) + lh_dev_thl_lim

        else 
          
          ! Mixture fraction 2
          lh_rt(k,sample) = rt_2(k) + (0.5_core_rknd/crt_2(k))*(chi(k,sample)-mu_chi_2(k)) +  &
                                 (0.5_core_rknd/crt_2(k))*eta(k,sample)

          ! Limit the quantity that temperature can vary by (in K)
          lh_dev_thl_lim = (-0.5_core_rknd/cthl_2(k))*(chi(k,sample)-mu_chi_2(k)) &
                         + (0.5_core_rknd/cthl_2(k))*eta(k,sample)

          lh_dev_thl_lim = max( min( lh_dev_thl_lim, thl_dev_lim ), -thl_dev_lim )

          lh_thl(k,sample) = thl_2(k) + lh_dev_thl_lim

        end if
      
      end do
    end do

    return
    
  end subroutine chi_eta_2_rtthl

end module transform_to_pdf_module
