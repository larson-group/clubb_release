!$Id$

module estimate_lh_micro_mod

  implicit none

  public :: estimate_lh_micro

  private :: autoconv_driver, rc_estimate, lh_microphys_driver

  private ! Default Scope

  contains

!------------------------------------------------------------------------

  subroutine estimate_lh_micro &
             ( dt, nnzp, n_micro_calls, d_variables, &
               X_u_all_levs, X_nl_all_levs, & 
               rt, thl, l_sample_flag, pdf_params, & 
               thlm, p_in_Pa, exner, rho, &
               wm, w_std_dev, dzq, rcm, rvm, &        
               cloud_frac, hydromet, &
               lh_hydromet_mc, lh_hydromet_vel, &
               lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, &
               lh_AKm, AKm, AKstd, AKstd_cld, & 
               AKm_rcm, AKm_rcc, lh_rcm_avg, &
               lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &
               lh_wm, lh_Ncp2_zt, lh_Nrp2_zt, lh_rrainp2_zt, lh_rcp2_zt, &
               lh_wp2_zt, lh_rtp2_zt, lh_thlp2_zt, &
               lh_cloud_frac, &
               microphys_sub )
! Description:
!   This subroutine computes microphysical grid box averages,
!   given a Latin Hypercube sample.
! References:
!   None
!------------------------------------------------------------------------

    use constants, only:  &
      pi,  & ! Variables(s)
      sstol, &
      zero_threshold

    use anl_erf, only:  &
      erf ! Procedure(s)

    use variables_prognostic_module, only:  &
      pdf_parameter  ! type

    use parameters_model, only: &
      hydromet_dim

    implicit none

    ! External
#include "../microphys_interface.inc"

    ! Input Variables

    real, intent(in) :: dt ! Model timestep     [s]

    integer, intent(in) :: &
      nnzp, &          ! Number of vertical levels
      n_micro_calls, & ! Number of calls to the microphysics
      d_variables      ! Number of variates

    ! Sample drawn from uniform distribution
    double precision, dimension(nnzp,n_micro_calls,d_variables+1), intent(in) :: &
      X_u_all_levs ! N x D+1 Latin hypercube sample from uniform dist

    double precision, dimension(nnzp,n_micro_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    double precision, dimension(nnzp,n_micro_calls), intent(in) :: &
      rt, thl ! Total water / Liquid potential temperature      [kg/kg] / [K]

    ! Flag that determines whether we have a special case (false)
    logical, dimension(nnzp), intent(in) :: l_sample_flag

    real, dimension(nnzp), intent(in) :: &
      cloud_frac, & ! Cloud fraction           [-]
      thlm,       & ! Liquid pot. temperature  [K]
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho           ! Density on thermo. grid  [kg/m^3]

    real, dimension(nnzp), intent(in) :: &
      wm, &        ! Mean w                             [m/s]
      w_std_dev, & ! Standard deviation of w            [m/s]
      dzq          ! Difference in height per gridbox   [m]

    real, dimension(nnzp), intent(in) :: &
      rcm, & ! Liquid water mixing ratio        [kg/kg]
      rvm    ! Vapor water mixing ratio         [kg/kg]

    type(pdf_parameter), intent(in) :: pdf_params

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Output Variables
    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      lh_rcm_mc, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc   ! LH estimate of time tendency of liquid potential temperature [K/s]

    real, dimension(nnzp), intent(out) :: &
      lh_AKm,   & ! Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
      AKm,       & ! Exact Kessler autoconversion, AKm,             [kg/kg/s]
      AKstd,     & ! Exact standard deviation of gba Kessler        [kg/kg/s]
      AKstd_cld, & ! Exact w/in cloud std of gba Kessler            [kg/kg/s]
      AKm_rcm,   & ! Exact local gba Kessler auto based on rcm      [kg/kg/s]
      AKm_rcc      ! Exact local gba Kessler based on w/in cloud rc [kg/kg/s]

    ! For comparison, estimate kth liquid water using Monte Carlo
    real, dimension(nnzp), intent(out) :: &
      lh_rcm_avg ! LH estimate of grid box avg liquid water [kg/kg]

    real, dimension(nnzp,hydromet_dim), intent(out) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp), intent(out) :: &
      lh_thlm,      & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,       & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,       & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,        & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt,    & ! Average value of the variance of the LH est. of vertical velocity [m^2/s^2]
      lh_rcp2_zt,   & ! Average value of the variance of the LH est. of rc             [(kg/kg)^2]
      lh_rtp2_zt,   & ! Average value of the variance of the LH est. of rt             [kg^2/kg^2]
      lh_thlp2_zt,  & ! Average value of the variance of the LH est. of thetal         [K^2]
      lh_rrainp2_zt,& ! Average value of the variance of the LH est. of rrain          [(kg/kg)^2]
      lh_Nrp2_zt,   & ! Average value of the variance of the LH est. of Nr             [#^2/kg^2]
      lh_Ncp2_zt,   & ! Average value of the variance of the LH est. of Nc             [#^2/kg^2]
      lh_cloud_frac   ! Average value of the latin hypercube est. of cloud fraction    [-]

    ! Local Variables

    ! Level on which calculations are occuring
    integer :: level

    ! PDF parameters
    real :: a
!   real :: w1, w2
!   real :: sw1, sw2
!   real :: thl1, thl2, sthl1, sthl2
!   real :: rt1,rt2
!   real :: srt1, srt2
    real :: ss1, ss2, s1, s2
    real :: cloud_frac1, cloud_frac2
!   real :: rc1, rc2

    ! Cloud fraction 0<cloud_frac<1, mean liquid water mix ratio [kg/kg]
!   real :: cloud_frac, rcm

    ! Double precision version of Monte Carlo Kessler ac
    double precision :: lh_AKm_dp
    ! Double precision version of Monte Carlo avg liquid water
    double precision :: lh_rcm_avg_dp

    double precision, dimension(n_micro_calls) :: &
      rcm_sample ! Sample points of rcm         [kg/kg]

    ! Variables needed for exact Kessler autoconversion, AKm
    real :: r_crit, K_one
    real :: sn1_crit, cloud_frac1_crit, sn2_crit, cloud_frac2_crit
    real :: AK1, AK2

    ! Variables needed for exact std of Kessler autoconversion, AKstd
    !      and within cloud standard deviation, AKstd_cld
    real :: AK1var, AK2var

    ! For comparison, compute within-cloud vertical velocity analytically.
    !real C_w_cld1, C_w_cld2, w_cld_avg

    ! ---- Begin Code ----

    ! Boundary condition
    lh_AKm(1)   = 0.0
    AKm(1)       = 0.0
    AKm_rcm(1)   = 0.0
    AKm_rcc(1)   = 0.0
    lh_rcm_avg(1)   = 0.0
    AKstd(1)     = 0.0
    AKstd_cld(1) = 0.0

    do level = 2, nnzp, 1
      ! Extract PDF parameters

      !w1         = pdf_params%w1(level)
      !w2         = pdf_params%w2(level)
      !sw1        = pdf_params%sw1(level)
      !sw2        = pdf_params%sw2(level)
      !rt1        = pdf_params%rt1(level)
      !rt2        = pdf_params%rt2(level)
      !srt1       = pdf_params%srt1(level)
      !srt2       = pdf_params%srt2(level)
      !thl1       = pdf_params%thl1(level)
      !thl2       = pdf_params%thl2(level)
      !sthl1      = pdf_params%sthl1(level)
      !sthl2      = pdf_params%sthl2(level)
      a           = pdf_params%a(level)
!     rc1         = pdf_params%rc1(level)
!     rc2         = pdf_params%rc2(level)
!     cloud_frac1 = pdf_params%cloud_frac1(level)
!     cloud_frac2 = pdf_params%cloud_frac2(level)
      cloud_frac1 = 1.0 ! For in and out of cloud sampling -dschanen 30 Jul 09
      cloud_frac2 = 1.0 !     "    "
      s1          = pdf_params%s1(level)
      s2          = pdf_params%s2(level)
      ss1         = pdf_params%ss1(level)
      ss2         = pdf_params%ss2(level)

      ! Compute mean cloud fraction and cloud water

!     cloud_frac = a * cloud_frac1 + (1-a) * cloud_frac2
!     rcm        = a * rc1 + (1-a) * rc2

      !------------------------------------------------------------------------
      ! Call Kessler autoconversion microphysics using Latin hypercube sample
      ! This acts as an interface between the boundary layer scheme
      !    and the microphysics.  To add a call to a microphysics scheme,
      !   alter two lines in autoconversion_driver.f.
      ! Then we compute Kessler ac analytically.
      !------------------------------------------------------------------------

      ! We prognose rt-thl-w,
      !    but we set means, covariance of hydromet mixing ratio's and number
      !    concentrations to constants.

      if ( .not. l_sample_flag(level) ) then
        rcm_sample(1:n_micro_calls) = dble( rcm(level) )
      else
        rcm_sample(1:n_micro_calls) = max( X_nl_all_levs(level,1:n_micro_calls,1), 0.d0)
      end if

      ! Call microphysics, i.e. Kessler autoconversion.
      ! A_K = (1e-3/s)*(rc-0.5kg/kg)*H(rc-0.5kg/kg)
      call autoconv_driver &
           ( n_micro_calls, d_variables, dble( a ), dble( cloud_frac1 ), dble( cloud_frac2 ), &
             rcm_sample, & 
               !X_nl(1:n,3), X_nl(1:n,4), X_nl(1:n,5),
             X_u_all_levs(level,:,:), lh_AKm_dp )

      ! Convert to real number
      lh_AKm(level) = real( lh_AKm_dp )

      ! Compute Monte Carlo estimate of liquid for test purposes.
      call rc_estimate &
           ( n_micro_calls, d_variables, dble( a ), dble( cloud_frac1 ), &
             dble( cloud_frac2 ), rcm_sample, & 
             ! X_nl(1:n,3), X_nl(1:n,4), X_nl(1:n,5),
             X_u_all_levs(level,:,:), lh_rcm_avg_dp )

      ! Convert to real number
      lh_rcm_avg(level) = real( lh_rcm_avg_dp )

      ! A test of the scheme:
      ! Compare exact (rcm) and Monte Carlo estimates (lh_rcm_avg) of
      !     specific liq water content, rcm.
      !      print*, 'rcm=', rcm
      !       print*, 'lh_rcm_avg=', lh_rcm_avg

      ! Exact Kessler autoconversion in units of (kg/kg)/s
      !        r_crit = 0.3e-3
      !        r_crit = 0.7e-3
      r_crit            = 0.2e-3
      K_one             = 1.e-3
      sn1_crit          = (s1-r_crit)/max( ss1, sstol )
      cloud_frac1_crit  = 0.5*(1+erf(sn1_crit/sqrt(2.0)))
      AK1               = K_one * ( (s1-r_crit)*cloud_frac1_crit  & 
                         + ss1*exp(-0.5*sn1_crit**2)/(sqrt(2*pi)) )
      sn2_crit          = (s2-r_crit)/max( ss2, sstol )
      cloud_frac2_crit  = 0.5*(1+erf(sn2_crit/sqrt(2.0)))
      AK2               = K_one * ( (s2-r_crit)*cloud_frac2_crit  & 
                         + ss2*exp(-0.5*sn2_crit**2)/(sqrt(2*pi)) )
      AKm(level)        = a * AK1 + (1-a) * AK2

      ! Exact Kessler standard deviation in units of (kg/kg)/s
      ! For some reason, sometimes AK1var, AK2var are negative
      AK1var   = max( zero_threshold, K_one * (s1-r_crit) * AK1  & 
               + K_one * K_one * (ss1**2) * cloud_frac1_crit  & 
               - AK1**2  )
      AK2var   = max( zero_threshold, K_one * (s2-r_crit) * AK2  & 
               + K_one * K_one * (ss2**2) * cloud_frac2_crit  & 
               - AK2**2  )
      ! This formula is for a grid box average:
      AKstd(level)  = sqrt(    a  * ( (AK1-AKm(level))**2 + AK1var ) & 
                  + (1-a) * ( (AK2-AKm(level))**2 + AK2var ) & 
                  )
      ! This formula is for a within-cloud average:
      if ( cloud_frac(level) > 0 ) then
        AKstd_cld(level) = sqrt( max( zero_threshold,   & 
                  (1./cloud_frac(level)) * ( a  * ( AK1**2 + AK1var ) & 
                            + (1-a) * ( AK2**2 + AK2var )  & 
                            ) & 
                 - (AKm(level)/cloud_frac(level))**2  ) & 
                        )
      else
        AKstd_cld(level) = zero_threshold
      end if

      ! Kessler autoconversion, using grid box avg liquid, rcm, as input
      AKm_rcm(level) = K_one * max( zero_threshold, rcm(level)-r_crit )

      ! Kessler ac, using within cloud liquid, rcm/cloud_frac, as input
      ! We found that for small values of cloud_frac this formula
      ! can still produce NaN values and therefore added this
      ! threshold of 0.001 here. -dschanen 3 June 2009
      if ( cloud_frac(level) > 0.001 ) then
        AKm_rcc(level) = cloud_frac(level) * K_one * &
                         max( zero_threshold, rcm(level)/cloud_frac(level)-r_crit )
      else
        AKm_rcc(level) = zero_threshold
      end if

!       print*, 'a=', a
!       print*, 's1=', s1
!       print*, 's2=', s2
!       print*, 'ss1=', ss1
!       print*, 'ss2=', ss2
!       print*, 'AKm =', AKm(level)
!       print*, 'lh_AKm (estimate) =', lh_AKm(level)
!       print*, 'AK1=', AK1
!       print*, 'AK2=', AK2
!       print*, 'AK1var=', AK1var
!       print*, 'AK2var=', AK2var
!       print*, 'AKstd =', AKstd
!       print*, 'AKstd_cld =', AKstd_cld
!       print*, 'AKm_rcc =', AKm_rcc
!       print*, 'AKm_rcm =', AKm_rcm

      !------------------------------------------------------------------------
      ! Compute within-cloud vertical velocity, avg over full layer
      ! This call is a kludge: I feed w values into rc variable
      ! in autoconv_driver.
      ! Only works if coeff=expn=1 in autoconversion_driver.
      !------------------------------------------------------------------------
      !     call autoconv_driver( n, d, a, cloud_frac1, cloud_frac2, X_nl( 1:n, 3 ), &
      !                           X_nl( 1:n, 3 ), X_nl( 1:n, 4 ), &
      !                           X_nl(1:n, 5), X_u, AKm2 )

      ! Another test:
      ! Compute within-cloud vertical velocity, avgd over full domain.
      !        C_w_cld1 =  cloud_frac1*w1
      !        C_w_cld2 =  cloud_frac2*w2
      !        w_cld_avg = a * C_w_cld1 + (1-a) * C_w_cld2

      ! The following two values should match
      !       print*, 'w_cld_avg=', w_cld_avg
      !       print*, 'ac_m2=', ac_m2

    end do ! level = 2, nnzp

    call lh_microphys_driver( dt, nnzp, n_micro_calls, d_variables, &
                              l_sample_flag, rt, thl, &
                              X_nl_all_levs, X_u_all_levs, &
                              thlm, p_in_Pa, exner, rho, wm, w_std_dev, &
                              dzq, rcm, rvm, pdf_params, hydromet, &
                              lh_rvm_mc, lh_rcm_mc, lh_hydromet_mc, &
                              lh_hydromet_vel, lh_thlm_mc, &
                              lh_hydromet, lh_thlm, lh_rcm, lh_rvm, lh_wm, &
                              lh_Ncp2_zt, lh_Nrp2_zt, lh_rrainp2_zt, lh_rcp2_zt, &
                              lh_wp2_zt, lh_rtp2_zt, lh_thlp2_zt, &
                              lh_cloud_frac, &
                              microphys_sub )

    return
  end subroutine estimate_lh_micro
!-----------------------------------------------------------------------
  subroutine autoconv_driver( n_micro_calls, d_variables, a, cloud_frac1, cloud_frac2, rc, &
                             !w, Nc, rr, &
                              X_u_one_lev, ac_m )
! Description:
!   Compute Kessler grid box avg autoconversion (kg/kg)/s.
! References:
!   None
!-----------------------------------------------------------------------

    use constants, only:  &
      fstderr, & ! Constant(s)
      g_per_kg

!   use error_code, only:  &
!     clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! External
    intrinsic :: epsilon

    ! Constant parameters
    logical, parameter :: &
      l_cloud_weighted_averaging = .false.

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    double precision, intent(in) :: &
      a,                       & ! Mixture fraction of Gaussians
      cloud_frac1, cloud_frac2   ! Cloud fraction associated w/ 1st, 2nd mixture component

    double precision, dimension(n_micro_calls), intent(in) :: &
      rc !, & ! n in-cloud values of spec liq water content (when positive) [kg/kg].
!     w,  & ! n in-cloud values of vertical velocity (m/s)
!     Nc, & ! n in-cloud values of droplet number (#/mg air)
!     rr    ! n in-cloud values of specific rain content (g/kg)

    double precision, dimension(n_micro_calls,d_variables+1), intent(in) :: &
      X_u_one_lev ! N x D+1 Latin hypercube sample from uniform dist

    ! Output Variables

    ! a scalar representing grid box average autoconversion;
    ! has same units as rc; divide by total cloud fraction to obtain
    ! within-cloud autoconversion
    double precision, intent(out) :: &
      ac_m

    ! Local Variables

    integer :: sample
    integer :: n1, n2
    double precision :: ac_m1, ac_m2
    double precision :: coeff, r_crit
    ! double precision expn
    double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of a, cloud_frac1, cloud_frac2.
    if ( a > 1.0d0 .or. a < 0.0d0 ) then
      write(fstderr,*) 'Error in autoconv_driver:  ',  &
                       'mixture fraction, a, does not lie in [0,1].'
      write(fstderr,*) 'a = ', a
      stop
    end if
    if ( cloud_frac1 > 1.0d0 .or. cloud_frac1 < 0.0d0 ) then
      write(fstderr,*) 'Error in autoconv_driver:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      write(fstderr,*) 'cloud_frac1 = ', cloud_frac1
      stop
    end if
    if ( cloud_frac2 > 1.0d0 .or. cloud_frac2 < 0.0d0 ) then
      write(fstderr,*) 'Error in autoconv_driver:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      write(fstderr,*) 'cloud_frac2 = ', cloud_frac2
      stop
    end if

    ! Make sure there is some cloud.
    ! Disable this for now, so we can loop over the whole domain.
    ! -dschanen 3 June 2009
!   if ( a*cloud_frac1 < 0.001d0 .and. (1-a)*cloud_frac2 < 0.001d0 ) then
!     if ( clubb_at_least_debug_level( 1 ) ) then
!       write(fstderr,*) 'Error in autoconv_driver:  ',  &
!                        'there is no cloud or almost no cloud!'
!     end if
!   end if

    ! Autoconversion formula prefactor and exponent.
    ! These are for Kessler autoconversion in (kg/kg)/s.
    coeff  = 1.d-3
    ! expn   = 1.d0
    ! r_crit = 0.3d0
    ! r_crit = 0.7d0
    r_crit = 0.2d0 / g_per_kg

    ! Initialize autoconversion in each mixture component
    ac_m1 = 0.d0
    ac_m2 = 0.d0

    ! Initialize numbers of sample points corresponding
    !    to each mixture component
    n1 = 0
    n2 = 0

    do sample = 1, n_micro_calls

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
      fraction_1 = a*cloud_frac1/max( a*cloud_frac1+(1-a)*cloud_frac2, epsilon( a ) )
!          print*, 'fraction_1= ', fraction_1

! V. Larson change to try to fix sampling
!          if ( X_u_one_lev(sample,d_variables+1) .lt. fraction_1 ) then
!          print*, '-1+2*int((sample+1)/2)= ', -1+2*int((sample+1)/2)
!          print*, '-1+2*int((sample+1)/2)= ', int(sample)
      if ( X_u_one_lev(sample,d_variables+1) < fraction_1 ) then
! End of V. Larson fix

! Use an idealized formula to compute autoconversion
!      in mixture comp. 1
! A_K = (1e-3/s)*(rc-0.5g/kg)*H(rc-0.5g/kg)
! This is the first of two lines where
!      a user must add a new microphysics scheme.
        ac_m1 = ac_m1 + coeff*max(0.d0,rc(sample)-r_crit)
        n1 = n1 + 1
      else
! Use an idealized formula to compute autoconversion
!      in mixture comp. 2
! A_K = (1e-3/s)*(rc-0.5g/kg)*H(rc-0.5g/kg)
! This is the second and last line where
!      a user must add a new microphysics scheme.
        ac_m2 = ac_m2 + coeff*max(0.d0,rc(sample)-r_crit)
        n2 = n2 + 1
      end if

      ! Loop to get new sample
    end do ! sample = 1, n_micro_calls

!! Convert sums to averages.
!! Old code that underestimates if a plume has no sample points
!       if (n1 .eq. 0) then
!          ac_m1 = 0.d0
!       else
!          ac_m1 = ac_m1/n1
!       end if

!       if (n2 .eq. 0) then
!         ac_m2 = 0.d0
!       else
!          ac_m2 = ac_m2/n2
!       end if

    if ( n1 == 0 .and. n2 == 0 ) then
      stop 'Error:  no sample points in autoconv_driver'
    end if

    if ( l_cloud_weighted_averaging ) then
       ! Convert sums to averages.
       ! If we have no sample points for a certain plume,
       !    then we estimate the plume liquid water by the
       !    other plume's value.
      if ( .not. (n1 == 0) ) then
        ac_m1 = ac_m1/n1
      end if

      if ( .not. (n2 == 0) ) then
        ac_m2 = ac_m2/n2
      end if

      if ( n1 == 0 ) then
        ac_m1 = ac_m2
      end if

      if ( n2 == 0 ) then
        ac_m2 = ac_m1
      end if

      ! Grid box average.
      ac_m = a*cloud_frac1*ac_m1 + (1-a)*cloud_frac2*ac_m2

    else
      ac_m = ( ac_m1 + ac_m2 ) / dble( n_micro_calls )

    end if

!   print*, 'autoconv_driver: acm=', ac_m

    return
  end subroutine autoconv_driver

!-------------------------------------------------------------------------------
  subroutine lh_microphys_driver &
             ( dt, nnzp, n_micro_calls, d_variables, &
               l_sample_flag, rt, thl, &
               X_nl_all_levs, X_u_all_levs, &
               thlm, p_in_Pa, exner, rho, wm, w_std_dev, &
               dzq, rcm, rvm, pdf_params, hydromet,  &
               lh_rvm_mc, lh_rcm_mc, lh_hydromet_mc, &
               lh_hydromet_vel, lh_thlm_mc, &
               lh_hydromet, lh_thlm, lh_rcm, lh_rvm, lh_wm, &
               lh_Ncp2_zt, lh_Nrp2_zt, lh_rrainp2_zt, lh_rcp2_zt, &
               lh_wp2_zt, lh_rtp2_zt, lh_thlp2_zt, &
               lh_cloud_frac, &
               microphys_sub )
! Description:
!   Estimate the tendency of a microphysics scheme via latin hypercube sampling
! References:
!   None
!-------------------------------------------------------------------------------

    use constants, only:  &
      fstderr, &  ! Constant(s)
      rc_tol, &
      cm3_per_m3

    use parameters_model, only: &
      hydromet_dim

!   use parameters_microphys, only: &
!     Ncm_initial

    use array_index, only: &
      iirrainm, &
      iiNrm, &
      iiNcm, &
      iiLH_rrain, &
      iiLH_Nr, &
      iiLH_Nc, &
      iiLH_rt, &
      iiLH_w

    use math_utilities, only: &
      compute_variance

    use variables_prognostic_module, only: &
      pdf_parameter

    implicit none

    ! External
#include "../microphys_interface.inc"

    intrinsic :: epsilon

    ! Constant parameters
    logical, parameter :: &
      l_compute_diagnostic_average = .true., &
      l_cloud_weighted_averaging   = .false.

    ! Input Variables
    real, intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    logical, dimension(nnzp), intent(in) :: &
      l_sample_flag  ! Whether we are sampling at this level

    double precision, dimension(nnzp,n_micro_calls), intent(in) :: &
      rt, & ! n_micro_calls values of total water mixing ratio     [kg/kg]
      thl   ! n_micro_calls values of liquid potential temperature [K]

    double precision, dimension(nnzp,n_micro_calls,d_variables+1), intent(in) :: &
      X_u_all_levs ! N x D+1 Latin hypercube sample from uniform dist

    double precision, target, dimension(nnzp,n_micro_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real, dimension(nnzp), intent(in) :: &
      thlm,       & ! Liquid pot. temperature  [K]
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho           ! Density on thermo. grid  [kg/m^3]

    real, dimension(nnzp), intent(in) :: &
      wm,        & ! Mean w                     [m/s]
      w_std_dev, & ! Standard deviation of w    [m/s]
      dzq          ! Difference in height per gridbox   [m]

    real, dimension(nnzp), intent(in) :: &
      rcm, & ! Liquid water mixing ratio        [kg/kg]
      rvm    ! Vapor water mixing ratio         [kg/kg]

    type(pdf_parameter), intent(in) :: pdf_params

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Output Variables

    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      lh_rcm_mc, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc   ! LH estimate of time tendency of liquid potential temperature [K/s]

    real, dimension(nnzp,hydromet_dim), intent(out) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp), intent(out) :: &
      lh_thlm,       & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,        & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,        & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,         & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt,     & ! Average value of the variance of the LH est. of vertical vel.  [m^2/s^2]
      lh_rcp2_zt,    & ! Average value of the variance of the LH est. of rc             [kg^2/kg^2]
      lh_rtp2_zt,    & ! Average value of the variance of the LH est. of rt             [kg^2/kg^2]
      lh_thlp2_zt,   & ! Average value of the variance of the LH est. of thetal         [K^2]
      lh_rrainp2_zt, & ! Average value of the variance of the LH est. of rrain          [kg^2/kg^2]
      lh_Nrp2_zt,    & ! Average value of the variance of the LH est. of Nr             [#^2/kg^2]
      lh_Ncp2_zt,    & ! Average value of the variance of the LH est. of Nc             [#^2/kg^2]
      lh_cloud_frac    ! Average value of the latin hypercube est. of cloud fraction    [-]

    ! Local Variables
    double precision, dimension(nnzp,hydromet_dim) :: &
      lh_hydromet_mc_m1,  & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_mc_m2,  & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_m1, & ! LH est of hydrometeor sedimentation velocity [m/s]
      lh_hydromet_vel_m2    ! LH est of hydrometeor sedimentation velocity [m/s]

    double precision, dimension(nnzp) :: &
      lh_rcm_mc_m1,  & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rcm_mc_m2,  & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_m1,  & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_rvm_mc_m2,  & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_m1, & ! LH est of time tendency of liquid potential temperature [K/s]
      lh_thlm_mc_m2    ! LH est of time tendency of liquid potential temperature [K/s]

    real, dimension(nnzp,hydromet_dim,n_micro_calls) :: &
      hydromet_tmp ! Hydrometeor species    [units vary]

    real, dimension(nnzp) :: &
      s_mellor_tmp    ! 's' (Mellor 1977)            [kg/kg]

    real, dimension(nnzp,n_micro_calls) :: &
      rv_tmp,  & ! Vapor water                  [kg/kg]
      thl_tmp, & ! Liquid potential temperature [K]
      rc_tmp,  & ! Liquid water                [kg/kg]
      w_tmp      ! Vertical velocity           [m/s]

    integer, dimension(nnzp) :: n1, n2, zero

    double precision, pointer, dimension(:,:) :: &
      s_mellor,  & ! n_micro_calls values of 's' (Mellor 1977)      [kg/kg]
      w            ! n_micro_calls values of vertical velocity      [m/s]

    double precision, dimension(nnzp) :: &
      cloud_frac1, cloud_frac2, & ! Cloud fraction for plume 1,2        [-]
      a, &  ! PDF parameter a
      fraction_1

    integer :: i, k, sample

    logical :: l_error

    ! ---- Begin Code ----

    a(:)  = dble( pdf_params%a(:) )

    if ( l_cloud_weighted_averaging ) then
      cloud_frac1(:) = dble( pdf_params%cloud_frac1(:) )
      cloud_frac2(:) = dble( pdf_params%cloud_frac2(:) )
      zero(:) = 0
    else
      cloud_frac1(:) = dble( 1.0 )
      cloud_frac2(:) = dble( 1.0 )
    end if

    ! Mellor's 's' is at the same index at iiLH_rt (total water mixing ratio)
    s_mellor => X_nl_all_levs(:,:,iiLH_rt)
    w        => X_nl_all_levs(:,:,iiLH_w)


    lh_hydromet_vel(:,:) = 0.

    ! Initialize microphysical tendencies for each mixture component
    lh_hydromet_mc_m1(:,:) = 0.d0
    lh_hydromet_mc_m2(:,:) = 0.d0

    lh_hydromet_vel_m1(:,:) = 0.d0
    lh_hydromet_vel_m2(:,:) = 0.d0

    lh_rcm_mc_m1(:) = 0.d0
    lh_rcm_mc_m2(:) = 0.d0

    lh_rvm_mc_m1(:) = 0.d0
    lh_rvm_mc_m2(:) = 0.d0

    lh_thlm_mc_m1(:) = 0.d0
    lh_thlm_mc_m2(:) = 0.d0

    ! Initialize numbers of sample points corresponding
    !    to each mixture component
    n1 = 0
    n2 = 0

    ! Choose which mixture fraction we are in.
    ! Account for cloud fraction.
    ! Follow M. E. Johnson (1987), p. 56.
    fraction_1(:) = a(:)*cloud_frac1(:)/( a(:)*cloud_frac1(:)+(1.-a(:))*cloud_frac2(:) )
!   print*, 'fraction_1= ', fraction_1

    do sample = 1, n_micro_calls

      where ( l_sample_flag )
        s_mellor_tmp(:)   = real( s_mellor(:,sample) ) 

        where( s_mellor(:,sample) > 0.0 )
          rc_tmp(:,sample)  = real( s_mellor(:,sample) ) 
        else where
          rc_tmp(:,sample)= 0.0
        end where

      else where ! Use the grid box mean
        s_mellor_tmp(:) = pdf_params%a(:) * pdf_params%s1(:) + 1.0-pdf_params%a(:)*pdf_params%s2(:)
        rc_tmp(:,sample) = rcm

      end where

      where ( l_sample_flag ) 
        w_tmp(:,sample)   = real( w(:,sample) )
        rv_tmp(:,sample)  = real( rt(:,sample) ) - rc_tmp(:,sample)
        thl_tmp(:,sample) = real( thl(:,sample) )
      else where ! Use the grid box mean
        w_tmp(:,sample)   = wm
        rv_tmp(:,sample)  = rvm
        thl_tmp(:,sample) = thlm
      end where

      ! Copy the sample points into the temporary arrays
      do i = 1, hydromet_dim, 1
        if ( i == iirrainm .and. iiLH_rrain > 0 ) then
          where ( l_sample_flag )
            hydromet_tmp(:,i,sample) = real( X_nl_all_levs(:,sample,iiLH_rrain) )
          else where
            hydromet_tmp(:,i,sample) = hydromet(:,i)
          end where
        else if ( i == iiNcm .and. iiLH_Nc > 0 ) then
          ! Kluge for when we don't have correlations between Nc, other variables
!         hydromet_tmp(:,iiNcm) = Ncm_initial * cm3_per_m3 / rho
          where ( l_sample_flag )
            hydromet_tmp(:,i,sample) = real( X_nl_all_levs(:,sample,iiLH_Nc) )
          else where
            hydromet_tmp(:,i,sample) = hydromet(:,i)
          end where
        else if ( i == iiNrm .and. iiLH_Nr > 0 ) then
          where ( l_sample_flag )
            hydromet_tmp(:,i,sample) = real( X_nl_all_levs(:,sample,iiLH_Nr) )
          else where
            hydromet_tmp(:,i,sample) = hydromet(:,i)
          end where
        else ! Use the mean field, rather than a sample point
          hydromet_tmp(:,i,sample) = hydromet(:,i)
        end if
      end do

      if ( l_compute_diagnostic_average ) then
        if ( sample == 1 ) then
          lh_hydromet(:,1:hydromet_dim) = hydromet_tmp(:,1:hydromet_dim,sample)
          lh_thlm = thl_tmp(:,sample)
          lh_rcm  = rc_tmp(:,sample)
          lh_rvm  = rv_tmp(:,sample)
          lh_wm   = w_tmp(:,sample)
          where ( l_sample_flag .and. s_mellor(:,sample) > 0. ) 
            lh_cloud_frac = 1.0
          else where
            lh_cloud_frac = 0.0
          end where
        else
          lh_hydromet(:,1:hydromet_dim) = lh_hydromet(:,1:hydromet_dim) &
            + hydromet_tmp(:,1:hydromet_dim,sample)
          lh_thlm = lh_thlm + thl_tmp(:,sample)
          lh_rcm  = lh_rcm  + rc_tmp(:,sample)
          lh_rvm  = lh_rvm  + rv_tmp(:,sample)
          lh_wm   = lh_wm   + w_tmp(:,sample)
          where ( l_sample_flag .and. s_mellor(:,sample) > 0. ) lh_cloud_frac = lh_cloud_frac + 1.0
        end if
      end if

      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub &
           ( dt, nnzp, .false., .true., thl_tmp(:,sample), p_in_Pa, exner, rho, pdf_params, &
             w_tmp(:,sample), w_std_dev, dzq, rc_tmp(:,sample), s_mellor_tmp, &
             rv_tmp(:,sample), hydromet_tmp(:,:,sample), lh_hydromet_mc, &
             lh_hydromet_vel, lh_rcm_mc, lh_rvm_mc, lh_thlm_mc )

      do i = 1, hydromet_dim
        where ( X_u_all_levs(1:nnzp,sample,d_variables+1) < fraction_1 )
          lh_hydromet_vel_m1(:,i) = lh_hydromet_vel_m1(:,i) + lh_hydromet_vel(:,i)
          lh_hydromet_mc_m1(:,i) = lh_hydromet_mc_m1(:,i) + lh_hydromet_mc(:,i)
        else where
          lh_hydromet_vel_m2(:,i) = lh_hydromet_vel_m2(:,i) + lh_hydromet_vel(:,i)
          lh_hydromet_mc_m2(:,i) = lh_hydromet_mc_m2(:,i) + lh_hydromet_mc(:,i)
        end where
      end do

      where ( X_u_all_levs(1:nnzp,sample,d_variables+1) < fraction_1 )
        lh_rcm_mc_m1(:) = lh_rcm_mc_m1(:) + lh_rcm_mc(:)
        lh_rvm_mc_m1(:) = lh_rvm_mc_m1(:) + lh_rvm_mc(:)
        lh_thlm_mc_m1(:) = lh_thlm_mc_m1(:) + lh_thlm_mc(:)
        n1(:) = n1(:) + 1

      else where
        lh_rcm_mc_m2(:) = lh_rcm_mc_m2(:) + lh_rcm_mc(:)
        lh_rvm_mc_m2(:) = lh_rvm_mc_m2(:) + lh_rvm_mc(:)
        lh_thlm_mc_m2(:) = lh_thlm_mc_m2(:) + lh_thlm_mc(:)
        n2(:) = n2(:) + 1

      end where

      ! Loop to get new sample
    end do ! sample = 1, n_micro_calls

    if ( l_compute_diagnostic_average ) then

      lh_hydromet(:,1:hydromet_dim) = lh_hydromet(:,1:hydromet_dim) / real( n_micro_calls )
      lh_thlm       = lh_thlm         / real( n_micro_calls )
      lh_rcm        = lh_rcm          / real( n_micro_calls )
      lh_rvm        = lh_rvm          / real( n_micro_calls )
      lh_wm         = lh_wm           / real( n_micro_calls )
      lh_cloud_frac = lh_cloud_frac   / real( n_micro_calls )

      ! Compute the variance of vertical velocity
      lh_wp2_zt = compute_variance( nnzp, n_micro_calls, w_tmp, lh_wm )

      ! Compute the variance of cloud water mixing ratio
      lh_rcp2_zt = compute_variance( nnzp, n_micro_calls, rc_tmp, lh_rcm )

      ! Compute the variance of total water
      lh_rtp2_zt = compute_variance( nnzp, n_micro_calls, rv_tmp+rc_tmp, lh_rvm+lh_rcm )

      ! Compute the variance of liquid potential temperature
      lh_thlp2_zt = compute_variance( nnzp, n_micro_calls, thl_tmp, lh_thlm )

      ! Compute the variance of rain water mixing ratio
      if ( iirrainm > 0 ) then
        lh_rrainp2_zt = compute_variance( nnzp, n_micro_calls, hydromet_tmp(:,iirrainm,:), &
                                          lh_hydromet(:,iirrainm) )
      end if

      ! Compute the variance of cloud droplet number concentration
      if ( iiNcm > 0 ) then
        lh_Ncp2_zt = compute_variance( nnzp, n_micro_calls, hydromet_tmp(:,iiNcm,:), &
                                       lh_hydromet(:,iiNcm) )
      end if

      ! Compute the variance of rain droplet number concentration
      if ( iiNrm > 0 ) then
        lh_Nrp2_zt = compute_variance( nnzp, n_micro_calls, hydromet_tmp(:,iiNrm,:), &
                                       lh_hydromet(:,iiNrm) )
      end if

    end if ! l_compute_diagnostic_average

    l_error = .false.
    do k = 1, nnzp
      if ( n1(k) == 0 .and. n2(k) == 0 ) then
        l_error = .true.
        write(0,*) 'Error:  no sample points in lh_microphys_driver, k =', k
      end if
    end do
    if ( l_error ) stop

    if ( l_cloud_weighted_averaging ) then
      ! Convert sums to averages.
      ! If we have no sample points for a certain plume, then we 
      ! estimate the plume liquid water by the other plume's value.
      do i = 1, hydromet_dim
        where ( n1 /= zero )
          lh_hydromet_vel_m1(:,i) = lh_hydromet_vel_m1(:,i) / dble( n1 )
          lh_hydromet_mc_m1(:,i) = lh_hydromet_mc_m1(:,i) / dble( n1 )
        end where
      end do

      where ( n1 /= zero )
        lh_rcm_mc_m1 = lh_rcm_mc_m1 / dble( n1 )
        lh_rvm_mc_m1 = lh_rvm_mc_m1 / dble( n1 )
        lh_thlm_mc_m1 = lh_thlm_mc_m1 / dble( n1 )
      end where

      do i = 1, hydromet_dim
        where ( n2 /= zero )
          lh_hydromet_vel_m2(:,i) = lh_hydromet_vel_m2(:,i) / dble( n2 )
          lh_hydromet_mc_m2(:,i) = lh_hydromet_mc_m2(:,i) / dble( n2 )
        end where
      end do

      where ( n2 /= zero )
        lh_rcm_mc_m2 = lh_rcm_mc_m2 / dble( n2 )
        lh_rvm_mc_m2 = lh_rvm_mc_m2 / dble( n2 )
        lh_thlm_mc_m2 = lh_thlm_mc_m2 / dble( n2 )
      end where

      ! Special cases
      do i = 1, hydromet_dim
        where ( n1 == zero )
          lh_hydromet_vel_m1(:,i) = lh_hydromet_vel_m2(:,i)
          lh_hydromet_mc_m1(:,i) = lh_hydromet_mc_m2(:,i)
        end where
      end do

      where ( n1 == zero )
        lh_rcm_mc_m1 = lh_rcm_mc_m2
        lh_rvm_mc_m1 = lh_rvm_mc_m2
        lh_thlm_mc_m1 = lh_thlm_mc_m2
      end where

      do i = 1, hydromet_dim
        where ( n2 == zero )
          lh_hydromet_vel_m2(:,i) = lh_hydromet_vel_m1(:,i)
          lh_hydromet_mc_m2(:,i) = lh_hydromet_mc_m1(:,i)
        end where
      end do

      where ( n2 == zero )
        lh_rcm_mc_m2 = lh_rcm_mc_m1
        lh_rvm_mc_m2 = lh_rvm_mc_m1
        lh_thlm_mc_m2 = lh_thlm_mc_m1
      end where

      ! Grid box average.
      forall( i = 1:hydromet_dim )
        lh_hydromet_vel(:,i) = real( a * cloud_frac1 * lh_hydromet_vel_m1(:,i) &
          + (1.d0-a) * cloud_frac2 * lh_hydromet_vel_m2(:,i) )
        lh_hydromet_mc(:,i)  = real( a * cloud_frac1 * lh_hydromet_mc_m1(:,i) &
          + (1.d0-a) * cloud_frac2 * lh_hydromet_mc_m2(:,i) )
      end forall

      lh_rcm_mc = real( a * cloud_frac1 * lh_rcm_mc_m1 + (1.d0-a) * cloud_frac2 * lh_rcm_mc_m2 )
      lh_rvm_mc = real( a * cloud_frac1 * lh_rvm_mc_m1 + (1.d0-a) * cloud_frac2 * lh_rvm_mc_m2 )
      lh_thlm_mc = real( a * cloud_frac1 * lh_thlm_mc_m1 + (1.d0-a) * cloud_frac2 * lh_thlm_mc_m2 )

    else

      ! Grid box average.
      forall( i = 1:hydromet_dim )
        lh_hydromet_vel(:,i) = real( lh_hydromet_vel_m1(:,i) + lh_hydromet_vel_m2(:,i) ) &
           / real( n_micro_calls )
        lh_hydromet_mc(:,i) = real( lh_hydromet_mc_m1(:,i) + lh_hydromet_mc_m2(:,i) ) &
           / real( n_micro_calls )
      end forall
      lh_rcm_mc = real( lh_rcm_mc_m1 + lh_rcm_mc_m2 ) / real( n_micro_calls )
      lh_rvm_mc = real( lh_rvm_mc_m1 + lh_rvm_mc_m2 ) / real( n_micro_calls )
      lh_thlm_mc = real( lh_thlm_mc_m1 + lh_thlm_mc_m2 ) / real( n_micro_calls )

    end if ! l_cloud_weighted_averaging

    return
  end subroutine lh_microphys_driver

!----------------------------------------------------------------------
  subroutine rc_estimate( n_micro_calls, d_variables, a, C1, C2, rc, & ! w,   & 
                         !N_pts, rr, 
                           X_u_one_lev, rc_m )
! Description:
!   Compute an Monte Carlo estimate of grid box avg liquid water.
! References:
!   None
!---------------------------------------------------------------

    use constants, only:  &
        fstderr  ! Constant(s)

!   use error_code, only:  &
!       clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_cloud_weighted_averaging   = .false.

    ! Input Variables
    integer, intent(in) :: &
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    double precision, intent(in) :: &
      a,      & ! Mixture fraction of Gaussians
      C1, C2    ! Cloud fraction associated w/ 1st, 2nd mixture component

    double precision, dimension(n_micro_calls), intent(in) :: &
      rc !, & ! n in-cloud values of spec liq water content [kg/kg].
!     w,  & ! n in-cloud values of vertical velocity (m/s)
!     Npts, & ! n in-cloud values of droplet number (#/kg air)
!     rr    ! n in-cloud values of specific rain content (kg/kg)

    double precision, dimension(n_micro_calls,d_variables+1), intent(in) :: &
      X_u_one_lev ! N x D+1 Latin hypercube sample from uniform dist

    ! Output Variables

    ! A scalar representing grid box avg specific liquid water;
    ! divide by total cloud fraction to obtain within-cloud liquid water
    double precision, intent(out) :: rc_m

    ! Local Variables

    integer :: sample
    integer :: n1, n2
    double precision :: rc_m1, rc_m2
    double precision :: coeff, expn
    double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of a, C1, C2.
    if ( a > 1.0d0 .or. a < 0.0d0 ) then
      write(fstderr,*) 'Error in rc_estimate:  ',  &
                       'mixture fraction, a, does not lie in [0,1].'
      stop
    end if
    if ( C1 > 1.0d0 .or. C1 < 0.0d0 ) then
      write(fstderr,*) 'Error in rc_estimate:  ',  &
                       'cloud fraction 1, C1, does not lie in [0,1].'
      stop
    end if
    if ( C2 > 1.0d0 .or. C2 < 0.0d0 ) then
      write(fstderr,*) 'Error in rc_estimate:  ',  &
                       'cloud fraction 2, C2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    ! Disable this for now, so we can loop over the whole domain.
    ! -dschanen 3 June 2009
!   if ( a*C1 < 0.001d0 .and. (1-a)*C2 < 0.001d0 ) then
!     if ( clubb_at_least_debug_level( 1 ) ) then
!       write(fstderr,*) 'Error in rc_estimate:  ',  &
!                        'there is no cloud or almost no cloud!'
!     end if
!   end if

    ! To compute liquid water, need to set coeff=expn=1.
    coeff = 1.d0
    expn  = 1.d0

    ! Initialize liquid in each mixture component
    rc_m1 = 0.d0
    rc_m2 = 0.d0

    ! Initialize numbers of sample points corresponding
    !    to each mixture component
    n1    = 0
    n2    = 0

    do sample = 1, n_micro_calls

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
      fraction_1 = a*C1/max( (a*C1+(1-a)*C2), epsilon( a ) )
      if ( X_u_one_lev(sample,d_variables+1) < fraction_1 ) then
        ! Use an idealized formula to compute liquid
        !      in mixture comp. 1
        rc_m1 = rc_m1 + coeff*(rc(sample))**expn
        n1    = n1 + 1
      else
        ! Use an idealized formula to compute liquid
        !      in mixture comp. 2
        rc_m2 = rc_m2 + coeff*(rc(sample))**expn
        n2    = n2 + 1
      end if

      ! Loop to get new sample
    end do

!! Convert sums to averages.
!! Old code that underestimates if n1 or n2 = 0.
!   if ( n1 == 0 ) then
!     rc_m1 = 0.d0
!   else
!     rc_m1 = rc_m1/n1
!   end if

!   if ( n2 == 0 ) then
!     rc_m2 = 0.d0
!   else
!     rc_m2 = rc_m2/n2
!   end if


    ! Convert sums to averages.
    ! If we have no sample points for a certain plume,
    !    then we estimate the plume liquid water by the
    !    other plume's value.
    if (n1 == 0 .and. n2 == 0) then
      stop 'Error:  no sample points in rc_estimate'
    end if

    if ( l_cloud_weighted_averaging ) then
      if ( .not. (n1 == 0) ) then
        rc_m1 = rc_m1/n1
      end if

     if ( .not. (n2 == 0) ) then
       rc_m2 = rc_m2/n2
     end if

      if (n1 == 0) then
        rc_m1 = rc_m2
      end if

      if (n2 == 0) then
        rc_m2 = rc_m1
      end if
      ! Grid box average.
      rc_m = a*C1*rc_m1 + (1.-a)*C2*rc_m2

    end if ! l_cloud_weighted_averaging

    ! Grid box average.
    rc_m = ( rc_m1 + rc_m2 ) / dble( n_micro_calls )

    return
  end subroutine rc_estimate
!---------------------------------------------------------------

end module estimate_lh_micro_mod
