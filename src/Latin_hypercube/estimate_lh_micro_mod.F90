!$Id$

module estimate_lh_micro_mod

  implicit none

  public :: estimate_lh_micro

  private :: autoconv_driver, rc_estimate, micro_driver

  private ! Default Scope

  contains

!------------------------------------------------------------------------

  subroutine estimate_lh_micro &
             ( dt, nnzp, n_micro_calls, d_variables, &
               X_u_all_levs, X_nl_all_levs, & 
               rt, thl, l_sample_flag, pdf_params, & 
               thlm, p_in_Pa, exner, rho, &
               wm, w_std_dev, dzq, rcm, rvm, &        
               cf, hydromet, &
               hydromet_mc_est, hydromet_vel_est, &
               rcm_mc_est, rvm_mc_est, thlm_mc_est, &
               AKm_est, AKm, AKstd, AKstd_cld, & 
               AKm_rcm, AKm_rcc, rcm_est, &
               lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &
               lh_wm, lh_wp2_zt, lh_cf, &
               microphys_sub )
! Description:
!   This subroutine computes microphysical grid box averages,
!   given a Latin Hypercube sample.
! References:
!   None
!------------------------------------------------------------------------

    use constants, only:  &
        pi,  & ! Variables(s)
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
      cf,         & ! Cloud fraction           [%]
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
      hydromet_mc_est, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_est   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      rcm_mc_est, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc_est, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc_est   ! LH estimate of time tendency of liquid potential temperature [K/s]

    real, dimension(nnzp), intent(out) :: &
      AKm_est,   & ! Monte Carlo estimate of Kessler autoconversion [kg/kg/s]
      AKm,       & ! Exact Kessler autoconversion, AKm,             [kg/kg/s]
      AKstd,     & ! Exact standard deviation of gba Kessler        [kg/kg/s]
      AKstd_cld, & ! Exact w/in cloud std of gba Kessler            [kg/kg/s]
      AKm_rcm,   & ! Exact local gba Kessler auto based on rcm      [kg/kg/s]
      AKm_rcc      ! Exact local gba Kessler based on w/in cloud rc [kg/kg/s]

    ! For comparison, estimate kth liquid water using Monte Carlo
    real, dimension(nnzp), intent(out) :: &
      rcm_est ! LH estimate of grid box avg liquid water [kg/kg]

    real, dimension(nnzp,hydromet_dim), intent(out) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp), intent(out) :: &
      lh_thlm, & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,  & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,  & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,   & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt, & ! Average value of the variance of the LH est. of vertical velocity [m^2/s^2]
      lh_cf      ! Average value of the latin hypercube est. of cloud fraction    [%]

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
    real :: R1, R2
!   real :: rc1, rc2

    ! Cloud fraction 0<cf<1, mean liquid water mix ratio [kg/kg]
!   real :: cf, rcm

    ! Double precision version of Monte Carlo Kessler ac
    double precision :: AKm_est_dp
    ! Double precision version of Monte Carlo avg liquid water
    double precision :: rcm_est_dp

    ! Variables needed for exact Kessler autoconversion, AKm
    real :: r_crit, K_one
    real :: sn1_crit, R1_crit, sn2_crit, R2_crit
    real :: AK1, AK2

    ! Variables needed for exact std of Kessler autoconversion, AKstd
    !      and within cloud standard deviation, AKstd_cld
    real :: AK1var, AK2var

    ! For comparison, compute within-cloud vertical velocity analytically.
    !real C_w_cld1, C_w_cld2, w_cld_avg

    ! ---- Begin Code ----

    ! Boundary condition
    AKm_est(1)   = 0.0
    AKm(1)       = 0.0
    AKm_rcm(1)   = 0.0
    AKm_rcc(1)   = 0.0
    rcm_est(1)   = 0.0
    AKstd(1)     = 0.0
    AKstd_cld(1) = 0.0

    do level = 2, nnzp, 1
      ! Extract PDF parameters

      !w1    = pdf_params%w1(level)
      !w2    = pdf_params%w2(level)
      !sw1   = pdf_params%sw1(level)
      !sw2   = pdf_params%sw2(level)
      !rt1   = pdf_params%rt1(level)
      !rt2   = pdf_params%rt2(level)
      !srt1  = pdf_params%srt1(level)
      !srt2  = pdf_params%srt2(level)
      !thl1  = pdf_params%thl1(level)
      !thl2  = pdf_params%thl2(level)
      !sthl1 = pdf_params%sthl1(level)
      !sthl2 = pdf_params%sthl2(level)
      a     = pdf_params%a(level)
!     rc1   = pdf_params%rc1(level)
!     rc2   = pdf_params%rc2(level)
      R1    = pdf_params%R1(level)
      R2    = pdf_params%R2(level)
      s1    = pdf_params%s1(level)
      s2    = pdf_params%s2(level)
      ss1   = pdf_params%ss1(level)
      ss2   = pdf_params%ss2(level)

      ! Compute mean cloud fraction and cloud water

!     cf    = a * R1 + (1-a) * R2
!     rcm   = a * rc1 + (1-a) * rc2

      !------------------------------------------------------------------------
      ! Call Kessler autoconversion microphysics using Latin hypercube sample
      ! This acts as an interface between the boundary layer scheme
      !    and the microphysics.  To add a call to a microphysics scheme,
      !   alter two lines in autoconversion_driver.f.
      ! Then we compute Kessler ac analytically.
      !------------------------------------------------------------------------

      ! We prognose rt-thl-w,
      !    but we set means, covariance of Nc, rr to constants.

      if ( .not. l_sample_flag(level) ) then

        ! In this case, sample points could not be constructed.
        ! Set autoconversion to zero.
        AKm_est(level)   = 0.0
        AKm(level)       = 0.0
        AKm_rcm(level)   = 0.0
        AKm_rcc(level)   = 0.0
        rcm_est(level)   = 0.0
        AKstd(level)     = 0.0
        AKstd_cld(level) = 0.0

      else

        ! Call microphysics, i.e. Kessler autoconversion.
        ! A_K = (1e-3/s)*(rc-0.5kg/kg)*H(rc-0.5kg/kg)
        call autoconv_driver &
             ( n_micro_calls, d_variables, dble( a ), dble( R1 ), dble( R2 ), &
               X_nl_all_levs(level,1:n_micro_calls,1), & 
               !X_nl(1:n,3), X_nl(1:n,4), X_nl(1:n,5),
               X_u_all_levs(level,:,:), AKm_est_dp )

        ! Convert to real number
        AKm_est(level) = real( AKm_est_dp )

        ! Compute Monte Carlo estimate of liquid for test purposes.
        call rc_estimate &
             ( n_micro_calls, d_variables, dble( a ), dble( R1 ), &
               dble( R2 ), X_nl_all_levs(level,1:n_micro_calls,1), & 
               ! X_nl(1:n,3), X_nl(1:n,4), X_nl(1:n,5),
               X_u_all_levs(level,:,:), rcm_est_dp )

        ! Convert to real number
        rcm_est(level) = real( rcm_est_dp )

        ! A test of the scheme:
        ! Compare exact (rcm) and Monte Carlo estimates (rcm_est) of
        !     specific liq water content, rcm.
        !      print*, 'rcm=', rcm
        !       print*, 'rcm_est=', rcm_est

        ! Exact Kessler autoconversion in units of (kg/kg)/s
        !        r_crit = 0.3e-3
        !        r_crit = 0.7e-3
        r_crit   = 0.2e-3
        K_one    = 1.e-3
        sn1_crit = (s1-r_crit)/ss1
        R1_crit  = 0.5*(1+erf(sn1_crit/sqrt(2.0)))
        AK1      = K_one * ( (s1-r_crit)*R1_crit  & 
                  + ss1*exp(-0.5*sn1_crit**2)/(sqrt(2*pi)) )
        sn2_crit = (s2-r_crit)/ss2
        R2_crit  = 0.5*(1+erf(sn2_crit/sqrt(2.0)))
        AK2      = K_one * ( (s2-r_crit)*R2_crit  & 
                  + ss2*exp(-0.5*sn2_crit**2)/(sqrt(2*pi)) )
        AKm(level) = a * AK1 + (1-a) * AK2

        ! Exact Kessler standard deviation in units of (kg/kg)/s
        ! For some reason, sometimes AK1var, AK2var are negative
        AK1var   = max( zero_threshold, K_one * (s1-r_crit) * AK1  & 
                 + K_one * K_one * (ss1**2) * R1_crit  & 
                 - AK1**2  )
        AK2var   = max( zero_threshold, K_one * (s2-r_crit) * AK2  & 
                 + K_one * K_one * (ss2**2) * R2_crit  & 
                 - AK2**2  )
        ! This formula is for a grid box average:
        AKstd(level)  = sqrt(    a  * ( (AK1-AKm(level))**2 + AK1var ) & 
                    + (1-a) * ( (AK2-AKm(level))**2 + AK2var ) & 
                    )
        ! This formula is for a within-cloud average:
        AKstd_cld(level) = sqrt( max( zero_threshold,   & 
                  (1./cf(level)) * ( a  * ( AK1**2 + AK1var ) & 
                            + (1-a) * ( AK2**2 + AK2var )  & 
                            ) & 
                 - (AKm(level)/cf(level))**2  ) & 
                        )

        ! Kessler autoconversion, using grid box avg liquid, rcm, as input
        AKm_rcm(level) = K_one * max( zero_threshold, rcm(level)-r_crit )

        ! Kessler ac, using within cloud liquid, rcm/cf, as input
        ! We found that for small values of cf this formula
        ! can still produce NaN values and therefore added this 
        ! threshold of 0.001 here. -dschanen 3 June 2009
        if ( cf(level) > 0.001 ) then
          AKm_rcc(level) = cf(level) * K_one * max( zero_threshold, rcm(level)/cf(level)-r_crit )
        else
          AKm_rcc(level) = zero_threshold
        end if

!       print*, 'a=', a
!       print*, 's1=', s1
!       print*, 's2=', s2
!       print*, 'ss1=', ss1
!       print*, 'ss2=', ss2
!       print*, 'AKm =', AKm(level)
!       print*, 'AKm_est (estimate) =', AKm_est(level)
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
        !     call autoconv_driver( n, d, a, R1, R2, X_nl( 1:n, 3 ), &
        !                           X_nl( 1:n, 3 ), X_nl( 1:n, 4 ), &
        !                           X_nl(1:n, 5), X_u, AKm2 )

        ! Another test:
        ! Compute within-cloud vertical velocity, avgd over full domain.
        !        C_w_cld1 =  R1*w1
        !        C_w_cld2 =  R2*w2
        !        w_cld_avg = a * C_w_cld1 + (1-a) * C_w_cld2

        ! The following two values should match
        !       print*, 'w_cld_avg=', w_cld_avg
        !       print*, 'ac_m2=', ac_m2

        ! End of overall if-then statement for Latin hypercube code
      end if

    end do ! level = 2, nnzp

    call micro_driver( dt, nnzp, n_micro_calls, d_variables, &
                       l_sample_flag, X_nl_all_levs(:,:,1), rt, thl, &
                       X_nl_all_levs(:,:,3), X_nl_all_levs(:,:,4), &
                       X_nl_all_levs(:,:,5), X_u_all_levs, &
                       thlm, p_in_Pa, exner, rho, wm, w_std_dev, &
                       dzq, rcm, rvm, pdf_params, hydromet, &
                       rvm_mc_est, rcm_mc_est, hydromet_mc_est, &
                       hydromet_vel_est, thlm_mc_est, &
                       lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &
                       lh_wm, lh_wp2_zt, lh_cf, &
                       microphys_sub )

    return
  end subroutine estimate_lh_micro
!-----------------------------------------------------------------------
  subroutine autoconv_driver( n_micro_calls, d_variables, a, R1, R2, rc, &
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

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    double precision, intent(in) :: &
      a,      & ! Mixture fraction of Gaussians
      R1, R2    ! Cloud fraction associated w/ 1st, 2nd mixture component

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

    ! Handle some possible errors re: proper ranges of a, R1, R2.
    if ( a > 1.0d0 .or. a < 0.0d0 ) then
      write(fstderr,*) 'Error in autoconv_driver:  ',  &
                       'mixture fraction, a, does not lie in [0,1].'
      write(fstderr,*) 'a = ', a
      stop
    end if
    if ( R1 > 1.0d0 .or. R1 < 0.0d0 ) then
      write(fstderr,*) 'Error in autoconv_driver:  ',  &
                       'cloud fraction 1, R1, does not lie in [0,1].'
      write(fstderr,*) 'R1 = ', R1
      stop
    end if
    if ( R2 > 1.0d0 .or. R2 < 0.0d0 ) then
      write(fstderr,*) 'Error in autoconv_driver:  ',  &
                       'cloud fraction 2, R2, does not lie in [0,1].'
      write(fstderr,*) 'R2 = ', R2
      stop
    end if

    ! Make sure there is some cloud.
    ! Disable this for now, so we can loop over the whole domain.
    ! -dschanen 3 June 2009
!   if ( a*R1 < 0.001d0 .and. (1-a)*R2 < 0.001d0 ) then
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
      fraction_1 = a*R1/max( a*R1+(1-a)*R2, epsilon( a ) )
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

! Convert sums to averages.
! If we have no sample points for a certain plume,
!    then we estimate the plume liquid water by the
!    other plume's value.
    if ( n1 == 0 .and. n2 == 0 ) then
      stop 'Error:  no sample points in autoconv_driver'
    end if

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
    ac_m = a*R1*ac_m1 + (1-a)*R2*ac_m2

!   print*, 'autoconv_driver: acm=', ac_m

    return
  end subroutine autoconv_driver

!-------------------------------------------------------------------------------
  subroutine micro_driver( dt, nnzp, n_micro_calls, d_variables, &
                           l_sample_flag, &
                           s_mellor, rt, thl, w, Nc, rr, X_u_all_levs, &
                           thlm, p_in_Pa, exner, rho, wm, w_std_dev, &
                           dzq, rcm, rvm, pdf_params, hydromet,  &
                           rvm_mc_est, rcm_mc_est, hydromet_mc_est, &
                           hydromet_vel_est, thlm_mc_est, &
                           lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &
                           lh_wm, lh_wp2_zt, lh_cf, &
                           microphys_sub )
! Description:
!   Estimate the tendency of a microphysics scheme via latin hypercube sampling
! References:
!   None
!-------------------------------------------------------------------------------

    use constants, only:  &
      fstderr  ! Constant(s)

    use parameters_model, only: &
      hydromet_dim

    use array_index, only: &
      iirrainm, &
      iiNcm

    use variables_prognostic_module, only: &
      pdf_parameter

    implicit none

    ! External
#include "../microphys_interface.inc"

    intrinsic :: epsilon

    ! Constant parameters
    logical, parameter :: &
      l_compute_diagnostic_average = .true.

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
      s_mellor,  & ! n_micro_calls values of 's' (Mellor 1977)            [kg/kg]
      rt,        & ! n_micro_calls values of total water mixing ratio     [kg/kg]
      thl,       & ! n_micro_calls values of liquid potential temperature [K]
      w,         & ! n_micro_calls values of vertical velocity            [m/s]
      Nc,        & ! n_micro_calls values of droplet number               [#/kg air]
      rr           ! n_micro_calls values of specific rain content        [kg/kg]

    double precision, dimension(nnzp,n_micro_calls,d_variables+1), intent(in) :: &
      X_u_all_levs ! N x D+1 Latin hypercube sample from uniform dist

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
      hydromet_mc_est, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_est   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      rcm_mc_est, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc_est, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc_est   ! LH estimate of time tendency of liquid potential temperature [K/s]

    real, dimension(nnzp,hydromet_dim), intent(out) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp), intent(out) :: &
      lh_thlm, & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,  & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,  & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,   & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt, & ! Average value of the variance of the LH est. of vertical velocity [m^2/s^2]
      lh_cf      ! Average value of the latin hypercube est. of cloud fraction    [%]

    ! Local Variables
    double precision, dimension(nnzp,hydromet_dim) :: &
      hydromet_mc_est_m1,  & ! LH est of hydrometeor time tendency          [(units vary)/s]
      hydromet_mc_est_m2,  & ! LH est of hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_est_m1, & ! LH est of hydrometeor sedimentation velocity [m/s]
      hydromet_vel_est_m2    ! LH est of hydrometeor sedimentation velocity [m/s]

    double precision, dimension(nnzp) :: &
      rcm_mc_est_m1,  & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      rcm_mc_est_m2,  & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc_est_m1,  & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      rvm_mc_est_m2,  & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc_est_m1, & ! LH est of time tendency of liquid potential temperature [K/s]
      thlm_mc_est_m2    ! LH est of time tendency of liquid potential temperature [K/s]


    real, dimension(nnzp,hydromet_dim) :: &
      hydromet_tmp ! Hydrometeor species    [units vary]

    real, dimension(nnzp) :: &
      rc_tmp,    & ! Liquid water                [kg/kg]
      rv_tmp,    & ! Vapor water                 [kg/kg]
      thl_tmp      ! Liquid potential temperature[K]

    real, dimension(nnzp,n_micro_calls) :: &
      w_tmp        ! Vertical velocity           [m/s]

    integer, dimension(nnzp) :: n1, n2, zero

    double precision, dimension(nnzp) :: &
      R1, R2, a, &
      fraction_1

    integer :: i, k, sample

    logical :: l_error
    ! ---- Begin Code ----

    a(:)  = dble( pdf_params%a(:) )
!   R1(:) = dble( pdf_params%R1(:) )
!   R2(:) = dble( pdf_params%R2(:) )
    R1(:) = dble( 1.0 )
    R2(:) = dble( 1.0 )

    zero(:) = 0

    ! Initialize microphysical tendencies for each mixture component
    hydromet_mc_est_m1(:,:) = 0.d0
    hydromet_mc_est_m2(:,:) = 0.d0

    hydromet_vel_est_m1(:,:) = 0.d0
    hydromet_vel_est_m2(:,:) = 0.d0

    rcm_mc_est_m1(:) = 0.d0
    rcm_mc_est_m2(:) = 0.d0

    rvm_mc_est_m1(:) = 0.d0
    rvm_mc_est_m2(:) = 0.d0

    thlm_mc_est_m1(:) = 0.d0
    thlm_mc_est_m2(:) = 0.d0

    ! Initialize numbers of sample points corresponding
    !    to each mixture component
    n1 = 0
    n2 = 0

    ! Choose which mixture fraction we are in.
    ! Account for cloud fraction.
    ! Follow M. E. Johnson (1987), p. 56.
    fraction_1(:) = a(:)*R1(:)/( a(:)*R1(:)+(1.-a(:))*R2(:) )
!   print*, 'fraction_1= ', fraction_1

    do sample = 1, n_micro_calls

      where ( l_sample_flag )
        where( s_mellor(:,sample) > 0.0 )
          rc_tmp  = real( s_mellor(:,sample) ) 
        else where
          rc_tmp = 0.0
        end where
      else where
        rc_tmp = rcm
      end where

      where ( l_sample_flag ) 
        w_tmp(:,sample) = real( w(:,sample) )
        rv_tmp  = real( rt(:,sample) ) - rc_tmp
        thl_tmp = real( thl(:,sample) )
      else where
        w_tmp(:,sample) = wm
        rv_tmp  = rvm
        thl_tmp = thlm
      end where

      do i = 1, hydromet_dim, 1
        if ( i == iirrainm ) then
          where ( l_sample_flag )
            hydromet_tmp(:,i) = real( rr(:,sample) )
          else where
            hydromet_tmp(:,i) = hydromet(:,i)
          end where
        else if ( i == iiNcm ) then
          where ( l_sample_flag )
            hydromet_tmp(:,i) = real( Nc(:,sample) )
          else where
            hydromet_tmp(:,i) = hydromet(:,i)
          end where
        else
          hydromet_tmp(:,i) = hydromet(:,i)
        end if
      end do

      if ( l_compute_diagnostic_average ) then
        if ( sample == 1 ) then
          lh_hydromet(:,1:hydromet_dim) = hydromet_tmp(:,1:hydromet_dim)
          lh_thlm = thl_tmp
          lh_rcm  = rc_tmp
          lh_rvm  = rv_tmp
          lh_wm   = w_tmp(:,sample)
          where ( l_sample_flag .and. s_mellor(:,sample) > 0 ) 
            lh_cf = 1.0
          else where
            lh_cf = 0.0
          end where
        else
          lh_hydromet(:,1:hydromet_dim) = lh_hydromet(:,1:hydromet_dim) &
            + hydromet_tmp(:,1:hydromet_dim)
          lh_thlm = lh_thlm + thl_tmp
          lh_rcm  = lh_rcm  + rc_tmp
          lh_rvm  = lh_rvm  + rv_tmp
          lh_wm   = lh_wm   + w_tmp(:,sample)
          where ( l_sample_flag .and. s_mellor(:,sample) > 0 ) lh_cf = lh_cf + 1.0
        end if
      end if

      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub &
           ( dt, nnzp, .false., .true., thl_tmp, p_in_Pa, exner, rho, pdf_params, &
             w_tmp(:,sample), w_std_dev, dzq, rc_tmp, rv_tmp, hydromet_tmp, hydromet_mc_est, &
             hydromet_vel_est, rcm_mc_est, rvm_mc_est, thlm_mc_est )

      do i = 1, hydromet_dim
        where ( X_u_all_levs(1:nnzp,sample,d_variables+1) < fraction_1 )
          hydromet_vel_est_m1(:,i) = hydromet_vel_est_m1(:,i) + hydromet_vel_est(:,i)
          hydromet_mc_est_m1(:,i) = hydromet_mc_est_m1(:,i) + hydromet_mc_est(:,i)
        else where
          hydromet_vel_est_m2(:,i) = hydromet_vel_est_m2(:,i) + hydromet_vel_est(:,i)
          hydromet_mc_est_m2(:,i) = hydromet_mc_est_m2(:,i) + hydromet_mc_est(:,i)
        end where
      end do

      where ( X_u_all_levs(1:nnzp,sample,d_variables+1) < fraction_1 )
        rcm_mc_est_m1(:) = rcm_mc_est_m1(:) + rcm_mc_est(:)
        rvm_mc_est_m1(:) = rvm_mc_est_m1(:) + rvm_mc_est(:)
        thlm_mc_est_m1(:) = thlm_mc_est_m1(:) + thlm_mc_est(:)
        n1(:) = n1(:) + 1

      else where
        rcm_mc_est_m2(:) = rcm_mc_est_m2(:) + rcm_mc_est(:)
        rvm_mc_est_m2(:) = rvm_mc_est_m2(:) + rvm_mc_est(:)
        thlm_mc_est_m2(:) = thlm_mc_est_m2(:) + thlm_mc_est(:)
        n2(:) = n2(:) + 1

      end where

      ! Loop to get new sample
    end do ! sample = 1, n_micro_calls

    if ( l_compute_diagnostic_average ) then

      lh_hydromet(:,1:hydromet_dim) = lh_hydromet(:,1:hydromet_dim) / real( n_micro_calls )
      lh_thlm = lh_thlm / real( n_micro_calls )
      lh_rcm  = lh_rcm  / real( n_micro_calls )
      lh_rvm  = lh_rvm  / real( n_micro_calls )
      lh_wm   = lh_wm   / real( n_micro_calls )
      lh_cf   = lh_cf   / real( n_micro_calls )

      ! Compute the variance of vertical velocity
      lh_wp2_zt = 0.0
      do sample = 1, n_micro_calls
        lh_wp2_zt = lh_wp2_zt + ( w_tmp(:,sample) - lh_wm )**2
      end do
      lh_wp2_zt = lh_wp2_zt / real( n_micro_calls )

    end if

!   rvm_tmp = 0
!   do sample = 1, n_micro_calls
!   do k = 1, nnzp, 1
!     if ( l_sample_flag(k) ) then
!       rvm_tmp(k) = rvm_tmp(k) + real( rt(k,sample) ) /1000.
!     else
!       rvm_tmp(k) = rvm_tmp(k) + rvm(k) + rcm(k)
!     end if
!   end do
!   end do
!   do k = 1, nnzp
!     write(6,'(2G12.4)') rvm_tmp(k) / real( n_micro_calls ), rvm(k)+rcm(k)
!   end do
!   pause

!   thlm_tmp = 0
!   do sample = 1, n_micro_calls
!   do k = 1, nnzp, 1
!     if ( l_sample_flag(k) ) then
!       thlm_tmp(k) = thlm_tmp(k) + real( thl(k,sample) )
!     else
!       thlm_tmp(k) = thlm_tmp(k) + thlm(k)
!     end if
!   end do
!   end do
!   do k = 1, nnzp
!     write(6,'(2G12.4)') thlm_tmp(k) / real( n_micro_calls ), thlm(k)
!   end do
!   pause

! Convert sums to averages.
! If we have no sample points for a certain plume,
!    then we estimate the plume liquid water by the
!    other plume's value.
    l_error = .false.
    do k = 1, nnzp
      if ( n1(k) == 0 .and. n2(k) == 0 ) then
        l_error = .true.
        write(0,*) 'Error:  no sample points in micro_driver, k =', k
      end if
    end do
    if ( l_error ) stop

    do i = 1, hydromet_dim
      where ( n1 /= zero )
        hydromet_vel_est_m1(:,i) = hydromet_vel_est_m1(:,i) / dble( n1 )
        hydromet_mc_est_m1(:,i) = hydromet_mc_est_m1(:,i) / dble( n1 )
      end where
    end do

    where ( n1 /= zero )
      rcm_mc_est_m1 = rcm_mc_est_m1 / dble( n1 )
      rvm_mc_est_m1 = rvm_mc_est_m1 / dble( n1 )
      thlm_mc_est_m1 = thlm_mc_est_m1 / dble( n1 )
    end where

    do i = 1, hydromet_dim
      where ( n2 /= zero )
        hydromet_vel_est_m2(:,i) = hydromet_vel_est_m2(:,i) / dble( n2 )
        hydromet_mc_est_m2(:,i) = hydromet_mc_est_m2(:,i) / dble( n2 )
      end where
    end do

    where ( n2 /= zero )
      rcm_mc_est_m2 = rcm_mc_est_m2 / dble( n2 )
      rvm_mc_est_m2 = rvm_mc_est_m2 / dble( n2 )
      thlm_mc_est_m2 = thlm_mc_est_m2 / dble( n2 )
    end where

    ! Special cases
    do i = 1, hydromet_dim
      where ( n1 == zero )
        hydromet_vel_est_m1(:,i) = hydromet_vel_est_m2(:,i)
        hydromet_mc_est_m1(:,i) = hydromet_mc_est_m2(:,i)
      end where
    end do

    where ( n1 == zero )
      rcm_mc_est_m1 = rcm_mc_est_m2
      rvm_mc_est_m1 = rvm_mc_est_m2
      thlm_mc_est_m1 = thlm_mc_est_m2
    end where

    do i = 1, hydromet_dim
      where ( n2 == zero )
        hydromet_vel_est_m2(:,i) = hydromet_vel_est_m1(:,i)
        hydromet_mc_est_m2(:,i) = hydromet_mc_est_m1(:,i)
      end where
    end do

    where ( n2 == zero )
      rcm_mc_est_m2 = rcm_mc_est_m1
      rvm_mc_est_m2 = rvm_mc_est_m1
      thlm_mc_est_m2 = thlm_mc_est_m1
    end where

    ! Grid box average.
    forall( i = 1:hydromet_dim )
      hydromet_vel_est(:,i) = real( a * R1 * hydromet_vel_est_m1(:,i) &
        + (1.d0-a) * R2 * hydromet_vel_est_m2(:,i) )
      hydromet_mc_est(:,i)  = real( a * R1 * hydromet_mc_est_m1(:,i) &
        + (1.d0-a) * R2 * hydromet_mc_est_m2(:,i) )
    end forall

    rcm_mc_est = real( a * R1 * rcm_mc_est_m1 + (1.d0-a) * R2 * rcm_mc_est_m2 )
    rvm_mc_est = real( a * R1 * rvm_mc_est_m1 + (1.d0-a) * R2 * rvm_mc_est_m2 )
    thlm_mc_est = real( a * R1 * thlm_mc_est_m1 + (1.d0-a) * R2 * thlm_mc_est_m2 )

    return
  end subroutine micro_driver

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

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    double precision, intent(in) :: &
      a,      & ! Mixture fraction of Gaussians
      C1, C2    ! Cloud fraction associated w/ 1st, 2nd mixture component

    double precision, dimension(n_micro_calls), intent(in) :: &
      rc !, & ! n in-cloud values of spec liq water content [g/kg].
!     w,  & ! n in-cloud values of vertical velocity (m/s)
!     Npts, & ! n in-cloud values of droplet number (#/mg air)
!     rr    ! n in-cloud values of specific rain content (g/kg)

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
      fraction_1 = a*C1/(a*C1+(1-a)*C2)
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
    rc_m = a*C1*rc_m1 + (1-a)*C2*rc_m2

    return
  end subroutine rc_estimate
!---------------------------------------------------------------

end module estimate_lh_micro_mod
