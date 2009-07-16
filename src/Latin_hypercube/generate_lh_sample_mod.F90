!$Id$
module generate_lh_sample_mod

  implicit none

  public :: generate_lh_sample

  private :: sample_points, latin_hyper_sample, gaus_mixt_points, & 
             truncate_gaus_mixt, ltqnorm, gaus_condt, & 
             st_2_rtthl, log_sqd_normalized


  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine generate_lh_sample &
             ( n_micro_calls, nt_repeat, d_variables, hydromet_dim, & 
               p_matrix, cf, pdf_params, level, & 
               hydromet, correlation_matrix, &
               rt, thl, & 
               X_u_one_lev, X_nl_one_lev, l_sample_flag )
! Description:
!   This subroutine generates a Latin Hypercube sample.

! References:
!   ``Supplying Local Microphysical Parameterizations with Information about
!     Subgrid Variability: Latin Hypercube Sampling'', JAS Vol. 62,
!     p. 4010--4026, Larson, et al. 2005.
!-------------------------------------------------------------------------------

    use KK_microphys_module, only: &
      PDF_TRIVAR_2G_LN_LN

    use constants, only:  &
      max_mag_correlation, &  ! Constant
      g_per_kg,  &  ! g/kg
      cm3_per_m3, & ! cm3 per m3
      rr_tol, &     ! rr tolerance in kg/kg
      Nr_tol, &     ! Nr tolerance in #/kg
      Nc_tol        ! Nc tolerance in #/kg

    use variables_prognostic_module, only:  &
      pdf_parameter  ! type

    use array_index, only: &
      iiNcm,    & ! Variables
      iiNrm,    &
      iirrainm, &
      iiLH_rrain, &
      iiLH_Nr, &
      iiLH_Nc, &
      iiLH_rt, &
      iiLH_thl, &
      iiLH_w

    implicit none

    ! External
    intrinsic :: dble, min, max, sqrt

    ! Constant Parameters
    logical, parameter :: &
      l_sample_out_of_cloud = .true.

    ! Input Variables
    integer, intent(in) :: &
      n_micro_calls, & ! `n'   Number of calls to microphysics (normally=2)
      nt_repeat,     & ! `n_t' Num. random samples before sequence repeats (normally=10)
      d_variables,   & ! `d'   Number of variates (normally=5)
      hydromet_dim     ! Number of hydrometeor species

    real, dimension(hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species [units vary]

    ! Cloud fraction
    real, intent(in) :: cf !  Cloud fraction, 0 <= cf <= 1

    integer, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      p_matrix !   N x D random matrix of integers that's fed to closure

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters output by closure_new [units vary]

    integer, intent(in) :: level  ! Level info. for PDF parameters.

    ! From the KK_microphys_module
    real, dimension(d_variables,d_variables), intent(in) :: &
      correlation_matrix ! Correlations

    ! Output Variables
    double precision, intent(out), dimension(n_micro_calls) :: &
      rt, & ! Total water mixing ratio          [kg/kg]
      thl   ! Liquid potential temperature      [K]

    double precision, intent(out), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev ! Sample drawn from uniform distribution from a particular grid level

    double precision, intent(out), dimension(n_micro_calls,d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal

    ! A true/false flag that determines whether
    ! the PDF allows us to construct a sample
    logical, intent(out) :: l_sample_flag

    ! Local Variables

    logical, dimension(d_variables) :: &
      l_d_variable_lognormal

    real :: a
    real :: w1, w2, sw1, sw2
    real :: thl1, thl2
    real :: sthl1, sthl2
    real :: rt1, rt2
    real :: srt1, srt2
    ! sub-plume correlation coefficient between rt, thl
    ! varies between -1 < rrtthl < 1
    real :: rrtthl

    real :: s1, s2, ss1, ss2
    real :: cloud_frac1, cloud_frac2
    real :: crt1, crt2
    real :: cthl1, cthl2

    ! Clip the magnitude of the correlation between rt and thl
    real :: rrtthl_reduced
    double precision :: rrtthl_reduced1, rrtthl_reduced2

    ! Means of s, t, w, Nc, Nr, rr for plumes 1 and 2
    double precision, dimension(d_variables) :: &
      mu1, mu2

    ! Covariance (not correlation) matrix of rt, thl, w, Nc, rr
    ! for plumes 1 and 2
    double precision, dimension(d_variables,d_variables) :: &
      Sigma_rtthlw_1,  & 
      Sigma_rtthlw_2

    double precision :: &
      Ncm,  & ! Cloud droplet number concentration.[number / kg air]
      Nc1,  & ! PDF parameter for mean of plume 1. [#/kg]
      Nc2,  & ! PDF parameter for mean of plume 2. [#/kg]
      sNc1, & ! PDF param for width of plume 1.    [(#/kg)^2]
      sNc2, & ! PDF param for width of plume 2.    [(#/kg^2]
      Nrm,  & ! Rain droplet number concentration. [number / kg air]
      Nr1,  & ! PDF parameter for mean of plume 1. [#/kg]
      Nr2,  & ! PDF parameter for mean of plume 2. [#/kg]
      sNr1, & ! PDF param for width of plume 1.    [(#/kg)^2]
      sNr2    ! PDF param for width of plume 2.    [(#/kg^2]

    real :: corr_srr, corr_sNr, corr_rrNr, srrNr1, srrNr2

    ! rr = specific rain content. [rr] = g rain / kg air
    double precision :: &
      rrainm, &  ! rain water mixing ratio            [g/kg]
      rr1, &  ! PDF parameter for mean of plume 1. [kg/kg]
      rr2, &  ! PDF parameter for mean of plume 2. [kg/kg]
      srr1, & ! PDF param for width of plume 1     [(kg/kg)^2]
      srr2    ! PDF param for width of plume 2.    [(kg/kg)^2]

    double precision :: &
      Ncp2_on_Ncm2, & ! = Ncp2 divided by Ncm^2    [-]
      Nrp2_on_Nrm2, & ! = Nrp2 divided by Nrm^2    [-]
      rrp2_on_rrainm2 ! = rrp2 divided by rrainm^2 [-]

    logical, dimension(hydromet_dim) :: l_small_lognormal

    integer :: i

    ! ---- Begin Code ----
    ! Determine which variables are a lognormal distribution
    i = max( iiLH_rt, iiLH_thl, iiLH_w )
    l_d_variable_lognormal(1:i) = .false. ! The 1st 3 variates
    l_d_variable_lognormal(i+1:d_variables) = .true.  ! Hydrometeors

    l_small_lognormal(:) = .false.

    ! Input pdf parameters.

    w1          = pdf_params%w1(level)
    w2          = pdf_params%w2(level)
    sw1         = pdf_params%sw1(level)
    sw2         = pdf_params%sw2(level)
    rt1         = pdf_params%rt1(level)
    rt2         = pdf_params%rt2(level)
    srt1        = pdf_params%srt1(level)
    srt2        = pdf_params%srt2(level)
    thl1        = pdf_params%thl1(level)
    thl2        = pdf_params%thl2(level)
    sthl1       = pdf_params%sthl1(level)
    sthl2       = pdf_params%sthl2(level)
    a           = pdf_params%a(level)
    cloud_frac1 = pdf_params%cloud_frac1(level)
    cloud_frac2 = pdf_params%cloud_frac2(level)
    s1          = pdf_params%s1(level)
    s2          = pdf_params%s2(level)
    ss1         = pdf_params%ss1(level)
    ss2         = pdf_params%ss2(level)
    rrtthl      = pdf_params%rrtthl(level)
    crt1        = pdf_params%crt1(level)
    crt2        = pdf_params%crt2(level)
    cthl1       = pdf_params%cthl1(level)
    cthl2       = pdf_params%cthl2(level)

    !-----------------------------------------------------------------------
    !
    ! Call Latin Hypercube sampler to compute microphysics.
    !    V. Larson Mar 2004.
    ! This acts as an interface between the boundary layer scheme
    !   and the microphysics.  To add a call to a microphysics scheme,
    !   alter two lines in autoconversion_driver.f.
    !-----------------------------------------------------------------------

    ! We prognose rt-thl-w,
    !    but we set means, covariance of Nc, qr to constants.

    l_sample_flag = .true.
    if ( l_sample_out_of_cloud ) then
      ! Sample non-cloudy grid boxes as well -dschanen 3 June 2009
      cloud_frac1 = 1.0
      cloud_frac2 = 1.0

    else if ( cf < 0.001 ) then
      ! In this case there are essentially no cloudy points to sample;
      ! Set sample points to zero.
      X_u_one_lev(:,:)    = 0.0
      X_nl_one_lev(:,:)   = 0.0
      l_sample_flag = .false.

    end if ! l_sample_out_of_cloud

    if ( srt1  == 0. .or. srt2  == 0. .or. & 
         sthl1 == 0. .or. sthl2 == 0. .or. & 
         sw1   == 0. .or. sw2   == 0. ) then

      ! In this case, Sigma_rtthlw matrix is ill-conditioned;
      !     then matrix operations will fail.

!         print*,'srt1=', srt1
!         print*,'srt2=', srt2
!         print*,'sthl1=', sthl1
!         print*,'sthl2=', sthl2
!         print*,'sw1=', sw1
!         print*,'sw2=', sw2

!         print*, 'Covariance matrix of r-thl-w is ill-conditioned'

      X_u_one_lev(:,:)    = 0.0
      X_nl_one_lev(:,:)   = 0.0
      l_sample_flag = .false.

    end if

    if ( l_sample_flag ) then


      ! Compute PDF parameters for Nc, rr.
      ! Assume that Nc, rr obey single-lognormal distributions

      ! Nc  = droplet number concentration.  [Nc] = number / kg air
      ! Ncm  = mean of N; must have Ncm>0
      ! Ncp2_on_Ncm2 = variance of Nc divided by Ncm^2; must have Ncp2>0.
      ! Nc1  = PDF parameter for mean of plume 1. [Nc1] = (#/kg)
      ! Nc2  = PDF parameter for mean of plume 2. [Nc2] = (#/kg)
      ! sNc1,2 = PDF param for width of plume 1,2. [sNc1,2] = (#/kg)**2

      if ( iiLH_Nc > 0 ) then 
        Ncm = dble( hydromet(iiNcm) )
        Ncp2_on_Ncm2 = dble( correlation_matrix(iiLH_Nc,iiLH_Nc) )

        call log_sqd_normalized( Ncm, Ncp2_on_Ncm2, dble( Nc_tol ), &
                                 Nc1, Nc2, sNc1, sNc2, l_small_lognormal(iiNcm) )
      end if

      ! rr = specific rain content. [rr] = g rain / kg air
      ! rrainm  = mean of rr; rrp2 = variance of rr, must have rrp2>0.
      ! rr1  = PDF parameter for mean of plume 1. [rr1] = (kg/kg)
      ! rr2  = PDF parameter for mean of plume 2. [rr2] = (kg/kg)
      ! srr1,2 = PDF param for width of plume 1,2. [srr1,2] = (kg/kg)**2

      if ( iiLH_rrain > 0 ) then 
        rrainm = dble( hydromet(iirrainm) )
        rrp2_on_rrainm2 = dble( correlation_matrix(iiLH_rrain,iiLH_rrain) )
        call log_sqd_normalized( rrainm, rrp2_on_rrainm2, dble( rr_tol ), &
                                 rr1, rr2, srr1, srr2, l_small_lognormal(iirrainm) )
      end if

      if ( iiLH_Nr > 0 ) then 
        Nrm = dble( hydromet(iiNrm) )
        Nrp2_on_Nrm2 = dble( correlation_matrix(iiLH_Nr,iiLH_Nr) )

        call log_sqd_normalized( Nrm, Nrp2_on_Nrm2, dble( Nr_tol ), &
                                 Nr1, Nr2, sNr1, sNr2, l_small_lognormal(iiNcm) )
      end if

      ! Means of s, t, w, Nc, Nr, rr for Gaussians 1 and 2

      mu1((/iiLH_rt,iiLH_thl,iiLH_w/)) &
        = (/ dble( s1 ), 0.d0, dble( w1 ) /)
      mu2((/iiLH_rt,iiLH_thl,iiLH_w/)) &
        = (/ dble( s2 ), 0.d0, dble( w2 ) /)

      if ( iiLH_rrain > 0 ) then
        mu1(iiLH_rrain) = rr1
        mu2(iiLH_rrain) = rr2
      end if
      if ( iiLH_Nc > 0 ) then
        mu1(iiLH_Nc) = Nc1
        mu2(iiLH_Nc) = Nc2
      end if

      if ( iiLH_Nr > 0 ) then
        mu1(iiLH_Nr) = Nr1
        mu2(iiLH_Nr) = Nr2
      end if

      ! An old subroutine, gaus_rotate, couldn't handle large correlations;
      !   I assume the replacement, gaus_condt, has equal trouble.
      !   Therefore we input smaller correlations
      ! max_mag_correlation = 0.99 in constants.F90
      rrtthl_reduced = min( max_mag_correlation, max( rrtthl, -max_mag_correlation ) )

      ! Within-plume rt-thl correlation terms with rt in kg/kg
      rrtthl_reduced1 = dble(rrtthl_reduced*sqrt( srt1*sthl1 ))
      rrtthl_reduced2 = dble(rrtthl_reduced*sqrt( srt2*sthl2 ))

      ! Covariance (not correlation) matrices of rt-thl-w-Nc-rr
      !    for Gaussians 1 and 2
      ! For now, assume no within-plume correlation of w,Nc,Nr,rr with
      !    any other variables.

      ! Sigma_rtthlw_1,2
      Sigma_rtthlw_1 = 0.d0 ! Start with no correlation, and add matrix elements
      Sigma_rtthlw_2 = 0.d0

      Sigma_rtthlw_1(iiLH_rt,(/iiLH_rt,iiLH_thl/))  = (/ dble( srt1 ), rrtthl_reduced1 /)

      Sigma_rtthlw_2(iiLH_rt,(/iiLH_rt,iiLH_thl/))  = (/ dble( srt2 ), rrtthl_reduced2 /)

      Sigma_rtthlw_1(iiLH_thl,(/iiLH_rt,iiLH_thl/)) = (/ rrtthl_reduced1, dble( sthl1 ) /)

      Sigma_rtthlw_2(iiLH_thl,(/iiLH_rt,iiLH_thl/)) = (/ rrtthl_reduced2, dble( sthl2 ) /)

      Sigma_rtthlw_1(iiLH_w,iiLH_w) = dble( sw1 )

      Sigma_rtthlw_2(iiLH_w,iiLH_w) = dble( sw2 )

      if ( iiLH_Nc > 0 ) then
        Sigma_rtthlw_1(iiLH_Nc,iiLH_Nc) = sNc1
        Sigma_rtthlw_2(iiLH_Nc,iiLH_Nc) = sNc2
      end if

      if ( iiLH_Nr > 0 ) then
        Sigma_rtthlw_1(iiLH_Nr,iiLH_Nr) = sNr1
        Sigma_rtthlw_2(iiLH_Nr,iiLH_Nr) = sNr2
      end if

      if ( iiLH_rrain > 0 ) then
        Sigma_rtthlw_1(iiLH_rrain,iiLH_rrain) = srr1
        Sigma_rtthlw_2(iiLH_rrain,iiLH_rrain) = srr2
      end if

      if ( iiLH_rrain > 0 .and. iiLH_Nr > 0 ) then 
        if ( rrainm > dble( rr_tol ) .and. Nrm > dble( Nr_tol ) ) then
          corr_srr = correlation_matrix(iiLH_rt,iiLH_rrain)
          corr_sNr = correlation_matrix(iiLH_rt,iiLH_rrain)
          corr_rrNr = correlation_matrix(iiLH_rrain,iiLH_Nr)
          srrNr1 = PDF_TRIVAR_2G_LN_LN( s1, real( rrainm ) , real( Nrm ), &
                  ss1, real( rrainm * sqrt( rrp2_on_rrainm2 ) ), &
                  real( Nrm * sqrt( Nrp2_on_Nrm2 ) ), &
                  corr_srr, corr_sNr, corr_rrNr, 1.0, 1./3., 2./3. )
          srrNr2 = PDF_TRIVAR_2G_LN_LN( s2, real( rrainm ), real( Nrm ), &
                  ss2, real( rrainm * sqrt( rrp2_on_rrainm2 ) ), &
                  real( Nrm * sqrt( Nrp2_on_Nrm2 ) ), &
                  corr_srr, corr_sNr, corr_rrNr, 1.0, 1./3., 2./3. )
          Sigma_rtthlw_1(iiLH_rrain,iiLH_Nr) = dble( srrNr1 )
          Sigma_rtthlw_1(iiLH_Nr,iiLH_rrain) = dble( srrNr1 )
          Sigma_rtthlw_2(iiLH_rrain,iiLH_Nr) = dble( srrNr2 )
          Sigma_rtthlw_2(iiLH_Nr,iiLH_rrain) = dble( srrNr2 )
        end if
      end if

      call sample_points( n_micro_calls, nt_repeat, d_variables, p_matrix, dble( a ), & 
                          dble( rt1 ), dble( thl1 ),  & 
                          dble( rt2 ), dble( thl2 ), & 
                          dble( crt1 ), dble( cthl1 ),  & 
                          dble( crt2 ), dble( cthl2 ), & 
                          dble( mu1 ), dble( mu2 ),  & 
                          Sigma_rtthlw_1, Sigma_rtthlw_2, & 
                          dble( cloud_frac1 ), dble( cloud_frac2 ), & 
                          l_d_variable_lognormal, &
                          rt, thl, X_u_one_lev, X_nl_one_lev ) 

      ! Kluge for lognormal variables
      ! Use the gridbox mean rather than a sampled value
      if ( iiLH_rrain > 0 .and. l_small_lognormal(iirrainm) ) then
        X_nl_one_lev(:,iiLH_rrain) = rrainm
      end if

      if ( iiLH_Nc > 0 .and. l_small_lognormal(iiNcm) ) then
        X_nl_one_lev(:,iiLH_Nc) = Ncm
      end if

      if ( iiLH_Nr > 0 .and. l_small_lognormal(iiNrm) ) then
        X_nl_one_lev(:,iiLH_Nr) = Nrm
      end if
      ! End of overall if-then statement for Latin hypercube code
    end if

    return
  end subroutine generate_lh_sample

!----------------------------------------------------------------------
! Description:
!   Generates n random samples from a d-dim Gaussian-mixture PDF.
!   Uses Latin hypercube method.
!   To be called from pdf_closure of CLUBB.

!   We take samples only from the cloudy part of the grid box.
!   We use units of kg/kg.
! References:
!   None
!----------------------------------------------------------------------

  subroutine sample_points( n_micro_calls, nt_repeat, d_variables, p_matrix, a, & 
                            rt1, thl1, rt2, thl2, & 
                            crt1, cthl1, crt2, cthl2, & 
                            mu1, mu2,  & 
                            Sigma_rtthlw_1, Sigma_rtthlw_2, & 
                            cloud_frac1, cloud_frac2, & 
                            l_d_variable_lognormal, &
                            rt, thl, X_u_one_lev, X_nl_one_lev )

    use constants, only:  &
        fstderr  ! Constant(s)

    use array_index, only: &
      iiLH_rt, & ! Variable(s)
      iiLH_thl

    implicit none

    ! Input variables

    integer, intent(in) :: &
      n_micro_calls, & ! `n'   Number of calls to microphysics (normally=2)
      nt_repeat,     & ! `n_t' Num. random samples before sequence repeats (normally=10)
      d_variables      ! Number of variates (normally=5)

    integer, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      p_matrix ! N x D+1 matrix of random integers.

    ! Weight of 1st Gaussian, 0 <= a <= 1
    double precision, intent(in) :: a

    !rt1, thl1 = mean of rt, thl for Gaus comp 1
    !rt2, thl2 = mean of rt, thl for Gaus comp 2
    double precision, intent(in) :: rt1, thl1, rt2, thl2

    ! Thermodynamic constants for plumes 1 and 2, units of kg/kg
    double precision, intent(in) :: &
      crt1,  & ! coefficient relating rt, s and t for Gaus comp 1
      cthl1, & ! coeff relating thl, s and t for component 1
      crt2,  & ! coefficient relating rt, s and t for component 2
      cthl2    ! coefficient relating thl, s and t for comp. 2

    ! Latin hypercube variables, i.e. s, t, w, etc.
    double precision, intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd components

    ! Covariance matrices of rt, thl, w for each Gaussian
    ! Columns of Sigma_rtthlw:     1   2   3   4   5   6
    !                              rt  thl w   Nc  Nr  rr
    !
    ! Columns of Sigma_stw, X_nl_one_lev:  1   2   3   4   5   6
    !                                      s   t   w   Nc  Nr  rr
    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma_rtthlw_1, &
      Sigma_rtthlw_2

    ! Cloud fractions for components 1 and 2
    double precision, intent(in) :: &
      cloud_frac1, cloud_frac2 ! cloud fraction associated w/ 1st, 2nd mixture component

    logical, intent(in), dimension(d_variables) :: &
      l_d_variable_lognormal ! Whether a given element of X_nl is lognormal

    ! Output Variables
    ! Total water, theta_l: mean plus perturbations
    double precision, intent(out), dimension(n_micro_calls) :: rt, thl

    double precision, intent(out), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev ! Sample drawn from uniform distribution from particular grid level

    double precision, intent(out), dimension(n_micro_calls,d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal


    ! Local Variables
    integer :: col, sample

    ! Covariance matrices for variables s, t, w for comps 1 & 2
    double precision, dimension(d_variables,d_variables) :: &
      Sigma_stw_1, Sigma_stw_2

    ! Sample of s points that is drawn only from normal distribution
    double precision, dimension(n_micro_calls) :: s_pts

    integer :: is_mellor, it_mellor
!   integer :: i, j
    ! ---- Begin Code ----

    is_mellor = iiLH_rt  ! Mellor's s is at the same index as rt in the Sigma arrays
    it_mellor = iiLH_thl ! Mellor's t is at the same index as thl "  "

    ! Convert each Gaussian from rt-thl-w variables to s-t-w vars.
    call rtpthlp_2_sptp( d_variables, Sigma_rtthlw_1, crt1, cthl1, Sigma_stw_1 )
    call rtpthlp_2_sptp( d_variables, Sigma_rtthlw_2, crt2, cthl2, Sigma_stw_2 )
!   do i = 1, d_variables
!     do j = 1, d_variables
!       write(6,'(e10.3)',advance='no') Sigma_rtthlw_1(i,j)
!     end do
!     write(6,*) ""
!   end do
!   pause
    ! Generate Latin hypercube sample, with one extra dimension
    !    for mixture component.
    call latin_hyper_sample( n_micro_calls, nt_repeat, d_variables+1, p_matrix, X_u_one_lev )

    ! Standard sample for testing purposes when n=2
    ! X_u_one_lev(1,1:(d+1)) = ( / 0.0001d0, 0.46711825945881d0, &
    !             0.58015016959859d0, 0.61894015386778d0, 0.1d0, 0.1d0  / )
    ! X_u_one_lev(2,1:(d+1)) = ( / 0.999d0, 0.63222458307464d0, &
    !             0.43642762850981d0, 0.32291562498749d0, 0.1d0, 0.1d0  / )


    ! Let s PDF (1st column) be a truncated Gaussian.
    ! Take sample solely from cloud points.
    col = is_mellor
    call truncate_gaus_mixt( n_micro_calls, d_variables, col, a, mu1, mu2, & 
                             Sigma_stw_1, Sigma_stw_2, cloud_frac1, cloud_frac2, X_u_one_lev, & 
                             s_pts )

    ! Generate n samples of a d-variate Gaussian mixture
    ! by transforming Latin hypercube points, X_u_one_lev.
    call gaus_mixt_points( n_micro_calls, d_variables, a, mu1, mu2,  & 
                           Sigma_stw_1, Sigma_stw_2, & 
                           cloud_frac1, cloud_frac2, X_u_one_lev, s_pts, X_nl_one_lev )

! Transform s (column 1) and t (column 2) back to rt and thl
! This is only needed if you need rt, thl in your microphysics.
!     call sptp_2_rtpthlp &
!          ( n_micro_calls, d_variables, a, crt1, cthl1, crt2, cthl2, &
!            cloud_frac1, cloud_frac2, X_nl_one_lev(1:n_micro_calls,1), &
!            X_nl_one_lev(1:n_micro_calls,2), &
!            X_u_one_lev, rtp, thlp )
      call st_2_rtthl( n_micro_calls, d_variables, a, rt1, thl1, rt2, thl2, &
                       crt1, cthl1, crt2, cthl2, &
                       cloud_frac1, cloud_frac2, mu1(1), mu2(1), &
                       X_nl_one_lev(1:n_micro_calls,is_mellor), &
                       X_nl_one_lev(1:n_micro_calls,it_mellor), &
                       X_u_one_lev, rt, thl )

! Compute some diagnostics
!       print*, 'C=', a*cloud_frac1 + (1-a)*cloud_frac2
!       print*, 'rtm_anl=', a*rt1+(1-a)*rt2, 'rtm_est=', mean(rt(1:n),n)
!       print*, 'thl_anl=',a*thl1+(1-a)*thl2, 'thlm_est=',mean(thl(1:n),n)
!       print*, 'rtpthlp_coef_est=', corrcoef(rt,thl,n)

    ! Convert lognormal variates (currently Nc and rr) to lognormal
    forall ( sample = 1:n_micro_calls )
      where ( l_d_variable_lognormal )
        X_nl_one_lev(sample,:) = exp( X_nl_one_lev(sample,:) )
      end where
    end forall

! Test diagnostics
!       print*, 'mean X_nl_one_lev(:,2)=', mean(X_nl_one_lev(1:n,2),n)
!       print*, 'mean X_nl_one_lev(:,3)=', mean(X_nl_one_lev(1:n,3),n)
!       print*, 'mean X_nl_one_lev(:,4)=', mean(X_nl_one_lev(1:n,4),n)
!       print*, 'mean X_nl_one_lev(:,5)=', mean(X_nl_one_lev(1:n,5),n)

!       print*, 'std X_nl_one_lev(:,2)=', std(X_nl_one_lev(1:n,2),n)
!       print*, 'std X_nl_one_lev(:,3)=', std(X_nl_one_lev(1:n,3),n)
!       print*, 'std X_nl_one_lev(:,4)=', std(X_nl_one_lev(1:n,4),n)
!       print*, 'std X_nl_one_lev(:,5)=', std(X_nl_one_lev(1:n,5),n)

!       print*, 'corrcoef X_nl_one_lev(:,2:3)=',
!     .            corrcoef( X_nl_one_lev(1:n,2), X_nl_one_lev(1:n,3), n )


    return
  end subroutine sample_points

!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine rtpthlp_2_sptp( d_variables, Sigma_rtthlw, crt, cthl, Sigma_stw )

! Description:
!   Transform covariance matrix from rt', theta_l' coordinates
!   to s', t' coordinates.
!   Use linear approximation for s', t'.
! References:
!   None
!-----------------------------------------------------------------------

    use constants, only: fstderr ! Constant(s)

    implicit none

    ! External

    intrinsic :: matmul

    ! Input Variables

    integer, intent(in) :: d_variables ! Number of variates

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma_rtthlw ! D x D dimensional covariance matrix of rt, thl, w, ...

    double precision, intent(in) :: &
      crt, cthl ! Coefficients that define s', t'

    ! Output Variables

    double precision, intent(out), dimension(d_variables,d_variables) :: &
      Sigma_stw ! Covariance matrix in terms of s', t' ordering of Sigma_stw is s, t, w, ...

    ! Local Variables

    integer j, k
    double precision, dimension(d_variables,d_variables) :: &
      T_transpose, T, M_int

    ! ---- Begin Code ----

    ! Check that matrix is large enough (at least 2x2)
    if (d_variables < 3) then
      write(fstderr,*) 'Error: Input matrix too small in rtpthlp_2_stpthlp.'
      stop
    end if

    ! Transform rt-thl-w matrix to s-t-w
    !    according to Sigma_stw = T*Sigma_rtthlw*T_transpose.

    ! Set up transformation matrix, T_transpose
    T_transpose(1,1) = crt
    T_transpose(2,1) = -1.d0*cthl
    T_transpose(1,2) = crt
    T_transpose(2,2) = 1.d0*cthl

    ! Put zeros in the two off-diagonal blocks
    if (d_variables > 2) then
      do j=1,2
        do k=3,d_variables
          T_transpose(j,k) = 0.d0
          T_transpose(k,j) = 0.d0
        end do
      end do
    end if

    ! Put identity matrix in the lower diagonal block.
    if (d_variables > 2) then
      do j = 3, d_variables
        do k = 3, d_variables
          if (j == k) then
            T_transpose(j,k) = 1.d0
          else
            T_transpose(j,k) = 0.d0
          end if
        end do
      end do
    end if

    ! Multiply to obtain M_int = Sigma_rtthlw * T_transpose
    M_int = matmul( Sigma_rtthlw, T_transpose )

    ! Set up other transformation matrix, T
    T      = T_transpose
    T(1,2) = -1.d0*cthl
    T(2,1) = crt

    ! Perform final matrix multiplication: Sigma_stw = T * M_int
    Sigma_stw = matmul( T, M_int )

    return
  end subroutine rtpthlp_2_sptp
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine latin_hyper_sample( n_micro_calls, nt_repeat, dp1, p_matrix, X)

! Description:
!   Generates a matrix X that contains a Latin Hypercube sample.
!   The sample is uniformly distributed.
! References:
!   See Art B. Owen (2003), ``Quasi-Monte Carlo Sampling,"
!      a chapter from SIGGRAPH 2003
!-------------------------------------------------------------------------------

    use mt95, only: genrand_real3 ! Procedure(s)

    use mt95, only: genrand_real ! Constants

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_micro_calls, & ! `n'   Number of calls to microphysics (normally=2)
      nt_repeat,     & ! `n_t' Num. random samples before sequence repeats (normally=10)
      dp1              !  d+1  Number of variates plus 1 (normally=6)

    integer, intent(in), dimension(n_micro_calls,dp1) :: &
      p_matrix    !n x dp1 array of permuted integers

    ! Output Variables

    double precision, intent(out), dimension(n_micro_calls,dp1) :: &
      X ! n by dp1 matrix, X, each row of which is a dp1-dimensional sample

    ! Local Variables

    real(kind=genrand_real) :: rand ! Random float with a range of (0,1)

    integer :: j, k

    ! ---- Begin Code ----

!  Compute random permutation row by row
!       do j=1,dp1
!       ! Generate a column vector of integers from 0 to n-1,
!       !    whose order is random.
!         call rand_permute( n, p_matrix(1:n,j) )
!       end do

    ! Choose values of sample using permuted vector and random number generator
    do j = 1,n_micro_calls
      do k = 1,dp1
        call genrand_real3( rand ) ! genrand_real3's range is (0,1)
        X(j,k) = (1.0d0/nt_repeat)*(p_matrix(j,k) + rand )
      end do
    end do

!        print*, 'p_matrix(:,1)= ', p_matrix(:,1)
!        print*, 'p_matrix(:,dp1)= ', p_matrix(:,dp1)
!        print*, 'X(:,1)= ', X(:,:)

    return
  end subroutine latin_hyper_sample

!----------------------------------------------------------------------
  subroutine gaus_mixt_points( n_micro_calls, d_variables, a, mu1, mu2, Sigma1, Sigma2, & 
                               cloud_frac1, cloud_frac2, X_u_one_lev, s_pts, X_gm )
! Description:
!   Generates n random samples from a d-dimensional Gaussian-mixture PDF.
!   Uses Latin hypercube method.
! References:
!   None
!----------------------------------------------------------------------

    use constants, only:  &
      fstderr  ! Constant(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, &  ! Number of calls to microphysics (normally=2) 
      d_variables       ! Number of variates (normally=5)

    double precision, intent(in) :: &
      a,     & ! Mixture fraction of Gaussians
      cloud_frac1, cloud_frac2   ! Cloud fraction associated w/ 1st, 2nd mixture component

    double precision, intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd Gaussians

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma1, Sigma2 ! dxd dimensional covariance matrices

    ! Latin hypercube sample from uniform distribution from a particular grid level
    double precision, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev 

    double precision, intent(in), dimension(n_micro_calls) :: &
     s_pts ! n-dimensional vector giving values of s

    ! Output Variables

    double precision, intent(out), dimension(n_micro_calls,d_variables) :: &
      X_gm ! [n by d] matrix, X_gm, each row of which is a d-dimensional sample

    ! Local Variables

    integer :: j, sample
!   double precision, dimension(n_micro_calls) :: std_normal
    double precision, dimension(d_variables) :: std_normal
    double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of a, cloud_frac1,
    ! cloud_frac2.
    if (a > 1.0d0 .or. a < 0.0d0) then
      write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                       'mixture fraction, a, does not lie in [0,1].'
      stop
    end if
    if (cloud_frac1 > 1.0d0 .or. cloud_frac1 < 0.0d0) then
      write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if (cloud_frac2 > 1.0d0 .or. cloud_frac2 < 0.0d0) then
      write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if (a*cloud_frac1 < 0.001d0 .and. (1-a)*cloud_frac2 < 0.001d0) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    do sample = 1, n_micro_calls

      ! From Latin hypercube sample, generate standard normal sample
      do j = 1, d_variables
        std_normal(j) = ltqnorm( X_u_one_lev(sample,j) )
      end do

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
      fraction_1 = ( a*cloud_frac1 ) / ( a*cloud_frac1 + (1-a)*cloud_frac2 )
      if ( X_u_one_lev(sample, d_variables+1) < fraction_1 ) then
        call gaus_condt( d_variables, &
                         std_normal, mu1, Sigma1, s_pts(sample), & 
                         X_gm(sample, 1:d_variables) )
      else
        call gaus_condt( d_variables, &
                         std_normal, mu2, Sigma2, s_pts(sample), & 
                         X_gm(sample, 1:d_variables) )
      end if

      ! Loop to get new sample
    end do

    return
  end subroutine gaus_mixt_points

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine truncate_gaus_mixt( n_micro_calls, d_variables, col, a, mu1, mu2, & 
                  Sigma1, Sigma2, cloud_frac1, cloud_frac2, X_u_one_lev, truncated_column )
! Description:
!   Converts sample points drawn from a uniform distribution
!    to truncated Gaussian points.
! References:
!   None
!-------------------------------------------------------------------------------

    use constants, only:  &
      fstderr  ! Constant(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, &  ! Number of calls to microphysics (normally=2) 
      d_variables,   &  ! Number of variates (normally=5)
      col               ! Scalar indicated which column of X_nl_one_lev to truncate

    double precision, intent(in) :: &
      a,    & ! Mixture fraction of Gaussians
      cloud_frac1, cloud_frac2  ! Cloud fraction associated w/ 1st, 2nd mixture component

    double precision, intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd Gaussians

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma1, Sigma2 ! dxd dimensional covariance matrices

    ! Latin hypercube sample from uniform distribution from a particular grid level
    double precision, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev 

    ! Output Variables

    ! A column vector of length n that is transformed from a Gaussian PDF to truncated Gaussian PDF.
    double precision, intent(out), dimension(n_micro_calls) :: truncated_column

    ! Local Variables

    integer :: sample
    double precision :: s_std
    double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of a, cloud_frac1,
    ! cloud_frac2.
    if ( (a > 1.0d0) .or. (a < 0.0d0) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'mixture fraction, a, does not lie in [0,1].'
      stop
    end if
    if ( (cloud_frac1 > 1.0d0) .or. (cloud_frac1 < 0.0d0) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if ( (cloud_frac2 > 1.0d0) .or. (cloud_frac2 < 0.0d0) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if (a*cloud_frac1 < 0.001d0 .and. (1-a) * cloud_frac2 < 0.001d0) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    ! Make s PDF (1st column) a truncated Gaussian.
    ! This allows us to sample solely from the cloud points.
    do sample = 1, n_micro_calls

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
      fraction_1 = a * cloud_frac1 / ( a * cloud_frac1 + (1.d0 - a) *cloud_frac2 )
      if ( X_u_one_lev( sample, d_variables+1 ) < fraction_1 ) then
        ! Replace first dimension (s) with
        !  sample from cloud (i.e. truncated standard Gaussian)
        s_std = ltqnorm( X_u_one_lev( sample, col ) * cloud_frac1 + (1.d0 - &
          cloud_frac1) )
        ! Convert to nonstandard normal with mean mu1 and variance Sigma1
        truncated_column(sample) =  & 
                   s_std * sqrt( Sigma1(col,col) ) + mu1(col)
      else

        ! Replace first dimension (s) with
        !   sample from cloud (i.e. truncated Gaussian)
        s_std = ltqnorm( X_u_one_lev( sample, col ) * cloud_frac2 + (1.d0 - &
          cloud_frac2) )

        ! Convert to nonstandard normal with mean mu2 and variance Sigma2
        truncated_column(sample) =  & 
                      s_std * sqrt( Sigma2(col,col) ) + mu2(col)
      end if

      ! Loop to get new sample
    end do

    return
  end subroutine truncate_gaus_mixt

!-----------------------------------------------------------------------
  double precision function ltqnorm( p )
! Description:
!   This function is ported to Fortran from the same function written in Matlab, see the following
!   description of this function.  Hongli Jiang, 2/17/2004
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
!   result has a relative error less than 1.15e-9.  A last refinement by
!   Halley's rational method is applied to achieve full machine precision.

!   Author:      Peter J. Acklam
!   Time-stamp:  2003-04-23 08:26:51 +0200
!   E-mail:      pjacklam@online.no
!   URL:         http://home.online.no/~pjacklam
!-----------------------------------------------------------------------

    use constants, only: Pi_DP ! Variable(s)

    implicit none

    ! External

    intrinsic :: log, sqrt

    ! Input Variable(s)

    double precision, intent(in) :: p


    ! Local Variable(s)

    double precision a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, & 
                     c1, c2, c3, c4, c5, c6, d1, d2, d3, d4

    double precision q, r, z, z1, plow, phigh

!       double preciseion e, erf_dp, u

!       Occurs in constants.F now.  Isn't actually used currently.
!        double precision, parameter :: pi=3.1415926d0

! Coefficients in rational approximations.
! equivalent: a(1)=a1, a(2)=a2, and etc, when a(1) is in Matlab.
! Similarly for b, c, and d's
    parameter (a1 = -3.969683028665376d+01,  & 
               a2 = 2.209460984245205d+02, & 
               a3 = -2.759285104469687d+02,  & 
               a4 = 1.383577518672690d+02, & 
               a5 = -3.066479806614716d+01,  & 
               a6 = 2.506628277459239d+00)
    parameter (b1 = -5.447609879822406d+01,  & 
               b2 = 1.615858368580409d+02, & 
               b3 = -1.556989798598866d+02,  & 
               b4 = 6.680131188771972d+01, & 
               b5 = -1.328068155288572d+01)
    parameter (c1 = -7.784894002430293d-03,  & 
               c2 = -3.223964580411365d-01, & 
               c3 = -2.400758277161838d+00,  & 
               c4 = -2.549732539343734d+00, & 
               c5 =  4.374664141464968d+00,  & 
               c6 =  2.938163982698783d+00)
    parameter (d1 =  7.784695709041462d-03,  & 
               d2 =  3.224671290700398d-01, & 
               d3 =  2.445134137142996d+00,  & 
               d4 =  3.754408661907416d+00)

    ! Default initialization
    z = 0.0

!  Define break-points.
    plow  = 0.02425d0
    phigh = 1.d0 - plow

!  Initialize output array. Don't need this in Fortran
!   z = zeros(size(p));

!  Rational approximation for lower region:
    if (p > 0.d0 .and. p < plow) then
      q = sqrt( -2 * log(p) )
      z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ & 
                ((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
!  Rational approximation for central region:
    else if (p >= plow .and. p <= phigh) then
      q = p - 0.5d0
      r = q * q
      z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q & 
                 /(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.d0)
! Rational approximation for upper region:
    else if (p > phigh .and. p < 1.d0) then
      q  = sqrt( -2.d0 * log(1.d0 - p) )
      z  = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) & 
                  /((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
    end if

!  Case when P = 0: z = -inf, to create inf z =-1./0.,
!     to create NaN's inf*inf.
    z1 = 0.d0
    if (p == 0.d0) then
      z = (-1.d0)/z1
    end if

! Case when P = 1:, z=inf
    if(p == 1.d0)then
      z = 1.d0/z1
    end if

!  Cases when output will be NaN:
!   k = p < 0 | p > 1 | isnan(p);
! usually inf*inf --> NaN's.
    if (p < 0.d0 .or. p > 1d0) then
      z = (1.d0/z1)**2
    end if

!  The relative error of the approximation has absolute value less
!  than 1.15e-9.  One iteration of Halley's rational method (third
!  order) gives full machine precision.
! V. Larson 20Feb04: Don't use the following if-end if loop.
!   The value of e is very different than what MATLAB produces,
!   possibly because of
!   poor values of erf from Numerical Recipes.
!   The value is close to MATLAB's
!   if I omit the following if-end if loop.
! End V. Larson comment
!!   k = 0 < p & p < 1;
!       if (p.gt.0 .and. p.lt.1)then
!         e = 0.5*(1.0 - erf_dp(-z/sqrt(2.))) - p          ! error
!         u = e * sqrt(2*pi_dp) * exp(z**2/2)       ! f(z)/df(z)
!         z = z - u/( 1 + z*u/2 )               ! Halley's method
!       end if

! return z as double precision:
    ltqnorm = z

    return
  end function ltqnorm

!----------------------------------------------------------------------
! Description:
!   Using Gaussian conditional distributions given s,
!   convert a standard, uncorrelated Gaussian to one with
!   mean mu and covariance structure Sigma.

! References:
!   Follow M. E. Johnson, ``Multivariate Statistical Simulation," p50.

!----------------------------------------------------------------------
  subroutine gaus_condt( d_variables, std_normal, mu, Sigma, s_pt, & 
                         nonstd_normal )

    use matrix_operations, only: linear_symm_upper_eqn_solve ! Procedure(s)

    implicit none

    ! External
    intrinsic :: sqrt, matmul, reshape

    ! Input Variables

    integer, intent(in) :: d_variables ! Number of variates (normally=5)

    double precision, intent(in), dimension(d_variables) :: &
      std_normal ! nxd matrix of n independent samples from d-variate standard normal distribution

    double precision, intent(in), dimension(d_variables,1) :: &
      mu ! d-dimensional column vector of means of Gaussian

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma ! dxd dimensional covariance matrix

    double precision, intent(in) :: s_pt ! Value of Mellor's s

    ! Output Variables

    ! nxd matrix of n samples from d-variate normal distribution
    !   with mean mu and covariance structure Sigma
    double precision, intent(out) :: &
      nonstd_normal(d_variables)

    ! Local Variables

    integer :: v

    double precision :: mu_two, mu_condt_const

    double precision, dimension(d_variables,1) :: mu_one

    double precision, dimension(d_variables,d_variables) :: &
      Sigma_oneone

    double precision, dimension(d_variables,1) :: &
      Sigma_onetwo

    double precision, dimension(1,d_variables) :: &
      Sigma_twoone

    double precision :: Sigma_twotwo, Sigma_condt

    ! Local intermediate quantities
    double precision, dimension(1,d_variables) :: Sigma_int
    double precision, dimension(1,1) :: dum

    ! First store value of s
    nonstd_normal(1) = s_pt

    ! Loop over variables t, w, . . .
    ! [Conditional distribution of
    !      X_two given X_one] = x_one = std_normal(1:n,v)
    ! is N [ mu_two + Sigma_twoone*Inverse(Sigma_oneone)*(x_one-mu_one),
    !    Sigma_twotwo - Sigma_twoone*Inverse(Sigma_oneone)*Sigma_onetwo]
    ! Here 'one' refers to given variables, 'two' refers to values to find.
    ! Loop over variables t, w, N, rr:
    do v=2, d_variables

      ! Means
      mu_one(1:(v-1),1) = mu(1:(v-1),1)
      mu_two            = mu(v,1)

      ! Matrices used to compute variances
      Sigma_twoone(1,1:(v-1)) = Sigma(v,1:(v-1))
      Sigma_onetwo(1:(v-1),1) = Sigma(1:(v-1),v)
      Sigma_oneone( 1:(v-1) , 1:(v-1) )  & 
           = Sigma( 1:(v-1) , 1:(v-1) )
      Sigma_twotwo = Sigma(v,v)

      ! Check whether input matrix is symmetric.
!         do j=1,(v-1)
!            if ( Sigma_twoone(1,j) .ne. Sigma_onetwo(j,1) ) then
!              print*,'Error: matrix Sigma not symmetric in gaus_condt.'
!               print*,'Sigma in gaus_condt=', Sigma
!            end if
!         end do

      ! Compute an intermediate matrix, Sigma_int(1,1:(v-1)),
      !    that is needed several times below.
      ! Solve A * X = B for X, where X here is sigma_int.
      call linear_symm_upper_eqn_solve &
           ( n = v-1, a = Sigma_oneone(1:v-1,1:v-1), b = Sigma_twoone(1,1:v-1), &
             x = Sigma_int(1,1:v-1) )

      ! Constant scalar needed to compute mean of non-standard normal.
      ! mu_condt_const = mu_two - Sigma_twoone*Sigma_oneone_inv*mu_one

      dum = matmul( Sigma_int(:,1:v-1), mu_one(1:v-1,:) )

      mu_condt_const = mu_two - dum(1,1)

      ! Variance of non-standard normal for variable v.
! Sigma_condt is a scalar.
! Sigma_condt = Sigma_twotwo -
!                    Sigma_twoone*Sigma_oneone_inv*Sigma_onetwo
      dum = matmul( Sigma_int(:,1:v-1), Sigma_onetwo(1:v-1,:) )
      Sigma_condt = Sigma_twotwo - dum(1,1)

      ! Compute next element of nonstd_normal from prior elements
! nonstd_normal(v) = sqrt(Sigma_condt)*std_normal(v)
!      + mu_condt_const ...
!      + Sigma_twoone*Sigma_oneone_inv*nonstd_normal(1:(v-1))
      dum = matmul( Sigma_int(:,1:v-1), &
                    reshape( nonstd_normal(1:v-1), (/v-1, 1/) ) )
      nonstd_normal(v) = sqrt( Sigma_condt )*std_normal(v)  & 
         + mu_condt_const  & 
         + dum(1,1)

      ! Loop to obtain new variable
    end do ! 2..d_variables

    return
  end subroutine gaus_condt

!-----------------------------------------------------------------------
! Description:
!   Converts from s, t variables to rt, thl
!-----------------------------------------------------------------------
  subroutine st_2_rtthl( n_micro_calls, d_variables, a, rt1, thl1, rt2, thl2, & 
                         crt1, cthl1, crt2, cthl2, & 
                         cloud_frac1, cloud_frac2, s1, s2, &
                         s_mellor, &
                         t_mellor, &
                         X_u_one_lev, rt, thl )

    use constants, only:  &
        fstderr  ! Constant(s)

    use error_code, only:  &
        clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input

    integer, intent(in) :: &
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    double precision, intent(in) :: &
      a,           & ! Mixture fraction of Gaussians 'a'        [-]
      rt1, rt2,    & ! n dimensional column vector of rt        [kg/kg]
      thl1, thl2,  & ! n dimensional column vector of thetal    [K]
      crt1, crt2,  & ! Constants from plumes 1 & 2 of rt
      cthl1, cthl2   ! Constants from plumes 1 & 2 of thetal

    double precision, intent(in) :: &
      cloud_frac1, cloud_frac2 ! Cloud fraction associated with 1st / 2nd mixture component

    double precision, intent(in) :: s1, s2

    ! n-dimensional column vector of Mellor's s and t, including mean and perturbation
    double precision, intent(in), dimension(n_micro_calls) :: &
      s_mellor, & 
      t_mellor 

    double precision, dimension(n_micro_calls,d_variables+1), intent(in) :: &
      X_u_one_lev ! n_micro_calls x d_var+1 Latin hypercube sample from uniform distribution 
    !       from a particular grid level

    ! Output variables

    double precision, dimension(n_micro_calls), intent(out) :: &
      rt, thl ! n-dimensional column vectors of rt and thl, including mean and perturbation

    ! Local

    integer :: sample
    double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of a, cloud_frac1,
    ! cloud_frac2.

    if ( a > 1.0d0 .or. a < 0.0d0 ) then
      write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                       'mixture fraction, a, does not lie in [0,1].'
      stop
    end if
    if ( cloud_frac1 > 1.0d0 .or. cloud_frac1 < 0.0d0 ) then
      write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if ( cloud_frac2 > 1.0d0 .or. cloud_frac2 < 0.0d0 ) then
      write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if ( a*cloud_frac1 < 0.001d0 .and. (1-a)*cloud_frac2 < 0.001d0 ) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    do sample = 1, n_micro_calls

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
      fraction_1     = a*cloud_frac1/(a*cloud_frac1+(1-a)*cloud_frac2)
      if ( X_u_one_lev(sample,d_variables+1) < fraction_1 ) then
        rt(sample)  = rt1 + (0.5d0/crt1)*(s_mellor(sample)-s1) +  & 
                           (0.5d0/crt1)*t_mellor(sample)
        thl(sample) = thl1 + (-0.5d0/cthl1)*(s_mellor(sample)-s1) +  & 
                           (0.5d0/cthl1)*t_mellor(sample)
      else
        rt(sample)  = rt2 + (0.5d0/crt2)*(s_mellor(sample)-s2) +  & 
                           (0.5d0/crt2)*t_mellor(sample)
        thl(sample) = thl2 + (-0.5d0/cthl2)*(s_mellor(sample)-s2) +  & 
                           (0.5d0/cthl2)*t_mellor(sample)
      end if

      ! Loop to get new sample
    end do

    return
  end subroutine st_2_rtthl
!-------------------------------------------------------------------------------
  subroutine log_sqd_normalized( Xm, Xp2_on_Xm2, X_tol, &
                                 X1, X2, sX1, sX2, l_small_X )
! Description:
!
! References:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: log

    ! Input Variables
    double precision, intent(in) :: &
      Xm,           & ! Mean X          [units vary]
      Xp2_on_Xm2,   & ! X'^2 / X^2      [-]
      X_tol           ! tolerance on X  [units vary]

    ! Output Variables
    double precision, intent(out) :: &
      X1, X2,  & ! PDF parameters for mean of plume 1, 2   [units vary]
      sX1, sX2   ! PDF parameters for width of plume 1, 2  [units vary]

    logical, intent(out) :: &
      l_small_X

    ! ---- Begin Code ----
    if ( Xm > X_tol ) then
      X1 = 0.5 * log( Xm**2 / ( 1. + Xp2_on_Xm2 ) )
      X2 = X1
      sX1 = log( 1. + Xp2_on_Xm2 )
      sX2 = sX1
      l_small_X = .false.
    else
      X1 = 0.0
      X2 = 0.0
!     sX1 = 0.0
!     sX2 = 0.0
      sX1 = X_tol**2
      sX2 = X_tol**2
      l_small_X = .true.
    end if

    return
  end subroutine log_sqd_normalized

end module generate_lh_sample_mod

