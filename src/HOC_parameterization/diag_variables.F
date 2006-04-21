!-----------------------------------------------------------------------
! $Id: diag_variables.F,v 1.1 2006-04-21 21:50:56 dschanen Exp $
        module diagnostic_variables

! This module contains definitions of all diagnostic
! arrays used in the single column model, as well as subroutines
! to allocate, deallocate and initialize them.

! Note that while these are all same dimension, there is a
! thermodynamic and momentum grid and they have different levels
!-----------------------------------------------------------------------

        implicit none

! Diagnostic variables

!       real, target, allocatable :: Scm(:)   ! PDF width parameter: m point
        real, target, allocatable :: Sct(:)   ! PDF width parameter: t point

        real, target, allocatable :: ug(:)    ! u geostrophic wind
        real, target, allocatable :: vg(:)    ! v geostrophic wind

        real, target, allocatable :: wmt(:)    ! imposed large scale w: t point
        real, target, allocatable :: wmm(:)    ! imposed large scale w: m point

        real, allocatable, target :: cf(:)        ! cloud fraction
        real, allocatable, target :: thvm(:)      ! virtual potential temp 
        real, allocatable, target :: rcm(:)       ! cloud water mixing ratio
        real, allocatable, target :: Ncm(:)       ! cloud droplet number conc. ! Brian
        real, allocatable, target :: sed_rcm(:)   ! sedimentation factor in rcm ! Brian
        real, allocatable, target :: rsm(:)       ! saturation mixing ratio  ! Brian
        real, allocatable, target :: rrm(:)       ! rain water mixing ratio  ! Brian
        real, allocatable, target :: cond_rrm(:)  ! condensation factor in rrm ! Brian
        real, allocatable, target :: auto_rrm(:)  ! autoconversion factor in rrm ! Brian
        real, allocatable, target :: accr_rrm(:)  ! accretion factor in rrm ! Brian
        real, allocatable, target :: sed_rrm(:)   ! sedimentation output for rrm ! Brian
        real, allocatable, target :: diff_rrm(:)  ! diffusion output for rrm ! Brian
        real, allocatable, target :: Nrm(:)       ! rain droplet number conc.  ! Brian
        real, allocatable, target :: Vrr(:)       ! sedimentation velocity of rain ! Brian
        real, allocatable, target :: mean_vol_rad(:) ! rain drop mean volume radius ! Brian
        real, allocatable, target :: rain_rate(:) ! rainfall rate in mm/day ! Brian
        real, allocatable, target :: Fprec(:)     ! precipitation flux in W/m^2 (momentum point) ! Brian
        real, allocatable, target :: Fcsed(:)     ! cloud water sedimentation flux (momentum point) ! Brian
                                                  ! Note: units are [kg(cloud water)]/[s*m^2]

        real, allocatable, target :: Frad(:)      ! Radiative flux (momentum point)
        real, allocatable, target :: Frad_LW(:)   ! Long_wave radiative flux (momentum point)
        real, allocatable, target :: Frad_SW(:)   ! Short-wave radiative flux (momentum point)

        real, allocatable, target :: radht(:)     ! SW + LW heating rate
        real, allocatable, target :: radht_SW(:)  ! Short-wave heating rate
        real, allocatable, target :: radht_LW(:)  ! Long-wave heating rate

        real, allocatable, target :: shear(:)     ! wind shear production

        ! Second order moments

        real, target, allocatable :: wprcp(:)     ! w'rc'
        real, target, allocatable :: thlprcp(:)   ! thl'rc'
        real, target, allocatable :: rtprcp(:)    ! rt'rc'

        ! Third order moments

        real, target, allocatable :: wpthlp2(:)   ! w'thl'^2
        real, target, allocatable :: wp2thlp(:)   ! w'^2thl'
        real, target, allocatable :: wprtp2(:)    ! w'rt'^2
        real, target, allocatable :: wp2rtp(:)    ! w'^2rt'
        real, target, allocatable :: wprtpthlp(:) ! w'rt'thl'
        real, target, allocatable :: wp2rcp(:)    ! w'^2rc'

        ! Fourth order moments

        real, allocatable, target :: wp4(:)

        ! Buoyancy related moments

        real, allocatable, target :: wpthvp(:)
        real, allocatable, target :: rtpthvp(:)
        real, allocatable, target :: thlpthvp(:)
        real, allocatable, target :: wp2thvp(:)

        real, allocatable, target :: Kht(:)        ! eddy diffusivity: zt grid
        real, allocatable, target :: Khm(:)        ! eddy diffusivity: zm grid

        ! mixing lengths
        real, allocatable, target :: Lscale(:), lup(:), ldown(:)

        real, allocatable, target :: em(:)
        real, allocatable, target :: taut(:)
!       real, allocatable, target :: taum(:)

        ! Surface data
        real :: ustar            ! average value of friction velocity (m/s)

#ifdef SCALARS
! New mixing scheme variables

        real, target, allocatable :: sclrpthvp(:,:)      ! sclr'th_v'
        real, target, allocatable :: sclrprtp(:,:)       ! sclr'rt'
        real, target, allocatable :: sclrp2(:,:)         ! sclr'^2
        real, target, allocatable :: sclrpthlp(:,:)      ! sclr'th_l'
        real, target, allocatable :: sclrprcp(:,:)       ! sclr'rc'

        real, target, allocatable :: wpsclrp(:,:)        ! w'sclr'
        real, target, allocatable :: wp2sclrp(:,:)       ! w'^2 sclr'
        real, target, allocatable :: wpsclrp2(:,:)       ! w'sclr'^2
        real, target, allocatable :: wpsclrprtp(:,:)     ! w'sclr'rt'
        real, target, allocatable :: wpsclrpthlp(:,:)    ! w'sclr'thl'

        real, target, allocatable :: wpedsclrp(:,:)
#endif /*SCALARS*/

! Variables for tuning
#ifdef STATS
        real, target, allocatable :: wp2zt(:)      ! w'^2 in zt
        real, target, allocatable :: thlp2zt(:)    
        real, target, allocatable :: wpthlpzt(:)   
        real, target, allocatable :: wprtpzt(:)
        real, target, allocatable :: rtp2zt(:)
        real, target, allocatable :: rtpthlpzt(:)
#endif /*STATS*/

! Variables needed for the pdf closure scheme
!
!       pdf_parms contains the parameters of the pdf:
!
!        pdf_parms(1) = w1
!        pdf_parms(2) = w2
!        pdf_parms(3) = sw1
!        pdf_parms(4) = sw2
!        pdf_parms(5) = rt1
!        pdf_parms(6) = rt2
!        pdf_parms(7) = srt1
!        pdf_parms(8) = srt2
!        pdf_parms(9) = thl1
!        pdf_parms(10) = thl2
!        pdf_parms(11) = sthl1
!        pdf_parms(12) = sthl2
!        pdf_parms(13) = a
!        pdf_parms(14) = rc1
!        pdf_parms(15) = rc2
!        pdf_parms(16) = rsl1
!        pdf_parms(17) = rsl2
!        pdf_parms(18) = R1
!        pdf_parms(19) = R2
!        pdf_parms(20) = s1
!        pdf_parms(21) = s2
!        pdf_parms(22) = ss1
!        pdf_parms(23) = ss2
!        pdf_parms(24) = rrtthl

        real, allocatable, target :: pdf_parms(:,:)

! Latin Hypercube arrays.  Vince Larson 22 May 2005
        real, allocatable, target :: AKm_est(:)   ! Kessler ac estimate
        real, allocatable, target :: AKm(:)       ! Exact Kessler ac
        real, allocatable, target :: AKstd(:)     ! St dev of exact Kessler ac
        real, allocatable, target :: AKstd_cld(:) ! Stdev of exact w/in cloud ac
        real, allocatable, target :: rcm_est(:)   ! Monte Carlo rcm estimate
        real, allocatable, target :: AKm_rcm(:)   ! Kessler ac based on rcm
        real, allocatable, target :: AKm_rcc(:)   ! Kessler ac based on rcm/cf
! End Latin hypercube arrays        

        contains

!-----------------------------------------------------------------------
!  Allocates and initializes prognostic scalar and array variables 
!  for the HOC model code
!-----------------------------------------------------------------------
        subroutine setup_diagnostic_variables( nzmax )
        use model_flags, only: kk_rain, LH_on
        use constants, only: sclrm_dimension, emin

        implicit none
 
        integer, parameter :: pdf_dimension = 26

        integer, intent(in) :: nzmax

!   --- Allocation ---

! Diagnostic variables

!       allocate( Scm(1:nzmax) )       ! PDF width parameter: m point
        allocate( Sct(1:nzmax) )       ! PDF width parameter: t point
        allocate( ug(1:nzmax) )        ! u geostrophic wind
        allocate( vg(1:nzmax) )        ! v geostrophic wind

        allocate( wmt(1:nzmax) )       ! imposed large scale w: t point
        allocate( wmm(1:nzmax) )       ! imposed large scale w: m point
 
        allocate( cf(1:nzmax) )        ! cloud fraction
        allocate( thvm(1:nzmax) )      ! virtual potential temperature
        allocate( rcm(1:nzmax) )       ! cloud water mixing ratio

        allocate( Ncm(1:nzmax) )       ! cloud droplet number conc. ! Brian
        allocate( rsm(1:nzmax) )       ! saturation mixing ratio  ! Brian
        if ( kk_rain ) then
          allocate( sed_rcm(1:nzmax) )   ! sedimentation factor in rcm ! Brian
          allocate( rrm(1:nzmax) )       ! rain water mixing ratio  ! Brian
          allocate( cond_rrm(1:nzmax) )  ! condensation factor in rrm ! Brian
          allocate( auto_rrm(1:nzmax) )  ! autoconversion factor in rrm ! Brian
          allocate( accr_rrm(1:nzmax) )  ! accretion factor in rrm ! Brian
          allocate( sed_rrm(1:nzmax) )   ! sedimentation output for rrm ! Brian
          allocate( diff_rrm(1:nzmax) )  ! diffusion output for rrm ! Brian
          allocate( Nrm(1:nzmax) )       ! rain droplet number conc.  ! Brian
          allocate( Vrr(1:nzmax) )       ! sedimentation velocity of rain ! Brian
          allocate( mean_vol_rad(1:nzmax) ) ! rain drop mean volume radius ! Brian
          allocate( rain_rate(1:nzmax) ) ! rainfall rate in mm/day ! Brian
          allocate( Fprec(1:nzmax) )     ! precipitation flux in W/m^2 (momentum point) ! Brian
          allocate( Fcsed(1:nzmax) )     ! cloud water sedimentation flux (momentum point) ! Brian
        end if ! kk_rain

        allocate( Frad(1:nzmax) )      ! radiative flux (momentum point)
        allocate( Frad_LW(1:nzmax) )   ! Long_wave radiative flux (momentum point)
        allocate( Frad_SW(1:nzmax) )   ! Short-wave radiative flux (momentum point)

        allocate( radht(1:nzmax) )     ! SW + LW heating rate
        allocate( radht_SW(1:nzmax) )  ! Short-wave heating rate
        allocate( radht_LW(1:nzmax) )  ! Long-wave heating rate

        allocate( shear(1:nzmax) )     ! wind shear production

        ! Second order moments

        allocate( wprcp(1:nzmax) )     ! w'rc'
        allocate( thlprcp(1:nzmax) )   ! thl'rc'
        allocate( rtprcp(1:nzmax) )    ! rt'rc'

        ! Third order moments

        allocate( wpthlp2(1:nzmax) )   ! w'thl'^2
        allocate( wp2thlp(1:nzmax) )   ! w'^2thl'
        allocate( wprtp2(1:nzmax) )    ! w'rt'^2
        allocate( wp2rtp(1:nzmax) )    ! w'^2rt'
        allocate( wprtpthlp(1:nzmax) ) ! w'rt'thl'
        allocate( wp2rcp(1:nzmax) )    ! w'^2rc'

        ! Fourth order moments

        allocate( wp4(1:nzmax) )

        ! Buoyancy related moments

        allocate( wpthvp(1:nzmax) )
        allocate( rtpthvp(1:nzmax) )
        allocate( thlpthvp(1:nzmax) )
        allocate( wp2thvp(1:nzmax) )

        allocate( Kht(1:nzmax) )
        allocate( Khm(1:nzmax) )

        allocate( em(1:nzmax) )
        allocate( Lscale(1:nzmax) )
        allocate( lup(1:nzmax) )
        allocate( ldown(1:nzmax) )
        allocate( taut(1:nzmax) )
!       allocate( taum(1:nzmax) )

#ifdef STATS
! Tuning Variables
        allocate( wp2zt(1:nzmax) )     ! w'^2 on thermo. grid
        allocate( thlp2zt(1:nzmax) )   ! thl'^2 on thermo. grid
        allocate( wpthlpzt(1:nzmax) )  ! w'thl' on thermo. grid
        allocate( wprtpzt(1:nzmax) )   ! w'rt' on thermo. grid
        allocate( rtp2zt(1:nzmax) )    ! rt'^2 on thermo. grid
        allocate( rtpthlpzt(1:nzmax) ) ! rt'thl' on thermo. grid
#endif /*STATS*/

! Array fpr pdf closure scheme

        allocate( pdf_parms(1:nzmax,1:pdf_dimension) )

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
        if ( LH_on ) then
          allocate( AKm_est(1:nzmax) )    ! Kessler ac estimate
          allocate( AKm(1:nzmax) )        ! Exact Kessler ac
          allocate( AKstd(1:nzmax) )      ! St dev of exact Kessler ac
          allocate( AKstd_cld(1:nzmax) )  ! St dev of exact w/in cloud Kessler ac
          allocate( rcm_est(1:nzmax) )      ! Monte Carlo rcm estimate
          allocate( AKm_rcm(1:nzmax) )      ! Kessler ac based on rcm
          allocate( AKm_rcc(1:nzmax) )      ! Kessler ac based on rcm/cf
        end if ! LH_on
! End of variables for Latin hypercube.

#ifdef SCALARS

! Variables for new mixing scheme

        allocate( sclrprtp(1:nzmax, 1:sclrm_dimension) )
        allocate( sclrp2(1:nzmax, 1:sclrm_dimension) )
        allocate( sclrpthvp(1:nzmax, 1:sclrm_dimension) )
        allocate( sclrpthlp(1:nzmax, 1:sclrm_dimension) )
        allocate( sclrprcp(1:nzmax, 1:sclrm_dimension) )

        allocate( wpsclrp(1:nzmax, 1:sclrm_dimension) )
        allocate( wp2sclrp(1:nzmax, 1:sclrm_dimension) )
        allocate( wpsclrp2(1:nzmax, 1:sclrm_dimension) )
        allocate( wpsclrprtp(1:nzmax, 1:sclrm_dimension) )
        allocate( wpsclrpthlp(1:nzmax, 1:sclrm_dimension) )

! Eddy Diff. Scalar

        allocate( wpedsclrp(1:nzmax, 1:sclrm_dimension) )
#endif /*SCALARS*/

!   --- Initializaton ---

! Diagnostic variables

!       Scm = 0.0      ! PDF width parameter: m point
        Sct = 0.0      ! PDF width parameter: t point
        ug  = 0.0      ! u geostrophic wind
        vg  = 0.0      ! v geostrophic wind

        wmt = 0.0      ! imposed large scale w: t point
        wmm = 0.0      ! imposed large scale w: m point
 
        cf           = 0.0  ! cloud fraction
        thvm         = 0.0  ! virtual potential temperature
        rcm          = 0.0  ! cloud water mixing ratio
        Ncm          = 0.0  ! cloud droplet number conc. ! Brian
        rsm          = 0.0  ! saturation mixing ratio  ! Brian
        if ( kk_rain ) then
          sed_rcm      = 0.0  ! sedimentation factor in rcm ! Brian
          rrm          = 0.0  ! rain water mixing ratio  ! Brian
          cond_rrm     = 0.0  ! condensation factor in rrm ! Brian
          auto_rrm     = 0.0  ! autoconversion factor in rrm ! Brian
          accr_rrm     = 0.0  ! accretion factor in rrm ! Brian
          sed_rrm      = 0.0  ! sedimentation output for rrm ! Brian
          diff_rrm     = 0.0  ! diffusion output for rrm ! Brian
          Nrm          = 0.0  ! rain droplet number conc.  ! Brian
          Vrr          = 0.0  ! sedimentation velocity of rain ! Brian
          mean_vol_rad = 0.0  ! rain drop mean volume radius ! Brian
          rain_rate    = 0.0  ! rainfall rate in mm/day ! Brian
          Fprec        = 0.0  ! precipitation flux in W/m^2 (momentum point) ! Brian
          Fcsed        = 0.0  ! cloud water sedimentation flux (momentum point) ! Brian
                              ! Note: units are [kg(cloud water)]/[s*m^2]
        end if ! kk_rain

        Frad    = 0.0
        Frad_LW = 0.0  ! Long_wave radiative flux (momentum point)
        Frad_SW = 0.0  ! Short-wave radiative flux (momentum point)


        shear = 0.0    ! wind shear production

        ! Second order moments

        wprcp   = 0.0
        thlprcp = 0.0
        rtprcp  = 0.0

        ! Third order moments

        wpthlp2   = 0.0
        wp2thlp   = 0.0
        wprtp2    = 0.0
        wp2rtp    = 0.0
        wp2rcp    = 0.0
        wprtpthlp = 0.0

        ! Fourth order moments

        wp4 = 0.0

        ! Buoyancy related moments

        wpthvp   = 0.0
        rtpthvp  = 0.0
        thlpthvp = 0.0
        wp2thvp  = 0.0

        Kht      = 0.0
        Khm      = 0.0

        em       = emin
        Lscale   = 0.0
        lup      = 0.0
        ldown    = 0.0
        taut     = 0.0
!       taum     = 0.0

! Array fpr pdf closure scheme

        pdf_parms(:,:) = 0.0

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
        if ( LH_on ) then
          AKm_est       = 0.0  ! Kessler ac estimate
          AKm           = 0.0  ! Exact Kessler ac
          AKstd         = 0.0  ! St dev of exact Kessler ac
          AKstd_cld     = 0.0  ! St dev of exact w/in cloud Kessler ac
          rcm_est       = 0.0  ! Monte Carlo rcm estimate
          AKm_rcm       = 0.0  ! Kessler ac based on rcm
          AKm_rcc       = 0.0  ! Kessler ac based on rcm/cf
        end if ! LH_on
! End of variables for Latin hypercube.

#ifdef SCALARS
! New mixing scheme Variables

        sclrprtp(:,:)      = 0.0
        sclrp2(:,:)        = 0.0
        sclrpthvp(:,:)     = 0.0
        sclrpthlp(:,:)     = 0.0
        sclrprcp(:,:)      = 0.0

        wpsclrp(:,:)       = 0.0
        wp2sclrp(:,:)      = 0.0
        wpsclrp2(:,:)      = 0.0
        wpsclrprtp(:,:)    = 0.0
        wpsclrpthlp(:,:)   = 0.0

        wpedsclrp(:,:)     = 0.0
#endif /*SCALARS*/

        return
        end subroutine setup_diagnostic_variables

!------------------------------------------------------------------------
! Subroutine to deallocate variables defined in module global
!
        subroutine cleanup_diagnostic_variables( )
        use model_flags

        implicit none

! Diagnostic variables

!       deallocate( Scm )       ! PDF width parameter: m point
        deallocate( Sct )       ! PDF width parameter: t point
        deallocate( ug )        ! u geostrophic wind
        deallocate( vg )        ! v geostrophic wind

        deallocate( wmt )       ! imposed large scale w: t point
        deallocate( wmm )       ! imposed large scale w: m point
 
        deallocate( cf )        ! cloud fraction
        deallocate( thvm )      ! virtual potential temperature
        deallocate( rcm )       ! cloud water mixing ratio
        deallocate( Ncm )       ! cloud droplet number conc. ! Brian
        deallocate( rsm )       ! saturation mixing ratio  ! Brian
        if ( kk_rain ) then
          deallocate( sed_rcm )   ! sedimentation factor in rcm ! Brian
          deallocate( rrm )       ! rain water mixing ratio  ! Brian
          deallocate( cond_rrm )  ! condensation factor in rrm ! Brian
          deallocate( auto_rrm )  ! autoconversion factor in rrm ! Brian
          deallocate( accr_rrm )  ! accretion factor in rrm ! Brian
          deallocate( sed_rrm )   ! sedimentation output for rrm ! Brian
          deallocate( diff_rrm )  ! diffusion output for rrm ! Brian
          deallocate( Nrm )       ! rain droplet number conc.  ! Brian
          deallocate( Vrr )       ! sedimentation velocity of rain ! Brian
          deallocate( mean_vol_rad ) ! rain drop mean volume radius ! Brian
          deallocate( rain_rate ) ! rainfall rate in mm/day ! Brian
          deallocate( Fprec )     ! precipitation flux in W/m^2 (momentum point) ! Brian
          deallocate( Fcsed )     ! cloud water sedimentation flux (momentum point) ! Brian
                                  ! Note: units are [kg(cloud water)]/[s*m^2]
        end if ! kk_rain

        deallocate( Frad )      ! radiative flux (momentum point)
        deallocate( Frad_LW )   ! Long_wave radiative flux (momentum point)
        deallocate( Frad_SW )   ! Short-wave radiative flux (momentum point)

        deallocate( radht )     ! SW + LW heating rate
        deallocate( radht_SW )  ! Short-wave heating rate
        deallocate( radht_LW )  ! Long-wave heating rate

        deallocate( shear )     ! wind shear production

        ! Second order moments

        deallocate( wprcp )     ! w'rc'
        deallocate( thlprcp )   ! thl'rc'
        deallocate( rtprcp )    ! rt'rc'

        ! Third order moments

        deallocate( wpthlp2 )   ! w'thl'^2
        deallocate( wp2thlp )   ! w'^2thl'
        deallocate( wprtp2 )    ! w'rt'^2
        deallocate( wp2rtp )    ! w'^2rt'
        deallocate( wprtpthlp ) ! w'rt'thl'
        deallocate( wp2rcp )    ! w'^2rc'

        ! Fourth order moments

        deallocate( wp4 )

        ! Buoyancy related moments

        deallocate( wpthvp )
        deallocate( rtpthvp )
        deallocate( thlpthvp )
        deallocate( wp2thvp )

        deallocate( Kht )
        deallocate( Khm )

        deallocate( em )
        deallocate( Lscale )
        deallocate( lup )
        deallocate( ldown )
        deallocate( taut )
!       deallocate( taum )

! Tuning Variables
#ifdef STATS
        deallocate( wp2zt )     ! w'^2 in t
        deallocate( thlp2zt )   ! th_l'^2 in t 
        deallocate( wpthlpzt )  ! w'th_l' in t
        deallocate( wprtpzt )   ! w'rt' in t
        deallocate( rtp2zt )    ! rt'^2 in t
        deallocate( rtpthlpzt ) ! rt'th_l' in t
#endif /*STATS*/


! Array for pdf closure scheme

        deallocate( pdf_parms )

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
        if ( LH_on ) then
          deallocate( AKm_est )   ! Kessler ac estimate
          deallocate( AKm )       ! Exact Kessler ac
          deallocate( AKstd )     ! St dev of exact Kessler ac
          deallocate( AKstd_cld ) ! St dev of exact w/in cloud Kessler ac
          deallocate( rcm_est )   ! Monte Carlo rcm estimate
          deallocate( AKm_rcm )   ! Kessler ac based on rcm
          deallocate( AKm_rcc )   ! Kessler ac based on rcm/cf
        end if ! LH_on
! End of variables for Latin hypercube.

#ifdef SCALARS
! New mixing scheme structures

        deallocate( sclrprtp )
        deallocate( sclrp2 )
        deallocate( sclrpthvp )
        deallocate( sclrpthlp )
        deallocate( sclrprcp )

        deallocate( wpsclrp )
        deallocate( wp2sclrp )
        deallocate( wpsclrp2 )
        deallocate( wpsclrprtp )
        deallocate( wpsclrpthlp )

        deallocate( wpedsclrp )
#endif /*SCALARS*/

        return
        end subroutine cleanup_diagnostic_variables

        end module diagnostic_variables
