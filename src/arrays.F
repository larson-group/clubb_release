!------------------------------------------------------------------------
! $Id: arrays.F,v 1.24 2005-10-27 19:45:19 dschanen Exp $

! This module contains definitions of all prognostic and diagnostic
! arrays used in the single column model, as well as subroutines
! to allocate and deallocate them.
!------------------------------------------------------------------------
        module model_arrays

        implicit none

! Prognostic variables

        real, target, allocatable :: um(:)      ! u wind
        real, target, allocatable :: vm(:)      ! v wind
        real, target, allocatable :: thlm(:)    ! liquid potential temperature
        real, target, allocatable :: rtm(:)     ! total water mixing ratio
        real, target, allocatable :: wp2(:)     ! w'^2
        real, target, allocatable :: rtp2(:)    ! rt'^2
        real, target, allocatable :: thlp2(:)   ! thl'^2
        real, target, allocatable :: wprtp(:)   ! w'rt'
        real, target, allocatable :: wpthlp(:)  ! w'thl'
        real, target, allocatable :: rtpthlp(:) ! rt'thl'
        real, target, allocatable :: wp3(:)     ! w'^3

#ifdef SCALARS
! New mixing scheme variables

        real, target, allocatable :: sclrm(:,:)          ! passive scalars
        real, target, allocatable :: sclrm_forcing(:,:)  ! scalars' forcing

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
#endif SCALARS

! Tendency arrays for prognostic variables

        real, allocatable :: umt(:)             ! u wind
        real, allocatable :: vmt(:)             ! v wind

! Diagnostic variables

!        real, target, allocatable :: Sc(:)        ! PDF width parameter
        real, target, allocatable :: Scm(:)       ! PDF width parameter: m point
        real, target, allocatable :: Sct(:)       ! PDF width parameter: t point
        real, target, allocatable :: ug(:)        ! u geostrophic wind
        real, target, allocatable :: vg(:)        ! v geostrophic wind
        real, target, allocatable :: thlm_forcing(:)   ! thlm ls forcing
        real, target, allocatable :: rtm_forcing(:)    ! rtm ls forcing

        real, target, allocatable :: wmt(:)        ! imposed large scale w: t point
        real, target, allocatable :: wmm(:)        ! imposed large scale w: m point

        real, allocatable, target :: p(:)
                         ! pressure (Pa) on thermodynamic points
        real, allocatable, target :: exner(:)        
                         ! exner = ( p / p0 ) ** kappa

        real, allocatable, target :: rhot(:)      ! density (kg/m^3)
        real, allocatable, target :: rhom(:)      ! density (kg/m^3)

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
        real, allocatable, target :: subs_rrm(:)  ! subsidence output for rrm ! Brian
        real, allocatable, target :: Nrm(:)       ! rain droplet number conc.  ! Brian
        real, allocatable, target :: cond_Nrm(:)  ! condensation factor in Nrm ! Brian
        real, allocatable, target :: auto_Nrm(:)  ! autoconversion factor in Nrm ! Brian
        real, allocatable, target :: sed_Nrm(:)   ! sedimentation output for Nrm ! Brian
        real, allocatable, target :: diff_Nrm(:)  ! diffusion output for Nrm ! Brian
        real, allocatable, target :: subs_Nrm(:)  ! subsidence output for Nrm ! Brian
        real, allocatable, target :: Vr(:)        ! sedimentation velocity of rain ! Brian
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

        real, allocatable, target :: upwp(:)      ! vertical u momentum flux
        real, allocatable, target :: vpwp(:)      ! vertical v momentum flux

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

        real, allocatable, target :: Lscale(:), lup(:), ldown(:)
                                                   ! mixing lengths
        real, allocatable, target :: em(:)
        real, allocatable, target :: taut(:), taum(:)

! Temporary array

!        real, allocatable :: tmp1(:), tmp2(:)

! Variables for tuning
#ifdef STATS
        real, target, allocatable :: wp2zt(:)      ! w'^2 in zt
        real, target, allocatable :: thlp2zt(:)    
        real, target, allocatable :: wpthlpzt(:)   
        real, target, allocatable :: wprtpzt(:)
        real, target, allocatable :: rtp2zt(:)
        real, target, allocatable :: rtpthlpzt(:)
#endif STATS

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


        end module model_arrays

!------------------------------------------------------------------------
! Subroutine to allocate variables defined in module global
!------------------------------------------------------------------------
        subroutine allocate_model_arrays( nzmax )

        use constants, only: sclrm_dimension
        use model_arrays

        implicit none

        integer, intent(in) :: nzmax

! Prognostic variables

        allocate( um(1:nzmax) )        ! u wind
        allocate( vm(1:nzmax) )        ! v wind
        allocate( thlm(1:nzmax) )      ! liquid potential temperature
        allocate( rtm(1:nzmax) )       ! total water mixing ratio
        allocate( wp2(1:nzmax) )       ! w'^2
        allocate( rtp2(1:nzmax) )      ! rt'^2
        allocate( thlp2(1:nzmax) )     ! thl'^2
        allocate( wprtp(1:nzmax) )     ! w'rt'
        allocate( wpthlp(1:nzmax) )    ! w'thl'
        allocate( rtpthlp(1:nzmax) )   ! rt'thlp'
        allocate( wp3(1:nzmax) )       ! w'^3

c Tendency arrays for prognostic variables

        allocate( umt(1:nzmax) )       ! u wind
        allocate( vmt(1:nzmax) )       ! v wind

c Diagnostic variables

!        allocate( Sc(1:nzmax) )        ! PDF width parameter
        allocate( Scm(1:nzmax) )       ! PDF width parameter: m point
        allocate( Sct(1:nzmax) )       ! PDF width parameter: t point
        allocate( ug(1:nzmax) )        ! u geostrophic wind
        allocate( vg(1:nzmax) )        ! v geostrophic wind
        allocate( thlm_forcing(1:nzmax) ) ! thlm ls forcing
        allocate( rtm_forcing(1:nzmax) )  ! rtm ls forcing

        allocate( wmt(1:nzmax) )       ! imposed large scale w: t point
        allocate( wmm(1:nzmax) )       ! imposed large scale w: m point
 
        allocate( p(1:nzmax) )         ! pressure (pascals)
        allocate( exner(1:nzmax) )     ! exner
        allocate( rhot(1:nzmax) )      ! density: t points
        allocate( rhom(1:nzmax) )      ! density: m points

        allocate( cf(1:nzmax) )        ! cloud fraction
        allocate( thvm(1:nzmax) )      ! virtual potential temperature
        allocate( rcm(1:nzmax) )       ! cloud water mixing ratio
        allocate( Ncm(1:nzmax) )       ! cloud droplet number conc. ! Brian
        allocate( sed_rcm(1:nzmax) )   ! sedimentation factor in rcm ! Brian
        allocate( rsm(1:nzmax) )       ! saturation mixing ratio  ! Brian
        allocate( rrm(1:nzmax) )       ! rain water mixing ratio  ! Brian
        allocate( cond_rrm(1:nzmax) )  ! condensation factor in rrm ! Brian
        allocate( auto_rrm(1:nzmax) )  ! autoconversion factor in rrm ! Brian
        allocate( accr_rrm(1:nzmax) )  ! accretion factor in rrm ! Brian
        allocate( sed_rrm(1:nzmax) )   ! sedimentation output for rrm ! Brian
        allocate( diff_rrm(1:nzmax) )  ! diffusion output for rrm ! Brian
        allocate( subs_rrm(1:nzmax) )  ! subsidence output for rrm ! Brian
        allocate( Nrm(1:nzmax) )       ! rain droplet number conc.  ! Brian
        allocate( cond_Nrm(1:nzmax) )  ! condensation factor in Nrm ! Brian
        allocate( auto_Nrm(1:nzmax) )  ! autoconversion factor in Nrm ! Brian
        allocate( sed_Nrm(1:nzmax) )   ! sedimentation output for Nrm ! Brian
        allocate( diff_Nrm(1:nzmax) )  ! diffusion output for Nrm ! Brian
        allocate( subs_Nrm(1:nzmax) )  ! subsidence output for Nrm ! Brian
        allocate( Vr(1:nzmax) )        ! sedimentation velocity of rain ! Brian
        allocate( mean_vol_rad(1:nzmax) ) ! rain drop mean volume radius ! Brian
        allocate( rain_rate(1:nzmax) ) ! rainfall rate in mm/day ! Brian
        allocate( Fprec(1:nzmax) )     ! precipitation flux in W/m^2 (momentum point) ! Brian
        allocate( Fcsed(1:nzmax) )     ! cloud water sedimentation flux (momentum point) ! Brian
                                       ! Note: units are [kg(cloud water)]/[s*m^2]

        allocate( Frad(1:nzmax) )      ! radiative flux (momentum point)
        allocate( Frad_LW(1:nzmax) )   ! Long_wave radiative flux (momentum point)
        allocate( Frad_SW(1:nzmax) )   ! Short-wave radiative flux (momentum point)

        allocate( radht(1:nzmax) )     ! SW + LW heating rate
        allocate( radht_SW(1:nzmax) )  ! Short-wave heating rate
        allocate( radht_LW(1:nzmax) )  ! Long-wave heating rate

        allocate( upwp(1:nzmax) )      ! vertical u momentum flux
        allocate( vpwp(1:nzmax) )      ! vertical v momentum flux

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
        allocate( taum(1:nzmax) )

! Temporary array

!        allocate( tmp1(1:nzmax) )
!        allocate( tmp2(1:nzmax) )

! Tuning Variables
#ifdef STATS
        allocate( wp2zt(1:nzmax) )     ! w'^2 in zt
        allocate( thlp2zt(1:nzmax) )
        allocate( wpthlpzt(1:nzmax) )
        allocate( wprtpzt(1:nzmax) )
        allocate( rtp2zt(1:nzmax) )
        allocate( rtpthlpzt(1:nzmax) )
#endif STATS

! Array fpr pdf closure scheme

        allocate( pdf_parms(1:nzmax,1:26) )
! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
        allocate( AKm_est(1:nzmax) )    ! Kessler ac estimate
        allocate( AKm(1:nzmax) )        ! Exact Kessler ac
        allocate( AKstd(1:nzmax) )      ! St dev of exact Kessler ac
        allocate( AKstd_cld(1:nzmax) )  ! St dev of exact w/in cloud Kessler ac
        allocate( rcm_est(1:nzmax) )      ! Monte Carlo rcm estimate
        allocate( AKm_rcm(1:nzmax) )      ! Kessler ac based on rcm
        allocate( AKm_rcc(1:nzmax) )      ! Kessler ac based on rcm/cf
! End of variables for Latin hypercube.

#ifdef SCALARS
! Variables for new mixing scheme
        allocate( sclrm(1:nzmax, 1:sclrm_dimension) )
        allocate( sclrm_forcing(1:nzmax, 1:sclrm_dimension) )

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
#endif SCALARS

        return
        end subroutine allocate_model_arrays

!------------------------------------------------------------------------
! Subroutine to deallocate variables defined in module global
!
        subroutine deallocate_model_arrays( )

        use model_arrays
        implicit none

! Prognostic variables

        deallocate( um )        ! u wind
        deallocate( vm )        ! v wind
        deallocate( thlm_forcing )
        deallocate( rtm_forcing )
        deallocate( thlm )      ! liquid potential temperature
        deallocate( rtm )       ! total water mixing ratio
        deallocate( wp2 )       ! w'^2
        deallocate( rtp2 )      ! rt'^2
        deallocate( thlp2 )     ! thl'^2
        deallocate( wprtp )     ! w'rt'
        deallocate( wpthlp )    ! w'thl'
        deallocate( rtpthlp )   ! rt'thl'
        deallocate( wp3 )       ! w'^3

! Tendency arrays for prognostic variables

        deallocate( umt )       ! u wind
        deallocate( vmt )       ! v wind

c Diagnostic variables

!        deallocate( Sc )        ! PDF width parameter
        deallocate( Scm )       ! PDF width parameter: m point
        deallocate( Sct )       ! PDF width parameter: t point
        deallocate( ug )        ! u geostrophic wind
        deallocate( vg )        ! v geostrophic wind

        deallocate( wmt )       ! imposed large scale w: t point
        deallocate( wmm )       ! imposed large scale w: m point
 
        deallocate( p )         ! pressure
        deallocate( exner )     ! exner
        deallocate( rhot )      ! density: t points
        deallocate( rhom )      ! density: m points

        deallocate( cf )        ! cloud fraction
        deallocate( thvm )      ! virtual potential temperature
        deallocate( rcm )       ! cloud water mixing ratio
        deallocate( Ncm )       ! cloud droplet number conc. ! Brian
        deallocate( sed_rcm )   ! sedimentation factor in rcm ! Brian
        deallocate( rsm )       ! saturation mixing ratio  ! Brian
        deallocate( rrm )       ! rain water mixing ratio  ! Brian
        deallocate( cond_rrm )  ! condensation factor in rrm ! Brian
        deallocate( auto_rrm )  ! autoconversion factor in rrm ! Brian
        deallocate( accr_rrm )  ! accretion factor in rrm ! Brian
        deallocate( sed_rrm )   ! sedimentation output for rrm ! Brian
        deallocate( diff_rrm )  ! diffusion output for rrm ! Brian
        deallocate( subs_rrm )  ! subsidence output for rrm ! Brian
        deallocate( Nrm )       ! rain droplet number conc.  ! Brian
        deallocate( cond_Nrm )  ! condensation factor in Nrm ! Brian
        deallocate( auto_Nrm )  ! autoconversion factor in Nrm ! Brian
        deallocate( sed_Nrm )   ! sedimentation output for Nrm ! Brian
        deallocate( diff_Nrm )  ! diffusion output for Nrm ! Brian
        deallocate( subs_Nrm )  ! subsidence output for Nrm ! Brian
        deallocate( Vr )        ! sedimentation velocity of rain ! Brian
        deallocate( mean_vol_rad ) ! rain drop mean volume radius ! Brian
        deallocate( rain_rate ) ! rainfall rate in mm/day ! Brian
        deallocate( Fprec )     ! precipitation flux in W/m^2 (momentum point) ! Brian
        deallocate( Fcsed )     ! cloud water sedimentation flux (momentum point) ! Brian
                                ! Note: units are [kg(cloud water)]/[s*m^2]

        deallocate( Frad )      ! radiative flux (momentum point)
        deallocate( Frad_LW )   ! Long_wave radiative flux (momentum point)
        deallocate( Frad_SW )   ! Short-wave radiative flux (momentum point)

        deallocate( upwp )      ! vertical u momentum flux
        deallocate( vpwp )      ! vertical v momentum flux

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
        deallocate( taum )

! Tuning Variables
#ifdef STATS
        deallocate( wp2zt )     ! w'^2 in t
        deallocate( thlp2zt )   ! th_l'^2 in t 
        deallocate( wpthlpzt )  ! w'th_l' in t
        deallocate( wprtpzt )   ! w'rt' in t
        deallocate( rtp2zt )    ! rt'^2 in t
        deallocate( rtpthlpzt ) ! rt'th_l' in t
#endif STATS

! Temporary array

!        deallocate( tmp1 )
!        deallocate( tmp2 )

! Array for pdf closure scheme

        deallocate( pdf_parms )

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
        deallocate( AKm_est )   ! Kessler ac estimate
        deallocate( AKm )       ! Exact Kessler ac
        deallocate( AKstd )     ! St dev of exact Kessler ac
        deallocate( AKstd_cld ) ! St dev of exact w/in cloud Kessler ac
        deallocate( rcm_est )   ! Monte Carlo rcm estimate
        deallocate( AKm_rcm )   ! Kessler ac based on rcm
        deallocate( AKm_rcc )   ! Kessler ac based on rcm/cf
! End of variables for Latin hypercube.

#ifdef SCALARS
! New mixing scheme structures
        deallocate( sclrm )
        deallocate( sclrm_forcing )

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
#endif SCALARS

        return
        end subroutine deallocate_model_arrays

!------------------------------------------------------------------------
! Subroutine to initialize variables defined in module global
!------------------------------------------------------------------------

        subroutine initialize_model_arrays( )

        use constants
        use model_arrays

        implicit none

! Prognostic variables

        um           = 0.      ! u wind
        vm           = 0.      ! v wind
        thlm_forcing = 0.
        rtm_forcing  = 0.
        thlm         = 0.      ! liquid potential temperature
        rtm          = 0.      ! total water mixing ratio
        wp2          = 2./3. * emin  ! w'^2
        rtp2         = 0.      ! rt'^2
        thlp2        = 0.      ! thl'^2
        wprtp        = 0.      ! w'rt'
        wpthlp       = 0.      ! w'thl'
        rtpthlp      = 0.      ! rt'thl'
        wp3          = 0.      ! w'^3

#ifdef SCALARS
! New mixing scheme Variables
        sclrm(:,:)         = 0. 
        sclrm_forcing(:,:) = 0.

        sclrprtp(:,:)      = 0.
        sclrp2(:,:)        = 0.
        sclrpthvp(:,:)     = 0. 
        sclrpthlp(:,:)     = 0.
        sclrprcp(:,:)      = 0.

        wpsclrp(:,:)       = 0. 
        wp2sclrp(:,:)      = 0. 
        wpsclrp2(:,:)      = 0. 
        wpsclrprtp(:,:)    = 0.
        wpsclrpthlp(:,:)   = 0.
#endif SCALARS

! Tendency arrays for prognostic variables

        umt = 0.               ! u wind
        vmt = 0.               ! v wind

c Diagnostic variables

!        Sc  = 0.        ! PDF width parameter
        Scm = 0.       ! PDF width parameter: m point
        Sct = 0.       ! PDF width parameter: t point
        ug  = 0.        ! u geostrophic wind
        vg  = 0.        ! v geostrophic wind

        wmt = 0.       ! imposed large scale w: t point
        wmm = 0.       ! imposed large scale w: m point
 
        p     = 0.         ! pressure
        exner = 0.     ! exner
        rhot  = 0.
        rhom  = 0.

        cf           = 0.   ! cloud fraction
        thvm         = 0.   ! virtual potential temperature
        rcm          = 0.   ! cloud water mixing ratio
        Ncm          = 0.   ! cloud droplet number conc. ! Brian
        sed_rcm      = 0.   ! sedimentation factor in rcm ! Brian
        rsm          = 0.   ! saturation mixing ratio  ! Brian
        rrm          = 0.   ! rain water mixing ratio  ! Brian
        cond_rrm     = 0.   ! condensation factor in rrm ! Brian
        auto_rrm     = 0.   ! autoconversion factor in rrm ! Brian
        accr_rrm     = 0.   ! accretion factor in rrm ! Brian
        sed_rrm      = 0.   ! sedimentation output for rrm ! Brian
        diff_rrm     = 0.   ! diffusion output for rrm ! Brian
        subs_rrm     = 0.   ! subsidence output for rrm ! Brian
        Nrm          = 0.   ! rain droplet number conc.  ! Brian
        cond_Nrm     = 0.   ! condensation factor in Nrm ! Brian
        auto_Nrm     = 0.   ! autoconversion factor in Nrm ! Brian
        sed_Nrm      = 0.   ! sedimentation output for Nrm ! Brian
        diff_Nrm     = 0.   ! diffusion output for Nrm ! Brian
        subs_Nrm     = 0.   ! subsidence output for Nrm ! Brian
        Vr           = 0.   ! sedimentation velocity of rain ! Brian
        mean_vol_rad = 0.   ! rain drop mean volume radius ! Brian
        rain_rate    = 0.   ! rainfall rate in mm/day ! Brian
        Fprec        = 0.   ! precipitation flux in W/m^2 (momentum point) ! Brian
        Fcsed        = 0.   ! cloud water sedimentation flux (momentum point) ! Brian
                            ! Note: units are [kg(cloud water)]/[s*m^2]

        Frad    = 0.
        Frad_LW = 0.   ! Long_wave radiative flux (momentum point)
        Frad_SW = 0.   ! Short-wave radiative flux (momentum point)

        upwp = 0.      ! vertical u momentum flux
        vpwp = 0.      ! vertical v momentum flux

        shear = 0.     ! wind shear production

        ! Second order moments

        wprcp   = 0.
        thlprcp = 0.
        rtprcp  = 0.

        ! Third order moments

        wpthlp2   = 0.
        wp2thlp   = 0.
        wprtp2    = 0.
        wp2rtp    = 0.
        wp2rcp    = 0.
        wprtpthlp = 0.

        ! Fourth order moments

        wp4 = 0.

        ! Buoyancy related moments

        wpthvp   = 0.
        rtpthvp  = 0.
        thlpthvp = 0.
        wp2thvp  = 0.

        Kht      = 0.
        Khm      = 0.

        em       = emin
        Lscale   = 0.
        lup      = 0.
        ldown    = 0.
        taut     = 0.
        taum     = 0.

! Temporary array

!        tmp1 = 0.
!        tmp2 = 0.

! Array fpr pdf closure scheme

        pdf_parms = 0.

! Variables for Latin hypercube microphysics.  Vince Larson 22 May 2005
        AKm_est       = 0.   ! Kessler ac estimate
        AKm           = 0.   ! Exact Kessler ac
        AKstd         = 0.   ! St dev of exact Kessler ac
        AKstd_cld     = 0.   ! St dev of exact w/in cloud Kessler ac
        rcm_est       = 0.   ! Monte Carlo rcm estimate
        AKm_rcm       = 0.   ! Kessler ac based on rcm
        AKm_rcc       = 0.   ! Kessler ac based on rcm/cf
! End of variables for Latin hypercube.

        return
        end subroutine initialize_model_arrays
!------------------------------------------------------------------------

!------------------------------------------------------------------------
! Function invalid_model_arrays()
! Checks for invalid floating point values in select model arrays.
!------------------------------------------------------------------------
        logical function invalid_model_arrays(  ) 

        use model_arrays

        implicit none

        invalid_model_arrays = .false.
           
!------------------------------------------------------------------------
!       Check whether any variable array contains a NaN for
!       um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp, 
!       wp2, & wp3.       
!------------------------------------------------------------------------

        if ( isNaN( um ) ) then
          print *, "NaN in um model array" 
!         print *, "um= ", um
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( vm ) ) then
          print *, "NaN in vm model array" 
!         print *, "vm= ", vm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( wp2 ) ) then
          print *, "NaN in wp2 model array" 
!         print *, "wp2= ", wp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( wp3 ) ) then
          print *, "NaN in wp3 model array" 
!         print *, "wp3= ", wp3
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( rtm ) ) then
          print *, "NaN in rtm model array" 
!         print *, "rtm= ", rtm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( thlm ) ) then
          print *, "NaN in thlm model array" 
!         print *, "rtm= ", rtm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( rtp2 ) ) then
          print *, "NaN in rtp2 model array" 
!         print *, "rtp2= ", rtp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( thlp2 ) ) then
          print *, "NaN in thlp2 model array" 
!         print *, "thlp2= ", thlp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( wprtp ) ) then
          print *, "NaN in wprtp model array" 
!         print *, "wprtp= ", wprtp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( wpthlp ) ) then
          print *, "NaN in wpthlp model array" 
!         print *, "wpthlp= ", wpthlp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( rtpthlp ) ) then
          print *, "NaN in rtpthlp model array" 
!         print *, "rtpthlp= ", rtpthlp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( p ) ) then
          print *, "NaN in p model array" 
!         print *, "p= ", p
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( Nrm ) ) then
          print *, "NaN in Nrm model array" 
!         print *, "Nrm= ", Nrm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( rrm ) ) then
          print *, "NaN in rrm model array" 
!         print *, "rrm= ", rrm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( wmt ) ) then
          print *, "NaN in wmt model array" 
!         print *, "wmt= ", wmt
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( wp2thvp ) ) then
          print *, "NaN in wp2thvp model array" 
!         print *, "wp2thvp = ", wp2thvp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isNaN( rtpthvp ) ) then
          print *, "NaN in rtpthvp model array" 
!         print *, "rtpthvp = ", rtpthvp
          invalid_model_arrays = .true.
        endif

        if ( isNaN( thlpthvp ) ) then
          print *, "NaN in thlpthvp model array" 
!         print *, "thlpthvp = ", thlpthvp
          invalid_model_arrays = .true.
        endif

        return
        contains
!------------------------------------------------------------------------
! function isNaN( )

! Checks if a given real is NaN, +inf or -inf

! Input: tVar: the real type array to be tested.

! Output: isNaN: true if a variable in tVar is of Not-a-Number type

! Notes: 
! Got a tip from Andy Vaught to use transfer( ) for portability.
! Ideally this would be done with intrinsic module ieee in fortran 2003
!  or with the isnan( ) function provided by some ( not pgf ) compilers
! Uses an assumed shape array, requiring it to be contained within
!   invalid_model_arrays, use an interface, or be part of the module.
!------------------------------------------------------------------------
        logical function isNaN( tVar )
        
        implicit none
      
!----- Arguments ------------------------------------------------------
        real, intent(in), dimension(:)    :: tVar

!---- Internal Variables ----------------------------------------------
        integer nanbits
        real    PosInf, NegInf

        nanbits = Z"7F800000"

        PosInf = transfer( nanbits, PosInf )
        NegInf = -( transfer( nanbits, NegInf ) )

        isNaN = .false.

!       This works on compilers with standardized floating point, because
!       the IEEE spec defines that subnormals should not equal themselves.
!       However, all compilers do not seem to follow this. 
        if ( any( (tVar /= tVar), dim=1) ) then
          isNaN = .true. 
          return
        endif

!       This a second check, assuming the above does not work as expected.
        if ( any( (tVar == PosInf), dim=1) ) then
          isNaN = .true. 
          return
        endif

!       This may never be needed, it's here just in case.
        if ( any( (tVar == NegInf), dim=1) ) then
          isNaN = .true. 
          return
        endif

        return
        end function isNaN
!------------------------------------------------------------------------

        end function invalid_model_arrays
!------------------------------------------------------------------------
