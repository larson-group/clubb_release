!-----------------------------------------------------------------------
! $Id: hoc.F,v 1.17 2005-06-23 20:07:53 dschanen Exp $

!      Module HOC
!-----------------------------------------------------------------------
       module hoc

       use outputfile_class
#ifdef STATS
       use statistics
#endif STATS

       implicit none

!      Setup hoc_model() as the sole external interface for now

       interface hoc_interface
         module procedure hoc_model, hoc_initialize, 
     .                    hoc_forcings_timestep, hoc_closure_timestep,
     .                    latin_hypercube_sampling
       end interface

       private hoc_forcings_timestep
       private latin_hypercube_sampling
       private hoc_initialize
       private hoc_closure_timestep
       public hoc_model

!      Model settings

       integer, private :: runtype     ! = 0 normal run
                                       ! = 1 GCSS ARM
                                       ! = 2 GCSS BOMEX
                                       ! = 3 GCSS FIRE
                                       ! = 4 Wangara day 33
                                       ! = 6 GCSS ATEX
                                       ! = 7 GCSS DYCOMS II (RF01)
                                       ! = 8 GCSS ASTEX with K & K changes 
                                       ! = 9 GCSS DYCOMS II (RF02)
                                       ! = 10 GCSS DYCOMS II (RF02) no drizzle
                                       ! = 11 Nov. 11 Altocumulus case.

       integer, private :: sfctype     ! = 0  fixed sfc sensible and
                                       !      latent heat fluxes as
                                       !      given in hoc.in
                                       ! = 1  bulk formula: uses
                                       !      given surface temperature
                                       !      and assumes over ocean

       real, private :: time_initial              ! time (s) of start of simulation
       real, private :: time_final                ! time (s) end of simulation
       real, private :: time_spinup               ! time (s) end of spin up period
       real, private :: time                      ! actual time of simulation

       real, private :: dtmain                    ! host model timestep (s)
       real, private :: dtclosure                 ! closure model timestep (s)
       real, private :: dt
       real, private :: dtoutput                  ! interval of file output (s)


!      Flag for Khairoutdinov and Kogan rain microphysics. - Brian
       logical                :: kk_rain

!      Surface data

       real, private :: Tsfc                       ! surface temperature (K)
       real, private :: psfc                       ! surface pressure (Pa)
       real, private :: SE                         ! sensible heat flux (K/s)
       real, private :: LE                         ! latent heat flux (1/s)
       real, private :: ustar                      ! - Brian DYCOMS_2

!     GrADS file variables

       character(len=200), private :: fdir, fname
       type (outputfile), private  :: fzt, fzm, fpdf, fsfc

!------- Local variables for Latin Hypercube sampling ------------------
       real, private    :: crt1, crt2, cthl1, cthl2   ! coeffs of s from pdf_closure_new
       integer, private :: i_rmd 

! Number of variables to sample
       integer, private, parameter :: d = 5

! n = number of calls to microphysics per timestep (normally=2)
       integer, private, parameter :: n = 12

! sequence_length = nt/n = number of timesteps before sequence repeats.
       integer, private, parameter :: sequence_length = 1

! nt = number of random samples before sequence of repeats (normally=10)
       integer, private, parameter :: nt = n * sequence_length

! Sample drawn from uniform distribution
       double precision, private :: X_u(1:n,1:(d+1))

! Sample that is transformed ultimately to normal-lognormal
       double precision, private :: X_nl(1:n,1:d)

! A true/false flag that determines whether
!     the PDF allows us to construct a sample
       logical, private :: sample_flag

!-------- End Latin hypercube section ----------------------------------

       contains

!-----------------------------------------------------------------------
!  subroutine hoc_model()

!  Subprogram to integrate the pde equations for pdf closure

!  Calls:  subroutine hoc_initialize (once)
!          subroutine hoc_forcings_timestep (ifinal times)
!          subroutine hoc_closure_timestep (ifinal*niterlong times)
!          subroutine deallocate_model_arrays (once)
!          function invalid_model_arrays

!  Input: Contants used in the model:
!         C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, C8_in, C10_in, C11_in
!         nu1_in, nu2_in, nu6_in, nu8_in, 
!         beta_in, gamma_coef_in, c_K_in
!         lmin_coef_in, taumin_in, taumax_in, mu_in

!         Other parameters:
!         runfile:  Namelist filename on disk to read in.
!         stdout:   Whether to print to stdout while running timesteps.

!  Output: 
!         isValid:  A false return value indicates something was unsound.
!-----------------------------------------------------------------------
       subroutine hoc_model( C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                       C8_in, C10_in, C11_in, 
     .                       nu1_in, nu2_in, nu6_in, nu8_in, beta_in,
     .                       gamma_coef_in, c_K_in, 
     .                       lmin_coef_in, taumin_in, taumax_in, mu_in,
     .                       runfile, isValid, stdout,
     .                       um_in, vm_in, rtm_in, thlm_in, 
     .                       wp2_in, wprtp_in, wpthlp_in, wp3_in )

       implicit none

!      External
       logical, external :: invalid_model_arrays

!      Subroutine Arguments (constants)
       real, intent(in) :: C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                     C8_in, C10_in, C11_in
       real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in
       real, intent(in) :: beta_in
       real, intent(in) :: gamma_coef_in
       real, intent(in) :: c_K_in, taumin_in, taumax_in, lmin_coef_in
       real, intent(in) :: mu_in

!      Subroutine Arguments (Model Setting)
       character(len=*), intent(in) :: runfile    ! Name of the nml file
       logical, intent(in)          :: stdout     ! whether to print to the terminal
       logical, intent(inout)       :: isValid    ! " " valid run ? (T/F)

!      Optional Input Variable arguments
       real, optional, dimension(:), intent(in) :: um_in
       real, optional, dimension(:), intent(in) :: vm_in
       real, optional, dimension(:), intent(in) :: rtm_in
       real, optional, dimension(:), intent(in) :: thlm_in
       real, optional, dimension(:), intent(in) :: wp2_in
       real, optional, dimension(:), intent(in) :: wprtp_in
       real, optional, dimension(:), intent(in) :: wpthlp_in
       real, optional, dimension(:), intent(in) :: wp3_in

!      Internal Timing Variables

       integer ifinal
       integer ioutput

       integer niterlong
       logical spunup                  ! Set but never used anywhere

!      Internal Loop Variables
       integer i, i1

! Initialize run to valid

       isValid = .true.
       
!-------------------- Initialize the Model run --------------------------
       !print*, "hoc_initialize"
       call hoc_initialize( C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                      C8_in, C10_in, C11_in, 
     .                      nu1_in, nu2_in, nu6_in, nu8_in, beta_in,
     .                      gamma_coef_in, c_K_in, 
     .                      lmin_coef_in, taumin_in, taumax_in, mu_in,
     .                      runfile, stdout ) 
       
!-------------------------- Time integration ---------------------------
       ifinal  = ( time_final - time_initial )  / dtmain
       ioutput = dtoutput / dtmain

       time    = time_initial
       spunup  = .false.

!<<<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<<<
       do i=1, ifinal
         isValid = .not. invalid_model_arrays ! Check for bad values 
                                              ! in the model arrays
        ! print*, "hoc_forcings_timestep"
         call hoc_forcings_timestep( i, stdout, isValid )
         if ( .not. isValid ) exit

!      Compute number of iterations for closure loop
         if ( time > time_spinup ) then
            spunup    = .true.
            niterlong = 1
            dt        = dtmain
         else
            spunup    = .false.
            niterlong = int( dtmain / dtclosure )
            dt        = dtclosure
         endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         do i1=1, niterlong
          ! print*, "hoc_closure_timestep"
           call hoc_closure_timestep( i1, isValid,
     .                                um_in, vm_in, rtm_in,      ! optional
     .                                thlm_in, wp2_in, wprtp_in, !  "   "
     .                                wpthlp_in, wp3_in )        !  "   "
           if ( .not. isValid ) exit
         enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          if ( .not. isValid ) exit
       enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<<<

! Free memory

        call deallocate_model_arrays( )

        return
        end subroutine hoc_model

!-----------------------------------------------------------------------
!   Subroutine hoc_initialize

!   Intuitively enough, this subroutine's purpose is to execute the
!   necessary steps for the initialization of the HOC model run. 

!   Calls: (all these are external)
!   subroutine setup_constants
!   subroutine allocate_model_arrays / initialize_model_arrays
!   subroutine gridsetup 
!   subroutine read_sounding
!   subroutine hydrostatic (twice)
!   subroutine sfc_var
!   subroutine length

!   function rsat
!   function sat_rcm
!   function zt2zm (from grid_class)

!   Fortran intrinsic functions:
!   sqrt, min, & max

!   If compiled with STATS:
!   subroutine stats_init
!   subroutine stats_begin_timestep
!   subroutine stats_end_timestep
!   subroutine stats_accumulate

!-----------------------------------------------------------------------
        subroutine hoc_initialize( C1_in, C2_in, C4_in, C5_in, C6_in, 
     .                             C7_in, C8_in, C10_in, C11_in, 
     .                             nu1_in, nu2_in, nu6_in, nu8_in, 
     .                             beta_in, gamma_coef_in, c_K_in, 
     .                             lmin_coef_in, taumin_in, taumax_in, 
     .                             mu_in, runfile, stdout ) 

        use constants
        use grid_class
        use model_arrays
#ifdef STATS
        use statistics
#endif STATS
        implicit none

!       External functions
        real, external    :: rsat
        real, external    :: sat_rcm

!       Subroutine Arguments (Constants)
        real, intent(in) :: C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                      C8_in, C10_in, C11_in
        real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in
        real, intent(in) :: beta_in
        real, intent(in) :: gamma_coef_in
        real, intent(in) :: c_K_in, taumin_in, taumax_in, lmin_coef_in
        real, intent(in) :: mu_in

!       Subroutine Arguments (Model Setting)
        logical, intent(in)          :: stdout        ! whether to print to the terminal
        character(len=*), intent(in) :: runfile       ! Name of the nml file

!       Grid definition
        integer nzmax
        real deltaz

!       Model settings
        integer day, month, year       ! day of start ot simulation

!       Internal
        integer k

!       Definition of namelists
        namelist /model_setting/ runtype, nzmax, deltaz, 
     .                           day, month, year, 
     .                           time_initial, time_final, time_spinup,
     .                           dtmain, dtclosure, dtoutput, 
     .                           sfctype, Tsfc, psfc, SE, LE, fcor

! Read namelist file
        open(10, file = runfile, status = 'old')
        read(10, nml = model_setting)
        close(10)

! Define constants
       call setup_constants( deltaz, C1_in, C2_in, C4_in, C5_in, C6_in, 
     .                       C7_in, C8_in, C10_in, C11_in,  
     .                       nu1_in, nu2_in, nu6_in, nu8_in, 
     .                       gamma_coef_in, mu_in, c_K_in,
     .                       beta_in, 
     .                       lmin_coef_in, taumin_in, taumax_in )

        call allocate_model_arrays( nzmax, n, nt, d ) 

        call initialize_model_arrays( )

! Set grid up
        call gridsetup( nzmax, deltaz, runtype )

! Read sounding information
        ! print *, "read_sounding"
        call read_sounding( thlm, rtm, um, vm, ug, vg, 
     .                      runfile, runtype )
!     .                      sclrm ) ! optional
! Debug
!        sclrm(1, :) = rtm
! End Debug

! At this point, thlm actually contains theta (except for DYCOMS). 
! We need to compute liquid water content, and initilialize thlm properly

! First, compute approximate pressure using theta
        call hydrostatic( thlm, psfc, p, exner, rhot, rhom )

! Second, use this pressure to compute liquid water
! from excess saturation

        do k = 1,gr%nnzp
           rcm(k) = max( rtm(k) - rsat(p(k),
     .                   thlm(k) * exner(k)), 0.0 )
        end do

! Compute initial theta-l

        select case (runtype)
        case ( 7, 8, 9, 10, 11 ) 
          ! thlm profile that is initially saturated at points.
          ! thlm profile remains the same as in the input sounding.
          ! use iterative method to find initial rcm.
          do k = 1, gr%nnzp, 1
             rcm(k) = sat_rcm( thlm(k), rtm(k), p(k), exner(k) )
          enddo

        case default
          ! Initial profile is non-saturated thlm or any type of theta.
          thlm = thlm - Lv/(Cp*exner) * rcm
        end select

! Now, compute initial thetav

        thvm = thlm + ep1 * T0 * rtm 
     .              + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

! Recompute more accurate initial exner function and pressure using thvm

        call hydrostatic( thvm, psfc, p, exner, rhot, rhom )

! Initialize imposed w

        do k = 1, gr%nnzp
          wmt(k) = 0.0
          wmm(k) = 0.0
        end do

! Initialize TKE and other fields as needed
        select case ( runtype )

        ! GCSS FIRE Sc
        case ( 3 ) 
          do k=1,gr%nnzp
            if ( gr%zm(k) < 700. ) then
               em(k) = 1.
            else
               em(k) = 0.
            end if
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS ATEX
        case ( 6 )
          um = max( um, -8. )

        ! GCSS DYCOMS II RF01
        case ( 7 ) 
          do k=1,gr%nnzp
            if ( gr%zm(k) < 800. ) then
              em(k) = 0.5
            else
              em(k) = 0.
            end if
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS DYCOMS II RF02
        case ( 9, 10 ) 
          em = 1.0

          wp2 = 2.0 / 3.0 * em

      ! Brian for Nov. 11 altocumulus case.
          case ( 11 )
            em = 1.0
            wp2 = (2.0/3.0) * em

        end select

! End Initialize TKE and other fields as needed

!       compute mixing length

        call length( thvm, thlm, rtm, rcm, 
     .               em, p, exner, Lscale, lup, ldown, mu )

!       dissipation time

        tmp1 = sqrt( zm2zt( em ) )
        taut = min( Lscale / tmp1, taumax )

!       modification to damp noise in stable region

        do k=1,gr%nnzp
          if ( wp2(k) <= 0.005 ) then
            taut(k) = taumin
          endif
        enddo

        taum = zt2zm( taut )

!       eddy diffusivity coefficient
!       c_K is 0.548 usually (Duynkerke and Driedonks 1987)

        Kht = c_K * Lscale * tmp1
        Khm = zt2zm( Kht )

!       surface values of second order moments

        call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),
     .                wp2(1), thlp2(1), rtp2(1), rtpthlp(1) )
!     .                sclrm(:,1), wpsclrp(:,1) )


#ifdef STATS
!       Initialize statistics output
        call stats_init( 50, runfile, gr%nnzp, gr%zt, gr%zm,
     .                   day, month, year, time, dtmain )
#endif STATS

#ifdef STATS
!      Write statistics for the initial timestep.
       call stats_begin_timestep( time, dtmain )
       if ( lstats_last .and. stdout ) then
         write(*,'(a,f8.1)') 'initial timestep; time = ', time
       endif

       call stats_accumulate( time, dtmain )
       call stats_end_timestep( time, dtmain, dtclosure )

#endif STATS
       end subroutine hoc_initialize
!----------------------------------------------------------------------
!    Subroutine hoc_forcings_timestep( )

!    Calculate tendency and surface variables

!    Calls: (* = model case)
!    subroutines *_sfclyr
!    subroutines *_tndncy
!    subroutine sfc_thermo_fluxes
!    subroutine sfc_momentum_fluxes
!    subroutine sfc_var

!    Fortran intrinsic functions:
!    none
!----------------------------------------------------------------------
        subroutine hoc_forcings_timestep( iteration, stdout, isValid )

! Modules to be included
        use constants
        use grid_class
        use model_arrays
#ifdef STATS
        use statistics, only: lstats_last
#endif STATS

        implicit none

!       Input / Output
        integer, intent(in)        :: iteration
        logical, intent(in)        :: stdout
        logical, intent(inout)     :: isValid

!       Internal 
        integer k

!       Set time
        time = time_initial + iteration * dtmain

!       Model sanity Check
        if (( beta < 0.0 ) .or. ( beta > 3.0 )) then
          print *, 'beta < 0 OR beta > 3'
          isValid = .false. 
        endif 

        if (lmin < 4.0) then
          print *, 'lmin < 4.0'
          isValid = .false.
        endif

!          print*,'hoc.F: just before iteration printout'

!          output to file every so often if called for
#ifdef STATS
        call stats_begin_timestep( time, dtmain )
        if ( lstats_last .and. stdout ) then
          write(*,'(a,i7,a,f8.1)') 'iteration = ', 
     .                              iteration, '; time = ', time
        endif
!       print*, 'hoc.F: just after iteration printout'
#endif STATS


!#######################################################################
!##############      FIND ALL DIAGNOSTIC VARIABLES        ##############
!#######################################################################

!-----------------------------------------------------------------------
!  For GCSS-ARM, set large-scale forcing
!-----------------------------------------------------------------------
!           print*, 'hoc.F: just before select runtype'

         select case ( runtype )
         case ( 0 ) ! Normal run
           kk_rain = .false.

         case ( 1 ) ! ARM Cu case
           call arm_tndcy( time, thlm_forcing, rtm_forcing )
           kk_rain = .false.

         case ( 2 ) ! BOMEX Cu case
           call bomex_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing )
           kk_rain = .false.

         case ( 3 ) ! FIRE Sc case
           call fire_tndcy( time, wmt, wmm, thlm_forcing, 
     .                      rtm_forcing)
           thlm_forcing = 0.

           call cloud_rad( rhot, rcm, Frad, thlm_forcing )
           kk_rain = .false.

         case ( 4 ) ! Wangara dry CBL
           call wangara_tndcy( time, wmt, wmm, 
     .                         thlm_forcing, rtm_forcing )
!     .                         sclrm_forcing ) ! optional variable
           kk_rain = .false.

         case ( 6 ) ! ATEX case
           call atex_tndcy( time, time_initial, rtm, 
     .                      wmt, wmm, thlm_forcing,
     .                      rtm_forcing, isValid )
           call cloud_rad( rhot, rcm, Frad, thlm_forcing )
           kk_rain = .false.

         case ( 7 ) ! DYCOMS2 RF01 case
           call dycoms2_rf01_tndcy( time, rhot, rhom, rtm, 
     .                              rcm, wmt, wmm, Frad,
     .                              thlm_forcing, rtm_forcing, 
     .                              isValid )
           kk_rain = .false.

         case ( 8 ) ! ASTEX Sc case for K & K
           call astex_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing )
           kk_rain = .true.

         case ( 9 ) ! DYCOMS2 RF02 case with drizzle
           call dycoms2_rf02_tndcy( rhot, rhom, rtm,
     .                              rcm, wmt, wmm, thlm_forcing, 
     .                              rtm_forcing, isValid ) 
           kk_rain = .true.

         case ( 10 ) ! DYCOMS2 RF02 case with no drizzle
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, 
     .                              rcm, wmt, wmm,
     .                              thlm_forcing, rtm_forcing, 
     .                              isValid ) 
           kk_rain = .false.

         case ( 11 ) ! Nov. 11 Altocumulus case.
           call nov11_altocu_tndcy( time, time_initial, dt, thlm,
     .                              rcm, p, exner, rhot, wmt, wmm,
     .                              thlm_forcing, rtm_forcing,
     .                              Frad, Frad_LW, Frad_SW )
           kk_rain = .false.

         case default
           stop "hoc_model: Don't know how to do this runtype"

         end select

!       print*, 'hoc.F: got past runtype selection'

!          Bc for the second order moments

         select case ( runtype )

         case ( 0, 3 )  ! Normal and FIRE
           call sfc_momentum_fluxes( um(2), vm(2),
     .                               upwp(1), vpwp(1) )
           ! sfctype = 0  fixed sfc sensible and latent heat fluxes 
           !                   as given in hoc.in
           ! sfctype = 1  bulk formula: uses given surface temperature 
           !                   and assumes over ocean
           if ( sfctype == 0 ) then
             wpthlp(1)     = SE
             wprtp(1)      = LE
           else if ( sfctype == 1 ) then
             call sfc_thermo_fluxes( um(2), vm(2),
     .                               Tsfc, psfc, 
     .                               thlm(2), rtm(2),
     .                               wpthlp(1), wprtp(1) )
           else
              write(*,*) 'Invalid value of sfctype = ', sfctype
              stop
           endif

!        For GCSS-ARM case
         case ( 1 )
           call arm_sfclyr( time, gr%zt(2), 1.1, 
     .                      thlm(2), um(2), vm(2),
     .                      upwp(1), vpwp(1), 
     .                      wpthlp(1), wprtp(1), ustar )

         case ( 2 ) 
           call bomex_sfclyr( um(2), vm(2),
     .                        upwp(1), vpwp(1),
     .                        wpthlp(1), wprtp(1) )

         case ( 4 )
           call wangara_sfclyr( time, um(2), vm(2),
     .                            upwp(1), vpwp(1),
     .                            wpthlp(1), wprtp(1) )
!     .                            wpsclrp(:,1) )

         case ( 6 )
           call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                       upwp(1), vpwp(1),
     .                       wpthlp(1), wprtp(1) )

         case ( 7 ) 
           call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, 
     .                               exner(1), um(2), vm(2), 
     .                               thlm(2), rtm(2), 
     .                               rhom(1), upwp(1), vpwp(1), 
     .                               wpthlp(1), wprtp(1) )
         case ( 8 )
           call astex_sfclyr( rhom(1), um(2), vm(2),
     .                        upwp(1), vpwp(1), wpthlp(1), 
     .                        wprtp(1) )
         case ( 9, 10 )
           call dycoms2_rf02_sfclyr( um(2), vm(2), 
     .                               thlm(2), rtm(2),
     .                               psfc, upwp(1), vpwp(1), 
     .                               wpthlp(1), wprtp(1) )

         case ( 11 )
             ! There are no surface momentum or heat fluxes
             ! for the Nov. 11 Altocumulus case.

         case default
           write(*,*) 'Invalid value of runtype = ', runtype
           stop

         end select

       end subroutine hoc_forcings_timestep
!-----------------------------------------------------------------------
!    subroutine hoc_closure_timestep

!    Calls:
!    subroutine diag_var
!    subroutine length
!    subroutine pdf_closure_new
!    subroutine rain (if kk_rain is true)
!    subroutine timestep_mixing
!    subroutine timestep_wp23
!    subroutine update_wind

!    function rsat
!    function zm2zt (from grid_class)
!    function zt2zm (from grid_class)

!    Fortran intrinsic functions:
!    present, sqrt, min, & max

!    If compiled with STATS:
!    subroutine stats_begin_timestep
!    subroutine stats_end_timestep
!    subroutine stats_accumulate

!    If LH_on is true:
!    permute_height_time (external)
!    latin_hypercube_sampling (Internal to module hoc)
!-----------------------------------------------------------------------
       subroutine hoc_closure_timestep( iteration, isValid,
     .                                  um_in, vm_in, rtm_in, 
     .                                  thlm_in, wp2_in, wprtp_in, 
     .                                  wpthlp_in, wp3_in )

! Modules to be included
       use constants
       use grid_class
       use model_arrays
#ifdef STATS
       use statistics
#endif STATS

       implicit none

!      External
       real, external    :: rsat

!      Input
       integer, intent(in)    :: iteration ! From i1 argument
       logical, intent(inout) :: isValid

!      Optional Input Variable arguments
       real, optional, dimension(:), intent(in) :: um_in
       real, optional, dimension(:), intent(in) :: vm_in
       real, optional, dimension(:), intent(in) :: rtm_in
       real, optional, dimension(:), intent(in) :: thlm_in
       real, optional, dimension(:), intent(in) :: wp2_in
       real, optional, dimension(:), intent(in) :: wprtp_in
       real, optional, dimension(:), intent(in) :: wpthlp_in
       real, optional, dimension(:), intent(in) :: wp3_in

!      Internal
       integer k

! Check the parameter list to see if any variables are being set manually
       if (present( um_in )) um         = um_in
       if (present( vm_in )) vm         = vm_in
       if (present( thlm_in )) thlm     = thlm_in
       if (present( wp2_in )) wp2       = wp2_in
       if (present( wpthlp_in )) wpthlp = wpthlp_in
       if (present( wp3_in )) wp3       = wp3_in
       if (present( rtm_in )) rtm       = rtm_in
       if (present( wprtp_in )) wprtp   = wprtp_in

!      Compute Sc with new formula from Vince

       Scm = gamma_coef * ( 1.0 - 
     .       min(
     .       max( ( wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) 
     .              + 0.01 * wtol * thltol ) )**2,
     .            ( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) 
     .              + 0.01 * wtol * rttol ) )**2 ),
     .       1.0      ) )

       Sct = zm2zt(Scm)

!      Compute thvm
       thvm = thlm + ep1 * T0 * rtm + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

!      Compute mixing length

       call length( thvm, thlm, rtm, rcm, 
     .              em, p, exner, Lscale, lup, ldown, mu )

!      Dissipation time

       tmp1 = sqrt( zm2zt( em ) )
       taut = min( Lscale / tmp1, taumax )

!      Modification to damp noise in stable region
       do k = 1, gr%nnzp
         if ( wp2(k) <= 0.005 ) then
           taut(k) = taumin
         endif
       enddo
!------------>For DYCOMS, we used this portion of experimental code below
!             instead of a call to length

!              call tau(thlm,thvm,rtm,rcm,cf,p,exner,taut)
!
!              Lscale = taut * tmp1
!              Lscale(1) = Lscale(2)
!              Lscale(gr%nnzp) = Lscale(gr%nnzp-1)

!<------------end of DYCOMS code

       taum = zt2zm( taut )

!      Eddy diffusivity coefficient
!      c_K is 0.548 usually (Duynkerke and Driedonks 1987)

       Kht = c_K * Lscale * tmp1
       Khm = zt2zm( Kht )

!      Diagnose variances

       call diag_var( taum, wmm, rtm, wprtp, 
     .                thlm, wpthlp, wp2, wp3, Scm, 
     .                rtp2, thlp2, rtpthlp,
     .                .false., dtclosure, isValid )
!     .                 sclrm, wpsclrp,       ! optional vars
!     .                 sclrprtp, sclrpthlp ) !   "       "

!       print*, 'hoc.F: got past diag_var'

       if ( .not. isValid ) return

!      Call closure scheme and boundary conditions

       call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),
     .               wp2(1), thlp2(1), rtp2(1), rtpthlp(1) )
!     .               wpsclrp(:,1) ) ! optional

! Latin hypercube sample generation 
! Generate p_height_time, an nnzp x nt x d array of random integers
       if ( LH_on ) then
         i_rmd = mod( iteration-1, sequence_length )
         if ( i_rmd == 0) then
           call permute_height_time( gr%nnzp, nt, d+1, p_height_time )
         endif
       endif
! End Latin hypercube generation

!       print*, 'hoc.F: i_rmd=', i_rmd

       wpthlp(gr%nnzp)  = 0.
       wprtp(gr%nnzp)   = 0.
       wp2(gr%nnzp)     = 0.
       thlp2(gr%nnzp)   = 0.
       rtp2(gr%nnzp)    = 0.
       rtpthlp(gr%nnzp) = 0.

       wp3(1)           = 0.
       wp3(2)           = 0.
       wp3(gr%nnzp)     = 0.

       do k = 2, gr%nnzp - 1
         call pdf_closure_new( p(k), exner(k), wmt(k), zm2zt(wp2, k), 
     .                         wp3(k), Sct(k), rtm(k), zm2zt(rtp2, k), 
     .                         zm2zt( wprtp, k ), thlm(k), 
     .                         zm2zt(thlp2, k), zm2zt(wpthlp, k),
     .                         wp4(k), wprtp2(k), wp2rtp(k),
     .                         wpthlp2(k), wp2thlp(k), 
     .                         zm2zt(rtpthlp, k), wprtpthlp(k),
     .                         cf(k), rcm(k), wpthvp(k), wp2thvp(k), 
     .                         rtpthvp(k), thlpthvp(k),
     .                         wprcp(k), wp2rcp(k),
     .                         rtprcp(k), thlprcp(k),
     .                         pdf_parms(k, :),
     .                         crt1, crt2, cthl1, cthl2 )

!        Latin hypercube sampling.
         if ( LH_on ) call latin_hypercube_sampling( k )
       enddo ! k = 2, nz-1

!            print*, 'hoc.F: AKm=', AKm
!            print*, 'hoc.F: AKm_est=', AKm_est
              

!      Interpolate variables back to m grid
       wp4         = zt2zm( wp4 )
       wpthvp      = zt2zm( wpthvp )
       wpthvp(1)   = wpthlp( 1 )
       thlpthvp    = zt2zm( thlpthvp )
       thlpthvp(1) = thlp2( 1 ) 
       rtpthvp     = zt2zm( rtpthvp )   
       wprcp       = zt2zm( wprcp )       
       rtprcp      = zt2zm( rtprcp )     
       thlprcp     = zt2zm( thlprcp )   

 

!      Alternate approach for lower wp4 boundary condition:
!      call pdf scheme at the 1st momentum level, zm(1)

       call pdf_closure_new
     .      ( p(1), exner(1),
     .        wmm(1), wp2(1), zt2zm( wp3, 1 ), Scm,
     .        zt2zm( rtm, 1), rtp2(1), wprtp(1),
     .        zt2zm(thlm,1), thlp2(1), wpthlp(1),
     .        wp4(1),
     .        wprtp2(1), wp2rtp(1),
     .        wpthlp2(1), wp2thlp(1),
     .        rtpthlp(1),
     .        wprtpthlp(1),
     .        cf(1), rcm(1),
     .        wpthvp(1), wp2thvp(1), 
     .        rtpthvp(1), thlpthvp(1),
     .        wprcp(1), wp2rcp(1),
     .        rtprcp(1), thlprcp(1),
     .        pdf_parms(1,:),
     .        crt1, crt2, cthl1, cthl2 )

       wp4(gr%nnzp)      = 0.
       wpthvp(gr%nnzp)   = 0.
       wp2thvp(gr%nnzp)  = 0.
       rtpthvp(gr%nnzp)  = 0.
       thlpthvp(gr%nnzp) = 0.                        

!      call Khairoutdinov and Kogan (2000) scheme for rain microphysics.

       if ( kk_rain ) then
         if ( LH_on ) then
           call rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .                rcm, rtm, rrm, Nrm, 
     .                cond_rrm, auto_rrm,
     .                accr_rrm, cond_Nrm, auto_Nrm, sed_rrm,
     .                diff_rrm, subs_rrm, Vr, mean_vol_rad,
     .                runtype, isValid, 
! LH Kessler autoconversion.  Vince Larson 22 May 2005:
     .                AKm_est, AKm )
         else 
           call rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .                rcm, rtm, rrm, Nrm, 
     .                cond_rrm, auto_rrm,
     .                accr_rrm, cond_Nrm, auto_Nrm, sed_rrm,
     .                diff_rrm, subs_rrm, Vr, mean_vol_rad,
     .                runtype, isValid )

         endif ! LH_on

         rain_rate = ( rrm * Vr ) * ( rhot / 1000.0 )     ! rain rate (m/s)
     .                                 * ( 86400.0 * 1000.0 )  ! rate: mm/day

         F_precip = ( rrm * Vr ) * ( rhot / 1000.0 )
     .                                * 1000.0 * Lv            ! Precip. flux (W/m^2)

       endif ! kk_rain

!#######################################################################
!##############         SAVE OUTPUT TO DISK               ##############
!#######################################################################

#ifdef STATS
!      Store values of surface fluxes for statistics
       if (lstats_samp) then
         if ( ish > 0 ) then
           sfc%x(1,ish) = sfc%x(1,ish) + wpthlp(1) * rhom(1) * cp
           sfc%n(1,ish) = sfc%n(1,ish) + 1
         endif
         if ( ilh > 0 ) then
           sfc%x(1,ilh) = sfc%x(1,ilh) + wprtp(1) * rhom(1) * Lv
           sfc%n(1,ilh) = sfc%n(1,ilh) + 1
         endif
         if ( iustar > 0 ) then   ! Added by Brian
           sfc%x(1,iustar) = sfc%x(1,iustar) + ustar
           sfc%n(1,iustar) = sfc%n(1,iustar) + 1
         endif
         if ( irain > 0 ) then   ! Added by Brian
           sfc%x(1,irain) = sfc%x(1,irain) + rain_rate(2)
           sfc%n(1,irain) = sfc%n(1,irain) + 1
         endif
         if ( ipflux > 0 ) then   ! Added by Brian
           sfc%x(1,ipflux) = sfc%x(1,ipflux) + F_precip(2)
           sfc%n(1,ipflux) = sfc%n(1,ipflux) + 1
         endif
       endif 
#endif STATS

! Added to allow tuning without using the zm file
#ifdef STATS
       wp2zt     = zm2zt( wp2 )
       thlp2zt   = zm2zt( thlp2 )
       wpthlpzt  = zm2zt( wpthlp )
       wprtpzt   = zm2zt( wprtp )
       rtp2zt    = zm2zt( rtp2 )
       rtpthlpzt = zm2zt( rtpthlp )
#endif STATS

#ifdef STATS
       call stats_accumulate( time, dtmain )
       call stats_end_timestep( time, dtmain, dtclosure )   ! Brian
#endif STATS

!#######################################################################
!############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
!#######################################################################

!      Advance rtm/wprtp and thlm/wpthlp one time step
       call timestep_mixing
     .      ( dt, taum, wmm, wmt, wp2, wp3, Scm, exner, 
     .        rtm, wprtp, rtpthvp, rtm_forcing,
     .        cond_rrm, auto_rrm, accr_rrm,             
     .        thlm, wpthlp, thlpthvp, thlm_forcing, isValid )
!     .        sclrm, wpsclrp, sclrpthvp, sclrm_forcing ) ! optional

!      Advance wp2/wp3 one timestep

       call timestep_wp23
     .     ( dtclosure, Lscale, taum, taut, Sct, 
     .       pdf_parms(:,13), Khm, ! Brian
     .       wmt, wmm,
     .       um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .       wp2, wp3, em, isValid)

       ! Set initial turbulence over first hour 
       ! of Nov. 11 altocumulus case.  Brian.
       select case ( runtype )
       case ( 11 )
         IF ( time < time_initial + 3600.0) THEN
           em = 1.0
           wp2 = (2.0/3.0) * em
         ENDIF
       end select

       if (present( wp2_in )) wp2       = wp2_in
       if (present( wpthlp_in )) wpthlp = wpthlp_in
       if (present( wp3_in )) wp3       = wp3_in

!      Update winds

       umt
     . = - wmt * ddzm( zt2zm( um ) )
     .   + fcor * ( vm - vg )
       call update_wind( um, umt, upwp(1), Khm, dtmain, upwp, isValid)
       um(1)       = 2. * um(2) - um(3)
       um(gr%nnzp) = 2. * um(gr%nnzp-1) - um(gr%nnzp-2 )

       if (present( um_in )) um       = um_in

       if ( .not. isValid ) return

       vmt 
     . = - wmt * ddzm( zt2zm( vm ) )
     .   - fcor * ( um - ug )
       call update_wind( vm, vmt, vpwp(1), Khm, dtmain, vpwp, isValid)
       vm(1)       = 2. * vm(2) - vm(3)
       vm(gr%nnzp) = 2. * vm(gr%nnzp-1) - vm(gr%nnzp-2 )
           
       if (present( vm_in )) vm       = vm_in

       if ( .not. isValid ) return

       ! Shear production.  Brian
       do k = 1, gr%nnzp-1, 1
         shear(k) = -upwp(k) * ( um(k+1) - um(k) ) * gr%dzm(k)
     .              - vpwp(k) * ( vm(k+1) - vm(k) ) * gr%dzm(k)
       enddo
       shear(gr%nnzp) = 0.0

       ! Store the saturation mixing ratio for output purposes.  Brian
       do k = 1, gr%nnzp, 1
         rsm(k) =
     .   rsat(p(k), (thlm(k)+(Lv/(Cp*exner(k)))*rcm(k))*exner(k))
       enddo

       return
       end subroutine hoc_closure_timestep
!-----------------------------------------------------------------------
!       Subroutine latin_hypercube_sampling()

!       Estimate using Latin Hypercubes.  This is usually disabled by default.
!       The actual generation of a random matrix is done in via a call from 
!       subroutine hoc_initialize to permute_height_time()

!       Calls:
!       subroutine lh_sampler
!       subroutine micro_calcs
!-----------------------------------------------------------------------
        subroutine latin_hypercube_sampling( k )
        use model_arrays
        implicit none  

!       Input 
        integer, intent(in) :: k  ! index

!       Internal
        integer             :: p_matrix(n, d+1)

! Choose which rows of LH sample to feed into closure.
        p_matrix(1:n,1:(d+1)) =
     .                 p_height_time( k,n*i_rmd+1:n*i_rmd+n, 1:(d+1) )

!       print*, 'hoc.F: got past p_matrix'

! Generate LH sample, represented by X_u and X_nl, for level k
        call lh_sampler( n, nt, d, p_matrix,
     .                   cf(k), pdf_parms(k, :),
     .                   crt1, crt2, cthl1, cthl2,
     .                   rrm(k), X_u, X_nl, sample_flag )

!       print *, 'hoc.F: got past lh_sampler'

! Perform LH and analytic microphysical calculations
        call micro_calcs( n, d, X_u, X_nl, sample_flag,
     .                    pdf_parms(k,:),
     .                    AKm_est(k), AKm(k), AKstd(k), AKstd_cld(k),
     .                    AKm_rcm(k), AKm_rcc(k), rcm_est(k) )

!       print*, 'k, AKm_est=', k, AKm_est(k)
!       print*, 'k, AKm=', k, AKm(k)

        end subroutine latin_hypercube_sampling
!-----------------------------------------------------------------------

        end module hoc
!------------------------------------------------------------------------
