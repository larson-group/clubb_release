!-----------------------------------------------------------------------
! $Id: hoc.F,v 1.69 2006-02-28 18:27:28 dschanen Exp $

!  MODULE HOC
!  Contains the necessary subroutines to execute individual HOC model
!  runs, using one of the driver programs (the simplest case being the
!  hoc_standalone program).
!-----------------------------------------------------------------------
       module hoc

       implicit none

!      Setup hoc_model() as the sole external interface for now

       private latin_hypercube_sampling
       private hoc_initialize
       private hoc_closure_timestep
       private hoc_forcings_timestep
       private hoc_model_standard, hoc_model_field_input
       public hoc_model

       interface hoc_model
         module procedure hoc_model_field_input, hoc_model_standard
       end interface


!  Directives for scalars;  
!  Theta and rt are arbitrary and used to test the code.

#ifdef SCALARS
#define SCLR_THETA 1
#define SCLR_RT 2
#endif

!      Model settings

!      Grid definition
       integer, private :: nzmax
       real, private    :: deltaz

       integer, private :: day, month, year    ! day of start of simulation


       character(len=50), private :: runtype

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

!      Surface data
       real, private :: Tsfc                       ! surface temperature (K)
       real, private :: psfc                       ! surface pressure (Pa)
       real, private :: SE                         ! sensible heat flux (K/s)
       real, private :: LE                         ! latent heat flux (1/s)
       real, private :: ustar                      ! - Brian DYCOMS_2

!      Coriolis parameter
       real, private :: fcor 

!      Flag for cloud water droplet sedimentation. - Brian
       logical, private :: cloud_sed

!      Flag for Khairoutdinov and Kogan rain microphysics. - Brian
       logical, private :: kk_rain

!      Flag for BUGsrad radiation scheme - dschanen
       logical, private :: lbugsrad

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

!      Definition of namelists
       namelist /model_setting/ runtype, nzmax, deltaz, 
     .                          day, month, year, 
     .                          time_initial, time_final, time_spinup,
     .                          dtmain, dtclosure,
     .                          sfctype, Tsfc, psfc, SE, LE, fcor,
     .                          cloud_sed, kk_rain, lbugsrad
       private model_setting

       contains
!-----------------------------------------------------------------------
!  subroutine hoc_model_field_input()

!  Subprogram to integrate the pde equations for pdf closure
!  This version is modified to set fields from pre-calculated data

!  Calls:  subroutine hoc_initialize (once)
!          subroutine hoc_forcings_timestep (ifinal times)
!          subroutine hoc_closure_timestep (ifinal*niterlong times)
!          subroutine deallocate_model_arrays (once)
!          function invalid_model_arrays
!          subroutine set_fields (passed as a parameter) 

!  Input: Contants used in the model:
!         C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, C7b_in, C7c_in,
!         C8_in, C8b_in, C10_in, C11_in, C11b_in, C11_c_in, C12_in
!         nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in
!         beta_in, gamma_coef_in, c_K_in
!         lmin_coef_in, taumin_in, taumax_in, mu_in

!         Other parameters:
!         runfile:  Namelist filename on disk to read in.
!         stdout:   Whether to print to stdout while running timesteps.

!  Output: 
!         isValid:  A false return value indicates something was unsound.
!-----------------------------------------------------------------------

       subroutine hoc_model_field_input( 
     .                       C1_in, C1b_in, C1c_in, 
     .                       C2rt_in, C2thl_in, C2rtthl_in, 
     .                       C4_in, C5_in, 
     .                       C6rt_in, C6rtb_in, C6rtc_in,
     .                       C6thl_in, C6thlb_in, C6thlc_in,
     .                       C7_in, C7b_in, C7c_in,
     .                       C8_in, C8b_in, C10_in, 
     .                       C11_in, C11b_in, C11c_in, C12_in,
     .                       nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .                       beta_in, gamma_coef_in, c_K_in, 
     .                       lmin_coef_in, taumin_in, taumax_in, mu_in,
     .                       runfile, isValid, stdout,
     .                       set_fields, sample_ratio )

       use constants, only: lmin, beta
       use model_arrays
       use grid_class
#ifdef STATS
       use statistics
#endif

       implicit none

       interface
         subroutine set_fields( timestep,
     .                          um, vm, rtm, thlm, 
     .                          wp2, wprtp, wpthlp, wp3,
     .                          rtp2, thlp2, rtpthlp,
     .                          upwp, vpwp )
         implicit none

         integer, intent(in)                       :: timestep
         real, optional, dimension(:), intent(out) :: um
         real, optional, dimension(:), intent(out) :: vm
         real, optional, dimension(:), intent(out) :: rtm
         real, optional, dimension(:), intent(out) :: thlm
         real, optional, dimension(:), intent(out) :: wp2
         real, optional, dimension(:), intent(out) :: wprtp
         real, optional, dimension(:), intent(out) :: wpthlp
         real, optional, dimension(:), intent(out) :: wp3
         real, optional, dimension(:), intent(out) :: rtp2
         real, optional, dimension(:), intent(out) :: thlp2
         real, optional, dimension(:), intent(out) :: rtpthlp
         real, optional, dimension(:), intent(out) :: upwp
         real, optional, dimension(:), intent(out) :: vpwp
         end subroutine set_fields
       end interface
     
!      External
       logical, external :: invalid_model_arrays
       intrinsic mod, real, int 

!      Subroutine Arguments (constants)
       real, intent(in) :: C1_in, C1b_in, C1c_in
       real, intent(in) :: C2rt_in, C2thl_in, C2rtthl_in
       real, intent(in) :: C4_in, C5_in
       real, intent(in) :: C6rt_in, C6rtb_in, C6rtc_in
       real, intent(in) :: C6thl_in, C6thlb_in, C6thlc_in
       real, intent(in) :: C7_in, C7b_in, C7c_in
       real, intent(in) :: C8_in, C8b_in, C10_in
       real, intent(in) :: C11_in, C11b_in, C11c_in, C12_in
       real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in
       real, intent(in) :: beta_in
       real, intent(in) :: gamma_coef_in
       real, intent(in) :: c_K_in, taumin_in, taumax_in, lmin_coef_in
       real, intent(in) :: mu_in

!      Subroutine Arguments (Model Setting)
       character(len=*), intent(in) :: runfile ! Name of the nml file
       logical, intent(in)          :: stdout  ! whether to print output

!      Subroutine Arguments (Output)
       logical, intent(inout)         :: isValid    ! valid run? (T/F)

!      sample_ratio is used to determine often to sample data from the GrADS 
!      input per iteration of the main time stepping loop.
!      e.g. if dtmain is 20.0 seconds and the data was written every 60.0 sec, 
!      then sample ratio should be equal to 3.0.
       real, intent(in)             :: sample_ratio 

!      Internal Timing Variables

       integer ifinal

       integer niterlong

       integer sample_time

!      Internal Loop Variables

       integer i, i1

!-------------------- Initialize the Model run --------------------------

!      Initialize run to valid
       isValid = .true.

!      Read namelist file
       open(10, file = runfile, status = 'old')
       read(10, nml = model_setting)
       close(10)

!      Define constants
       call setup_constants( deltaz, C1_in, C1b_in, C1c_in, 
     .                       C2rt_in, C2thl_in, C2rtthl_in,
     .                       C4_in, C5_in, 
     .                       C6rt_in, C6rtb_in, C6rtc_in,
     .                       C6thl_in, C6thlb_in, C6thlc_in,
     .                       C7_in, C7b_in, C7c_in,
     .                       C8_in, C8b_in, C10_in, 
     .                       C11_in, C11b_in, C11c_in, C12_in,
     .                       nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .                       gamma_coef_in, mu_in, c_K_in,
     .                       beta_in, lmin_coef_in, 
     .                       taumin_in, taumax_in, lbugsrad )

!      Model sanity Check
!      When tuning, we will allow these invalid results
!       if (( beta < 0.0 ) .or. ( beta > 3.0 )) then
!         print *, 'beta < 0 OR beta > 3'
!         isValid = .false. 
!       endif 

       if (lmin < 4.0) then
         print *, 'lmin < 4.0'
         isValid = .false.
       endif
       if ( .not. isValid ) return

!      Allocate arrays

       call allocate_model_arrays( nzmax )

       call initialize_model_arrays( )

!      Setup grid

       call gridsetup( nzmax, deltaz )

       time = time_initial

       call hoc_initialize( thlm, rtm, um, vm, ug, vg, wp2,
     .                      rcm, wmt, wmm,
     .                      em, exner, taut, taum, thvm,
     .                      p, rhot, rhom, Lscale, lup,
     .                      ldown, Kht, Khm, upwp, vpwp,
     .                      wprtp, thlp2, rtp2, rtpthlp, wpthlp,
     .                      runfile, stdout
#ifdef SCALARS
     .                      ,sclrm, wpsclrp, sclrp2,
     .                      sclrpthlp, sclrprtp, edsclrm
#endif /* SCALARS */
     .                     )

#ifdef STATS
!      Initialize statistics output
       call stats_init( 50, runfile, gr%nnzp, gr%zt, gr%zm,
     .                  day, month, year, time, dtmain )
#endif /*STATS*/

!-------------------------- Time integration ---------------------------
!     Call hoc_closure_timestep once per each GrADS output time 
       ifinal  = ( time_final - time_initial ) / dtmain

       sample_time = 1

!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

       ! Run from iteration 1 to iteration ifinal + 1 in order to output
       ! information to the GrADS file for the final model output time.
       ! Brian.
       do i = 1, ifinal-1, 1

#ifdef STATS
!        When this time step is over, the time will be time + dtmain

         call stats_begin_timestep( time + dtmain, dtmain )
         if ( lstats_last .and. stdout ) then
           write(*,'(a,i7,a,f8.1)') 'iteration = ', 
     .                              i, '; time = ', time
         endif
#endif /*STATS*/

!        Sample from the data file
         
         if ( mod( real( i ), sample_ratio) == 0.0 ) then
           sample_time = int( i / sample_ratio ) + 1
         endif
         call set_fields( sample_time,
     .                    um, vm, rtm,           ! optional
     .                    thlm, wp2, wprtp,      !  "   "
     .                    wpthlp, wp3,           !  "   "
     .                    rtp2, thlp2, rtpthlp,  !  "   "
     .                    upwp, vpwp )           !  "   "
         isValid = .not. invalid_model_arrays( ) ! Check for bad values 
                                                 ! in the model arrays
         if ( .not. isValid ) exit

         call hoc_forcings_timestep( i, stdout, isValid )

         if ( .not. isValid ) exit

!      Compute number of iterations for closure loop
         if ( time > time_spinup ) then
            niterlong = 1
            dt        = dtmain
         else
            niterlong = int( dtmain / dtclosure )
            dt        = dtclosure
         endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         do i1=1, niterlong
           call hoc_closure_timestep( i1, niterlong, isValid )
#ifdef STATS
           call stats_end_timestep( time + dtmain, dtmain )
#endif /*STATS*/

           ! SET TIME
           ! Advance time here, not in hoc_closure_timestep,
           ! in order to facilitate use of stats.
           ! A host model, e.g. WRF, would advance time outside
           ! of hoc_closure_timestep.  Vince Larson 7 Feb 2006
           IF ( i1 < niterlong ) THEN
             time = time_initial + (i-1) * dtmain + i1 * dtclosure
           ELSEIF ( i1 == niterlong ) THEN
              time = time_initial + i * dtmain
           ENDIF

           if ( .not. isValid ) exit
         enddo

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         if ( .not. isValid ) exit
       enddo ! i=1, ifinal

!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

! Free memory

       call deallocate_model_arrays( )
#ifdef STATS
       call stats_finalize( )
#endif

       return
       end subroutine hoc_model_field_input
!-----------------------------------------------------------------------
!  subroutine hoc_model_standard()

!  Subprogram to integrate the pde equations for pdf closure.
!    This is the standard call.

!  Calls:  subroutine hoc_initialize (once)
!          subroutine hoc_forcings_timestep (ifinal times)
!          subroutine hoc_closure_timestep (ifinal*niterlong times)
!          subroutine deallocate_model_arrays (once)
!          function invalid_model_arrays

!  Input: Contants used in the model:
!         C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, C7b_in, C7c_in,
!         C8_in, C8b_in, C10_in, C11_in, C11b_in, C11c_in, C12_in
!         nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in 
!         beta_in, gamma_coef_in, c_K_in
!         lmin_coef_in, taumin_in, taumax_in, mu_in

!         Other parameters:
!         runfile:  Namelist filename on disk to read in.
!         stdout:   Whether to print to stdout while running timesteps.

!  Output: 
!         isValid:  A false return value indicates something was unsound.
!-----------------------------------------------------------------------
       subroutine hoc_model_standard( C1_in, C1b_in, C1c_in,
     .                                C2rt_in, C2thl_in, C2rtthl_in,
     .                                C4_in, C5_in, 
     .                                C6rt_in, C6rtb_in, C6rtc_in,
     .                                C6thl_in, C6thlb_in, C6thlc_in,
     .                                C7_in, C7b_in, C7c_in,
     .                                C8_in, C8b_in, C10_in, 
     .                                C11_in, C11b_in, C11c_in, C12_in,
     .                                nu1_in, nu2_in, nu6_in, nu8_in,
     .                                nu_r_in, beta_in, gamma_coef_in,
     .                                c_K_in, lmin_coef_in, 
     .                                taumin_in, taumax_in, mu_in,
     .                                runfile, isValid, stdout )
       use constants, only: lmin, beta
       use grid_class
       use model_arrays
#ifdef STATS
       use statistics
#endif

       implicit none

!      External
       logical, external :: invalid_model_arrays
       intrinsic mod, real, int 

!      Subroutine Arguments (constants)
       real, intent(in) :: C1_in, C1b_in, C1c_in
       real, intent(in) :: C2rt_in, C2thl_in, C2rtthl_in
       real, intent(in) :: C4_in, C5_in
       real, intent(in) :: C6rt_in, C6rtb_in, C6rtc_in
       real, intent(in) :: C6thl_in, C6thlb_in, C6thlc_in
       real, intent(in) :: C7_in, C7b_in, C7c_in
       real, intent(in) :: C8_in, C8b_in, C10_in
       real, intent(in) :: C11_in, C11b_in, C11c_in
       real, intent(in) :: C12_in
       real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in
       real, intent(in) :: beta_in
       real, intent(in) :: gamma_coef_in
       real, intent(in) :: c_K_in, taumin_in, taumax_in, lmin_coef_in
       real, intent(in) :: mu_in

!      Subroutine Arguments (Model Setting)
       character(len=*), intent(in) :: runfile    ! Name of the nml file
       logical, intent(in)          :: stdout     ! whether to print to the terminal
       logical, intent(inout)       :: isValid    ! " " valid run ? (T/F)

!      Internal
       integer k

!      Internal Timing Variables

       integer ifinal

       integer niterlong

!      Internal Loop Variables
       integer i, i1


!-------------------- Initialize the Model run --------------------------

!      Initialize run to valid
       isValid = .true.

!      Read namelist file
       open(10, file = runfile, status = 'old')
       read(10, nml = model_setting)
       close(10)

!      Define constants
       call setup_constants( deltaz, C1_in, C1b_in, C1c_in,
     .                       C2rt_in, C2thl_in, C2rtthl_in,
     .                       C4_in, C5_in, 
     .                       C6rt_in, C6rtb_in, C6rtc_in,
     .                       C6thl_in, C6thlb_in, C6thlc_in,
     .                       C7_in, C7b_in, C7c_in,
     .                       C8_in, C8b_in, C10_in, 
     .                       C11_in, C11b_in, C11c_in, C12_in,
     .                       nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .                       gamma_coef_in, mu_in, c_K_in,
     .                       beta_in, lmin_coef_in, 
     .                       taumin_in, taumax_in, lbugsrad )

!      Model sanity Check
       if (( beta < 0.0 ) .or. ( beta > 3.0 )) then
         print *, 'beta < 0 OR beta > 3'
         isValid = .false. 
       endif 

       if (lmin < 4.0) then
         print *, 'lmin < 4.0'
         isValid = .false.
       endif
       if ( .not. isValid ) return

!      Allocate arrays
       call allocate_model_arrays( nzmax ) 

       call initialize_model_arrays( )

!      Setup grid

       call gridsetup( nzmax, deltaz )
 
       call hoc_initialize( thlm, rtm, um, vm, ug, vg, wp2,
     .                      rcm, wmt, wmm,
     .                      em, exner, taut, taum, thvm,
     .                      p, rhot, rhom, Lscale, lup,
     .                      ldown, Kht, Khm, upwp, vpwp,
     .                      wprtp, thlp2, rtp2, rtpthlp, wpthlp,
     .                      runfile, stdout
#ifdef SCALARS
     .                      ,sclrm, wpsclrp, sclrp2,
     .                      sclrpthlp, sclrprtp, edsclrm
#endif /* SCALARS */
     .                    )

       time = time_initial

#ifdef STATS
!      Initialize statistics output
       call stats_init( 50, runfile, gr%nnzp, gr%zt, gr%zm,
     .                  day, month, year, time, dtmain )
#endif /*STATS*/

!-------------------------- Time integration ---------------------------
       ifinal  = ( time_final - time_initial )  / dtmain

!<<<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<<<

       do i = 1, ifinal, 1

#ifdef STATS
!        When this time step is over, the time will be time + dtmain

         call stats_begin_timestep( time + dtmain, dtmain )
         if ( lstats_last .and. stdout ) then
           write(*,'(a,i7,a,f8.1)') 'iteration = ', 
     .                               i, '; time = ', time
         endif
#endif /*STATS*/

         isValid = .not. invalid_model_arrays( ) ! Check for bad values 
                                                 ! in the model arrays
         if ( .not. isValid ) exit

         call hoc_forcings_timestep( i, stdout, isValid )

         if ( .not. isValid ) exit

!      Compute number of iterations for closure loop
         if ( time > time_spinup ) then
           niterlong = 1
           dt        = dtmain
         else
           niterlong = int( dtmain / dtclosure )
           dt        = dtclosure
         endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         do i1=1, niterlong
           call hoc_closure_timestep( i1, niterlong, isValid )
#ifdef STATS
           call stats_end_timestep( time + dtmain, dtmain )
#endif /*STATS*/

           ! SET TIME
           ! Advance time here, not in hoc_closure_timestep,
           ! in order to facilitate use of stats.
           ! A host model, e.g. WRF, would advance time outside
           ! of hoc_closure_timestep.  Vince Larson 7 Feb 2006
           IF ( i1 < niterlong ) THEN
             time = time_initial + (i-1) * dtmain + i1 * dtclosure
           ELSEIF ( i1 == niterlong ) THEN
              time = time_initial + i * dtmain
           ENDIF

           if ( .not. isValid ) exit
         enddo

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          if ( .not. isValid ) exit
       enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<<<

! Free memory

        call deallocate_model_arrays( )
#ifdef STATS
        call stats_finalize( )
#endif

        return
        end subroutine hoc_model_standard

!-----------------------------------------------------------------------
!   Subroutine hoc_initialize

!   Intuitively enough, this subroutine's purpose is to execute the
!   necessary steps for the initialization of the HOC model run. 

!   Calls: (all these are external)
!   subroutine read_sounding
!   subroutine hydrostatic (twice)
!   subroutine sfc_var
!   subroutine length

!   function rsat
!   function sat_rcm
!   function zt2zm (from grid_class)

!-----------------------------------------------------------------------
        subroutine hoc_initialize( thlm, rtm, um, vm, ug, vg, wp2,
     .                             rcm, wmt, wmm,
     .                             em, exner, taut, taum, thvm,
     .                             p, rhot, rhom, Lscale, lup,
     .                             ldown, Kht, Khm, upwp, vpwp,
     .                             wprtp, thlp2, rtp2, rtpthlp, wpthlp,
     .                             runfile, stdout,
     .                             sclrm, wpsclrp, sclrp2,
     .                             sclrpthlp, sclrprtp, edsclrm )

        use constants
        use grid_class
        use sounding
        use surface
#ifdef STATS
        use statistics
#endif /*STATS*/
        implicit none

!       External functions
        real, external    :: rsat
        real, external    :: sat_rcm
        intrinsic min, max, trim, sqrt

!       Output
        real, dimension(:), intent(out) :: thlm, rtm
        real, dimension(:), intent(out) :: um, vm 
        real, dimension(:), intent(out) :: ug, vg 
        real, dimension(:), intent(out) :: wp2
        real, dimension(:), intent(out) :: rcm, wmt, wmm
        real, dimension(:), intent(out) :: exner, em
        real, dimension(:), intent(out) :: p, rhot, rhom
        real, dimension(:), intent(out) :: Lscale, lup
        real, dimension(:), intent(out) :: ldown, Kht, Khm, upwp, vpwp
        real, dimension(:), intent(out) :: wprtp, thlp2, rtp2 
        real, dimension(:), intent(out) :: rtpthlp, wpthlp
        real, dimension(:), intent(out) :: taum, taut
        real, dimension(:), intent(out) :: thvm

!       Output (optional)
        real, optional, dimension(:,:), intent(out) :: sclrm
        real, optional, dimension(:,:), intent(out) :: wpsclrp, sclrp2
        real, optional, dimension(:,:), intent(out) :: sclrpthlp
        real, optional, dimension(:,:), intent(out) :: sclrprtp

        real, optional, dimension(:,:), intent(out) :: edsclrm

!       Subroutine Arguments (Model Setting)
        logical, intent(in)          :: stdout        ! whether to print to the terminal
        character(len=*), intent(in) :: runfile       ! Name of the nml file

!       Internal
        real, dimension(size( em )) :: tmp1
        integer k

! Read sounding information

        call read_sounding( thlm, rtm, um, vm, ug, vg, 
     .                      runfile, runtype
#ifdef SCALARS               
     .                      ,sclrm, edsclrm ! optional
#endif /*SCALARS*/
     .                    )


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

        select case ( trim( runtype ) )
        case ( "dycoms2_rf01", "astex_a209", "nov11_altocu",
     .         "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .         "dycoms2_rf02_nd", "dycoms2_rf02_so" )
          ! thlm profile that is initially saturated at points.
          ! thlm profile remains the same as in the input sounding.
          ! use iterative method to find initial rcm.
          do k = 1, gr%nnzp, 1
             rcm(k) = sat_rcm( thlm(k), rtm(k), p(k), exner(k) )
          enddo

        case default
          ! Initial profile is non-saturated thlm or any type of theta.
          thlm = thlm - Lv/(Cp*exner) * rcm
#ifdef SCALARS
          sclrm(:,SCLR_THETA) = sclrm(:,SCLR_THETA) - Lv/(Cp*exner)*rcm
#endif
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
        select case ( trim( runtype ) )

        ! GCSS BOMEX
        case ( "bomex" )
          do k = 1, gr%nnzp
             em(k) = 1.0 - (gr%zm(k)/3000.0)
             if ( em(k) < emin ) then
                em(k) = emin
             endif
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS ARM
        case ( "arm" )
          do k = 1, gr%nnzp
             if ( gr%zm(k) < 150.0 ) then
                em(k) = ( 0.15 * (1.0 - gr%zm(k)/150.0) ) / rhom(k)
             else
                em(k) = emin
             endif
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS FIRE Sc
        case ( "fire" ) 
          do k=1,gr%nnzp
            if ( gr%zm(k) < 700. ) then
               em(k) = 1.
            else
               em(k) = emin
            end if
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS ATEX
        case ( "atex" )
          um = max( um, -8. )

          do k = 1, gr%nnzp
             em(k) = 1.0 - (gr%zm(k)/3000.0)
             if ( em(k) < emin ) then
                em(k) = emin
             endif
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS DYCOMS II RF01
        case ( "dycoms2_rf01" ) 
          do k=1,gr%nnzp
            if ( gr%zm(k) < 800. ) then
              em(k) = 0.5
            else
              em(k) = emin
            end if
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS DYCOMS II RF02
        case ( "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .         "dycoms2_rf02_nd", "dycoms2_rf02_so" ) 
          em = 1.0

          wp2 = 2.0 / 3.0 * em

        ! Brian for Nov. 11 altocumulus case.
        case ( "nov11_altocu" )

! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
          do k=1,gr%nnzp
            if ( gr%zm(k) < 1400. ) then
               em(k) = 0.1
            else
               em(k) = emin
            end if
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)
! End Vince Larson's change.

          wp2 = 2.0 / 3.0  * em

        end select

! End Initialize TKE and other fields as needed

!       compute mixing length

        call length( thvm, thlm, rtm, rcm, 
     .               em, p, exner, Lscale, lup, ldown )

!       dissipation time

        tmp1 = SQRT( MAX( emin, zm2zt( em ) ) )
        taut = MIN( Lscale / tmp1, taumax )
        taum = MIN( ( zt2zm( Lscale ) / SQRT( em ) ), taumax )
!        taum = zt2zm( taut )

!       modification to damp noise in stable region

        do k=1,gr%nnzp
          if ( wp2(k) <= 0.005 ) then
            taut(k) = taumin
            taum(k) = taumin
          endif
        enddo

!       eddy diffusivity coefficient
!       c_K is 0.548 usually (Duynkerke and Driedonks 1987)

        Kht = c_K * Lscale * tmp1
        Khm = c_K * zt2zm( Lscale ) * SQRT( em )
!        Khm = zt2zm( Kht )

       return
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

!----------------------------------------------------------------------
        subroutine hoc_forcings_timestep( iteration, stdout, isValid )

! Modules to be included
        use constants
        use grid_class
        use model_arrays
        use gcss
        use surface

        implicit none

!       Input / Output
        integer, intent(in)        :: iteration
        logical, intent(in)        :: stdout
        logical, intent(inout)     :: isValid

!#######################################################################
!##############      FIND ALL DIAGNOSTIC VARIABLES        ##############
!#######################################################################

!-----------------------------------------------------------------------
!  For GCSS-ARM, set large-scale forcing
!-----------------------------------------------------------------------
!           print*, 'hoc.F: just before select runtype'

         select case ( runtype )
         case ( "normal" ) ! Normal run

         case ( "arm" ) ! ARM Cu case
           call arm_tndcy( time, thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing 
#endif /* SCALARS */
     .                   )

         case ( "bomex" ) ! BOMEX Cu case
           call bomex_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                       ,sclrm_forcing 
#endif /* SCALARS */
     .                     )

         case ( "fire" ) ! FIRE Sc case
           call fire_tndcy( time, wmt, wmm, thlm_forcing, 
     .                      rtm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing 
#endif /* SCALARS */
     . )
           thlm_forcing = 0.
#ifdef SCALARS
           sclrm_forcing(:, SCLR_THETA) = 0.
#endif /* SCALARS */

           call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing(:,SCLR_THETA)
#endif /* SCALARS */
     .                   )

         case ( "wangara" ) ! Wangara dry CBL
           call wangara_tndcy( time, wmt, wmm, 
     .                         thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                        ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                        )

         case ( "atex" ) ! ATEX case
           call atex_tndcy( time, time_initial, rtm, 
     .                      wmt, wmm, thlm_forcing,
     .                      rtm_forcing, isValid
#ifdef SCALARS
     .                     ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                     )

           call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing(:, SCLR_THETA)
#endif /* SCALARS */
     .                   )

         case ( "dycoms2_rf01" ) ! DYCOMS2 RF01 case
           call dycoms2_rf01_tndcy( time, rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, Frad, radht, thlm_forcing,
     .                              rtm_forcing, isValid  
#ifdef SCALARS
     .                             ,sclrm, sclrm_forcing  ! optional variables
#endif /* SCALARS */
     .                     )

         case ( "astex_a209" ) ! ASTEX Sc case for K & K
           call astex_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                             ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                     )

         case ( "dycoms2_rf02_do" ) ! DYCOMS2 RF02 case with drizzle only.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, isValid
#ifdef SCALARS
     .                              ,sclrm_forcing  ! optional variable
#endif /* SCALARS */
     .                             )

         case ( "dycoms2_rf02_ds" ) ! DYCOMS2 RF02 case with drizzle and cloud sedimentation.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, isValid
#ifdef SCALARS
     .                              ,sclrm_forcing  ! optional variable
#endif /* SCALARS */
     .                             )

         case ( "dycoms2_rf02_nd" ) ! DYCOMS2 RF02 case with no drizzle and no cloud sedimentation.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, isValid
#ifdef SCALARS
     .                              ,sclrm_forcing  ! optional variable
#endif /*SCALARS*/
     .                             )

         case ( "dycoms2_rf02_so" ) ! DYCOMS2 RF02 case with cloud water sedimentation only.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, isValid
#ifdef SCALARS
     .                             ,sclrm_forcing   ! optional variable
#endif /*SCALARS*/
     .                             )

         case ( "nov11_altocu" ) ! Nov. 11 Altocumulus case.
           call nov11_altocu_tndcy( time, time_initial, dt, thlm,
     .                              rcm, p, exner, rhot, wmt, wmm,
     .                              thlm_forcing, rtm_forcing,
     .                              Frad, Frad_LW, Frad_SW, radht
#ifdef SCALARS
     .                             ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                             )

         case default
           stop "hoc_model: Don't know how to do this runtype"

         end select

!        Bc for the second order moments

         select case ( trim( runtype ) )

         case ( "normal", "fire" )  ! Normal and FIRE
           call sfc_momentum_fluxes( um(2), vm(2),
     .                               upwp(1), vpwp(1) )
           ! sfctype = 0  fixed sfc sensible and latent heat fluxes 
           !                   as given in hoc.in
           ! sfctype = 1  bulk formula: uses given surface temperature 
           !                   and assumes over ocean
           if ( sfctype == 0 ) then
             wpthlp(1) = SE
             wprtp(1)  = LE
#ifdef SCALARS
             wpsclrp(:,SCLR_THETA) = SE
             wpsclrp(:,SCLR_RT)    = LE
#endif /* SCALARS */
           else if ( sfctype == 1 ) then
             call sfc_thermo_fluxes( um(2), vm(2),
     .                               Tsfc, psfc, 
     .                               thlm(2), rtm(2),
     .                               wpthlp(1), wprtp(1)
#ifdef SCALARS
     .                               ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                              )

           else
              write(*,*) 'Invalid value of sfctype = ', sfctype
              stop
           endif

!        For GCSS-ARM case
         case ( "arm" )
           call arm_sfclyr( time, gr%zt(2), 1.1, 
     .                      thlm(2), um(2), vm(2),
     .                      upwp(1), vpwp(1), 
     .                      wpthlp(1), wprtp(1), ustar
#ifdef SCALARS
     .                      ,wpsclrp(1,:)
#endif /* SCALARS */
     .                    )

         case ( "bomex" ) 
           call bomex_sfclyr( um(2), vm(2),
     .                        upwp(1), vpwp(1),
     .                        wpthlp(1), wprtp(1)
#ifdef SCALARS
     .                        ,wpsclrp(1,:)
#endif /* SCALARS */
     .                    )

         case ( "wangara" )
           call wangara_sfclyr( time, um(2), vm(2),
     .                           upwp(1), vpwp(1),
     .                           wpthlp(1), wprtp(1)
#ifdef SCALARS
     .                           ,wpsclrp(1,:)
#endif /* SCALARS */
     .                        )

         case ( "atex" )
           call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                       upwp(1), vpwp(1),
     .                       wpthlp(1), wprtp(1)
#ifdef SCALARS
     .                      ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                    )

         case ( "dycoms2_rf01" ) 
           call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, 
     .                               exner(1), um(2), vm(2), 
     .                               thlm(2), rtm(2), 
     .                               rhom(1), upwp(1), vpwp(1), 
     .                               wpthlp(1), wprtp(1)
#ifdef SCALARS
     .                              ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                              )
         case ( "astex_a209" )
           call astex_sfclyr( rhom(1), um(2), vm(2),
     .                        upwp(1), vpwp(1), wpthlp(1), 
     .                        wprtp(1)
#ifdef SCALARS
     .                        ,wpsclrp(1,:)
#endif /* SCALARS */
     .                      )
         case ( "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .          "dycoms2_rf02_nd", "dycoms2_rf02_so" )
           call dycoms2_rf02_sfclyr( um(2), vm(2),
     .                               upwp(1), vpwp(1),
     .                               wpthlp(1), wprtp(1), ustar
#ifdef SCALARS
     .                               ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                             )

         case ( "nov11_altocu" )
             ! There are no surface momentum or heat fluxes
             ! for the Nov. 11 Altocumulus case.

         case default
           write(*,*) 'Invalid value of runtype = ', runtype
           stop

         end select

         if ( lbugsrad ) then
           call bugsrad_hoc( deltaz*gr%nnzp, gr%nnzp, thlm, rcm, rtm,
     .                       rrm, cf, p, rhom, Tsfc,
     .                       radht, radht_SW, radht_LW,
     .                       Frad, Frad_SW, Frad_LW, thlm_forcing )
         endif

       return
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
!    subroutine update_ed_var

!    function rsat
!    function zm2zt (from grid_class)
!    function zt2zm (from grid_class)

!    If compiled with STATS:
!    subroutine stats_begin_timestep
!    subroutine stats_end_timestep
!    subroutine stats_accumulate

!    If LH_on is true:
!    permute_height_time (external)
!    latin_hypercube_sampling (Internal to module hoc)
!-----------------------------------------------------------------------
      subroutine hoc_closure_timestep( i1, niterlong, isValid )

!      Modules to be included
       use constants
       use grid_class
       use model_arrays
       use mixing
       use diagnose_variances
       use surface
       use pdf_closure
       use rain_equations, only: rain
#ifdef STATS
       use statistics
#endif /*STATS*/

       implicit none

!      External
       real, external :: rsat
       intrinsic present, sqrt, min, max

!      Input
       integer, intent(in)    :: i1        ! closure iteration number
       integer, intent(in)    :: niterlong ! total # of closure iterations
       logical, intent(inout) :: isValid

!      Internal
       integer i, j, k
       real, dimension(size( em )) :: tmp1

!      Internal for latin hypercubes
       integer, dimension(1:gr%nnzp, 
     .                    1:nt, 1:(d+1)) :: p_height_time ! matrix of rand ints

#ifdef SCALARS
       real, dimension(1:gr%nnzp,
     .                 1:sclrm_dimension) :: sclr_tmp1, sclr_tmp2,
     .                                       sclr_tmp3, sclr_tmp4
#endif /* SCALARS */


!      Set Surface variances
!      Surface variances should be set here, before the call to diag_var.
!      The reasons that surface variances can be set here are because the
!      only variables that are the input into surface variances are the
!      surface values of wpthlp, wprtp, upwp, and vpwp.  The surface values
!      of all those variables are set in the surface forcings section of the
!      GCSS cases listed in the main timestep above.  Even if they weren't
!      set there, the updates to wpthlp, wprtp, upwp, and vpwp are at the
!      end of the closure loop, right before the code loops back around to
!      this point at the top of the closure loop.
!      Surface variances need to be set here for two reasons.  One reason is
!      that the values of rtp2, thlp2, and rtpthlp at the surface will be
!      used to find the diffusional term and the mean advection term in each
!      predictive equation for those respective terms.  The other reason is
!      that if the correct surface variances are not set here and diag_var
!      outputs it's own value for them, it will results in a faulty value for
!      Scm at the surface.  Brian Griffin.  December 18, 2005.

!      Surface effects should not be included with any case where the lowest
!      level is not the ground level.  Brian Griffin.  December 22, 2005.
       IF ( TRIM( runtype ) /= "nov11_altocu" ) THEN
          call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),
     .                  wp2(1), thlp2(1), rtp2(1), rtpthlp(1)
#ifdef SCALARS
     .                  ,wpsclrp(1,:), sclrp2(1,:),
     .                  sclrprtp(1,:), sclrpthlp(1,:)
#endif /* SCALARS */
     .                )
       ELSE
          ! Nov. 11 Lowest Level Variances.
          wp2(1)     = (2.0/3.0) * emin
          thlp2(1)   = 0.0
          rtp2(1)    = 0.0
          rtpthlp(1) = 0.0
#ifdef SCALARS
          DO i = 1, sclrm_dimension, 1
             sclrp2(1,i)    = 0.0
             sclrprtp(1,i)  = 0.0
             sclrpthlp(1,i) = 0.0
          ENDDO
#endif /* SCALARS */
       ENDIF

!      Diagnose variances
!      We found that if we call diag_var first, we can use a longer timestep.
!      This requires us to store in memory Scm and taum between timesteps.
       call diag_var( taum, wmm, rtm, wprtp, 
     .                thlm, wpthlp, wp2, wp3, Scm, 
     .                rtp2, thlp2, rtpthlp,
     .                .false., dtclosure, isValid ! diag_var needs to be modified if liter=.true.
#ifdef SCALARS
     .                ,sclrm, wpsclrp, sclrp2,    ! optional vars
     .                sclrprtp, sclrpthlp         !   "       "
#endif /* SCALARS */
     .               )
       if ( .not. isValid ) return

!      Set Boundary Conditions on thlp2, rtp2, and rtpthlp output.
       thlp2(gr%nnzp)   = 0.0
       rtp2(gr%nnzp)    = 0.0
       rtpthlp(gr%nnzp) = 0.0
#ifdef SCALARS
!      Set Boundary Conditions for optional passive scalars
       do i = 1, sclrm_dimension
          sclrprtp(gr%nnzp,i)  = 0.0
          sclrp2(gr%nnzp,i)    = 0.0
          sclrpthlp(gr%nnzp,i) = 0.0
       enddo ! i = 1, sclrm_dimension
#endif /* SCALARS */


!      Compute Sc with new formula from Vince

       Scm = gamma_coef * ( 1.0 - 
     .       min(
     .       max( ( wpthlp / ( sqrt( wp2 ) * sqrt( thlp2 ) 
     .              + 0.01 * wtol * thltol ) )**2,
     .            ( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) 
     .              + 0.01 * wtol * rttol ) )**2 ),
     .       1.0      ) )

       Sct = zm2zt( Scm )


!    Latin hypercube sample generation
!    Generate p_height_time, an nnzp x nt x d array of random integers
       if ( LH_on ) then
         i_rmd = mod( i1-1, sequence_length )
         if ( i_rmd == 0) then
           call permute_height_time( gr%nnzp, nt, d+1, p_height_time )
         endif
       endif
!    End Latin hypercube sample generation

!       print*, 'hoc.F: i_rmd=', i_rmd


!      Call closure scheme

#ifdef SCALARS
!      Put closure scheme input on the m grid
       do i = 1, sclrm_dimension
          sclr_tmp1(:,i) = zm2zt( wpsclrp(:,i) ) 
          sclr_tmp2(:,i) = zm2zt( sclrprtp(:,i) ) 
          sclr_tmp3(:,i) = zm2zt( sclrp2(:,i) ) 
          sclr_tmp4(:,i) = zm2zt( sclrpthlp(:,i) ) 
       enddo ! i = 1, sclrm_dimension
#endif /* SCALARS */

       do k = 2, gr%nnzp, 1
           call pdf_closure_new( p(k), exner(k),
     .                           wmt(k), zm2zt(wp2, k), wp3(k), Sct(k),
     .                           rtm(k), zm2zt(rtp2, k), 
     .                           zm2zt( wprtp, k ),
     .                           thlm(k), zm2zt( thlp2, k ),
     .                           zm2zt( wpthlp, k ),
     .                           zm2zt(rtpthlp, k),
     .                           wp4(k), wprtp2(k), wp2rtp(k),
     .                           wpthlp2(k), wp2thlp(k), 
     .                           wprtpthlp(k),
     .                           cf(k), rcm(k), wpthvp(k), wp2thvp(k), 
     .                           rtpthvp(k), thlpthvp(k),
     .                           wprcp(k), wp2rcp(k),
     .                           rtprcp(k), thlprcp(k),
     .                           pdf_parms(k, :),
     .                           crt1, crt2, cthl1, cthl2,
     .                           Ncm(k), runtype
#ifdef SCALARS
     .                           ,sclrm(k,:), sclr_tmp1(k,:),      ! optional
     .                           wpsclrprtp(k,:), wpsclrp2(k,:),   !   "  "
     .                           sclr_tmp2(k,:), sclr_tmp3(k,:),   !   "  " 
     .                           sclr_tmp4(k,:), sclrpthvp(k,:),   !   "  "
     .                           wpsclrpthlp(k,:), sclrprcp(k,:),  !   "  "
     .                           wp2sclrp(k,:)                     !   "  "
#endif /* SCALARS */
     .                          )
!        Latin hypercube sampling.
         if ( LH_on ) call latin_hypercube_sampling( k, p_height_time )

       enddo ! k = 2, nz-1

!            print*, 'hoc.F: AKm=', AKm
!            print*, 'hoc.F: AKm_est=', AKm_est

!      Interpolate momentum variables back to momentum grid.
!      Since top momentum level is higher than top thermo level,
!      set variables at top momentum level to 0.

       wp4               = zt2zm( wp4 )
       wp4(gr%nnzp)      = 0.0
       wpthvp            = zt2zm( wpthvp )
       wpthvp(gr%nnzp)   = 0.0
       thlpthvp          = zt2zm( thlpthvp )
       thlpthvp(gr%nnzp) = 0.0
       rtpthvp           = zt2zm( rtpthvp )
       rtpthvp(gr%nnzp)  = 0.0
       wprcp             = zt2zm( wprcp )
       wprcp(gr%nnzp)    = 0.0
       rtprcp            = zt2zm( rtprcp )
       rtprcp(gr%nnzp)   = 0.0
       thlprcp           = zt2zm( thlprcp )
       thlprcp           = 0.0
#ifdef SCALARS
!      Interpolate passive scalars back onto the m grid
       do i = 1, sclrm_dimension
          sclrpthvp(:,i)       = zt2zm( sclrpthvp(:,i) )
          sclrpthvp(gr%nnzp,i) = 0.0
          sclrprcp(:,i)        = zt2zm( sclrprcp(:,i) )
          sclrprcp(gr%nnzp,i)  = 0.0
       enddo ! i=1, sclrm_dimension
#endif /* SCALARS */

!      Alternate approach for lower wp4 boundary condition:
!      call pdf scheme at the 1st momentum level, zm(1)

#ifdef SCALARS
       do i = 1, sclrm_dimension
          sclr_tmp1(:,i) = zt2zm( sclrm(:,i) )
       enddo
#endif /* SCALARS */

        call pdf_closure_new
     .      ( p(1), exner(1),
     .        wmm(1), wp2(1), zt2zm( wp3, 1 ), Scm(1),
     .        zt2zm( rtm, 1 ), rtp2(1), wprtp(1),
     .        zt2zm( thlm, 1 ), thlp2(1), wpthlp(1),
     .        rtpthlp(1),
     .        wp4(1),
     .        wprtp2(1), wp2rtp(1),
     .        wpthlp2(1), wp2thlp(1),
     .        wprtpthlp(1),
     .        cf(1), rcm(1),
     .        wpthvp(1), wp2thvp(1), 
     .        rtpthvp(1), thlpthvp(1),
     .        wprcp(1), wp2rcp(1),
     .        rtprcp(1), thlprcp(1),
     .        pdf_parms(1,:),
     .        crt1, crt2, cthl1, cthl2,
     .        Ncm(1), runtype
#ifdef SCALARS
     .        ,sclr_tmp1(1,:), wpsclrp(1,:),
     .        wpsclrprtp(1,:), wpsclrp2(1,:),
     .        sclrprtp(1,:), sclrp2(1,:),
     .        sclrpthlp(1,:), sclrpthvp(1,:),
     .        wpsclrpthlp(1,:), sclrprcp(1,:),
     .        wp2sclrp(1,:)
#endif /* SCALARS */
     .      )

!      Set other specific boundary conditions.
       wpthvp(1)         = wpthlp( 1 )
       thlpthvp(1)       = thlp2( 1 )
       wp2thvp(gr%nnzp)  = 0.0
#ifdef SCALARS
       sclrpthvp(1,SCLR_THETA) = sclrpthlp(1,SCLR_THETA)
#endif /* SCALARS */

!      Compute thvm

       thvm = thlm + ep1 * T0 * rtm + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

!      tke is assumed to be 3/2 of wp2

       em = 1.5 * wp2

!      Compute mixing length

       call length( thvm, thlm, rtm, rcm, 
     .              em, p, exner, Lscale, lup, ldown )

!      Dissipation time

       tmp1 = SQRT( MAX( emin, zm2zt( em ) ) )
       taut = MIN( Lscale / tmp1, taumax )
       taum = MIN( ( zt2zm( Lscale ) / SQRT( max(em, emin) ) ), taumax )
!       taum = zt2zm( taut )

!      Modification to damp noise in stable region
       do k = 1, gr%nnzp
         if ( wp2(k) <= 0.005 ) then
           taut(k) = taumin
           taum(k) = taumin
         endif
       enddo

!      Eddy diffusivity coefficient
!      c_K is 0.548 usually (Duynkerke and Driedonks 1987)

       Kht = c_K * Lscale * tmp1
       Khm = c_K * zt2zm( Lscale ) * SQRT( max( em, emin ) )
!       Khm = zt2zm( Kht )


!      call Khairoutdinov and Kogan (2000) scheme for rain microphysics.

       IF ( kk_rain ) THEN

          CALL rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .               rtm, rcm, Ncm, rrm, Nrm, cond_rrm,
     .               auto_rrm, accr_rrm, sed_rrm, diff_rrm,
     .               subs_rrm, cond_Nrm, auto_Nrm, sed_Nrm,
     .               diff_Nrm, subs_Nrm, Vr, mean_vol_rad,
     .               AKm_est, AKm, isValid )

          ! Rainfall rate (mm/day) should be defined on thermodynamic
          ! levels.  Brian.
          rain_rate = ( rrm * Vr ) * ( rhot / rho_lw )   ! rainfall rate (m/s)
     .                                  * ( 86400.0 * 1000.0 )

          ! Precipitation Flux (W/m^2) should be defined on
          ! momentum levels.  Brian.
          Fprec = ( zt2zm(rrm) * zt2zm(Vr) ) * ( rhom / rho_lw )
     .                                  * rho_lw * Lv

       ENDIF ! kk_rain

       if ( .not. isValid ) return

!      Store the saturation mixing ratio for output purposes.  Brian
       do k = 1, gr%nnzp, 1
         rsm(k) =
     .   rsat(p(k), (thlm(k)+(Lv/(Cp*exner(k)))*rcm(k))*exner(k))
       enddo

!#######################################################################
!############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
!#######################################################################


!      Advance rtm/wprtp and thlm/wpthlp one time step

       call timestep_mixing
     .      ( dt, taum, wmm, wmt, wp2, wp3, Scm, 
     .        rtm, wprtp, rtpthvp, rtm_forcing,
     .        thlm, wpthlp, thlpthvp, thlm_forcing, 
     .        rhot, rhom, exner, rcm, Ncm, 
     .        cond_rrm, auto_rrm, accr_rrm,             
     .        Fcsed, sed_rcm,
     .        cloud_sed, kk_rain, isValid 
#ifdef SCALARS
     .        ,sclrm, wpsclrp, sclrpthvp, sclrm_forcing   ! optional
#endif /* SCALARS */
     .      )
       if ( .not. isValid ) return

!      Set Boundary Conditions on wpthlp and wprtp output.
       wpthlp(gr%nnzp)  = 0.0
       wprtp(gr%nnzp)   = 0.0
#ifdef SCALARS
!      Set Boundary Conditions for optional passive scalars
       do i = 1, sclrm_dimension
          wpsclrp(gr%nnzp,i)   = 0.0
       enddo ! i = 1, sclrm_dimension
#endif /* SCALARS */

!      Special Nov. 11 condition for total water above cloud.
       IF ( TRIM( runtype ) == "nov11_altocu" ) THEN
          IF ( time == time_initial + 3600.0 ) THEN
             DO k = 1, gr%nnzp, 1
                IF ( gr%zt(k) > 1555.0 ) THEN
                   rtm(k) = 0.89 * rtm(k)
                ENDIF
             ENDDO
          ENDIF
       ENDIF


!      Advance wp2/wp3 one timestep

       call timestep_wp23
     .     ( dtclosure, Lscale, taum, taut, Sct, 
     .       pdf_parms(:, 13), Khm,
     .       wmt, wmm,
     .       um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .       wp2, wp3, isValid, runtype )

       if ( .not. isValid ) return

!      Set Boundary Conditions on wp2 and wp3 output.
       wp2(gr%nnzp)     = 0.0
       wp3(gr%nnzp)     = 0.0


#ifdef SCALARS
!      Eddy Diffusivity Code
       do i=1, sclrm_dimension
         edsclrmt(1:gr%nnzp,i) = - wmt * ddzm( zt2zm( edsclrm(:,i) ) )
         if ( i == SCLR_THETA) then
           call update_ed_var( edsclrm(:,i), edsclrmt(:,i), wpthlp(1), 
     .                         Khm, dtmain, wpedsclrp(:,i), isValid )
         else ! i == SCLR_RT
           call update_ed_var( edsclrm(:,i), edsclrmt(:,i), wprtp(1),
     .                         Khm, dtmain, wpedsclrp(:,i), isValid )
         end if
       end do
       ! Set boundary condition as in rt
       edsclrm(1,1:sclrm_dimension) = edsclrm(2,1:sclrm_dimension)
       if ( .not. isValid ) return
#endif /*SCALARS*/


!      Update winds

       umt
     . = - wmt * ddzm( zt2zm( um ) )
     .   + fcor * ( vm - vg )
       call update_ed_var( um, umt, upwp(1), Khm, dtmain, upwp, isValid)
       um(1)       = ( ( um(3)-um(2) )/( gr%zt(3)-gr%zt(2) ) )
     .              * ( gr%zt(1)-gr%zt(2) ) + um(2)
       um(gr%nnzp) = ( ( um(gr%nnzp-1)-um(gr%nnzp-2) )
     .                /( gr%zt(gr%nnzp-1)-gr%zt(gr%nnzp-2) ) )
     .              * ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
     .               + um(gr%nnzp-1)

       if ( .not. isValid ) return

       vmt 
     . = - wmt * ddzm( zt2zm( vm ) )
     .   - fcor * ( um - ug )
       call update_ed_var( vm, vmt, vpwp(1), Khm, dtmain, vpwp, isValid)
       vm(1)       = ( ( vm(3)-vm(2) )/( gr%zt(3)-gr%zt(2) ) )
     .              * ( gr%zt(1)-gr%zt(2) ) + vm(2)
       vm(gr%nnzp) = ( ( vm(gr%nnzp-1)-vm(gr%nnzp-2) )
     .                /( gr%zt(gr%nnzp-1)-gr%zt(gr%nnzp-2) ) )
     .              * ( gr%zt(gr%nnzp)-gr%zt(gr%nnzp-1) )
     .               + vm(gr%nnzp-1)
           
       if ( .not. isValid ) return


!      Shear production.  Brian
       do k = 1, gr%nnzp-1, 1
         shear(k) = -upwp(k) * ( um(k+1) - um(k) ) * gr%dzm(k)
     .              -vpwp(k) * ( vm(k+1) - vm(k) ) * gr%dzm(k)
       enddo
       shear(gr%nnzp) = 0.0

!#######################################################################
!#############            ACCUMULATE STATISTICS            #############
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
           sfc%x(1,ipflux) = sfc%x(1,ipflux) + Fprec(1)
           sfc%n(1,ipflux) = sfc%n(1,ipflux) + 1
         endif
       endif 

!      Added to allow tuning without using the zm file

       wp2zt     = zm2zt( wp2 )
       thlp2zt   = zm2zt( thlp2 )
       wpthlpzt  = zm2zt( wpthlp )
       wprtpzt   = zm2zt( wprtp )
       rtp2zt    = zm2zt( rtp2 )
       rtpthlpzt = zm2zt( rtpthlp )

       call stats_accumulate( time + dtmain, dtmain )

#endif /*STATS*/

       return
       end subroutine hoc_closure_timestep
!-----------------------------------------------------------------------
!       Subroutine latin_hypercube_sampling()

!       Estimate using Latin Hypercubes.  This is usually disabled by default.
!       The actual generation of a random matrix is done in a call from the
!       subroutine hoc_initialize to permute_height_time()

!       Calls:
!       subroutine lh_sampler
!       subroutine micro_calcs
!-----------------------------------------------------------------------
        subroutine latin_hypercube_sampling( k, p_height_time )
        use model_arrays

        implicit none  

!       Input 
        integer, intent(in) :: k  ! index

        integer, dimension(:,:,:), intent(in) :: p_height_time

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

        return
        end subroutine latin_hypercube_sampling
      end module hoc
