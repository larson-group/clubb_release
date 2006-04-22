!-----------------------------------------------------------------------
! $Id: hoc.F,v 1.75 2006-04-22 22:17:41 dschanen Exp $

!  MODULE HOC
!  Contains the necessary subroutines to execute individual HOC model
!  runs, using one of the driver programs (the simplest case being the
!  hoc_standalone program).
!-----------------------------------------------------------------------
       module hoc

       implicit none

!      Setup hoc_model() as the sole external interface for now

       private hoc_initialize
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
!         beta_in, gamma_coef_in, gamma_coefb_in, gamma_coefc_in, c_K_in
!         lmin_coef_in, taumin_in, taumax_in, mu_in

!         Other parameters:
!         runfile:  Namelist filename on disk to read in.
!         stdout:   Whether to print to stdout while running timesteps.

!  Output: 
!         isValid:  A false return value indicates something was unsound.
!-----------------------------------------------------------------------

       subroutine hoc_model_field_input( 
     .                       C1_in, C1b_in, C1c_in, 
     .                       C2_in, C2b_in, C2c_in, 
     .                       C4_in, C5_in, 
     .                       C6rt_in, C6rtb_in, C6rtc_in,
     .                       C6thl_in, C6thlb_in, C6thlc_in,
     .                       C7_in, C7b_in, C7c_in,
     .                       C8_in, C8b_in, C10_in, 
     .                       C11_in, C11b_in, C11c_in, C12_in,
     .                       nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .                       beta_in, gamma_coef_in, 
     .                       gamma_coefb_in, gamma_coefc_in, c_K_in, 
     .                       lmin_coef_in, taumin_in, taumax_in, mu_in,
     .                       runfile, isValid, stdout,
     .                       set_fields, sample_ratio )

       use constants, only: lmin, beta
       use diagnostic_variables
       use prognostic_variables
       use grid_class
       use numerical_check
       use hoc_parameterization_interface
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
       intrinsic mod, real, int 

!      Subroutine Arguments (constants)
       real, intent(in) :: C1_in, C1b_in, C1c_in
       real, intent(in) :: C2_in, C2b_in, C2c_in
       real, intent(in) :: C4_in, C5_in
       real, intent(in) :: C6rt_in, C6rtb_in, C6rtc_in
       real, intent(in) :: C6thl_in, C6thlb_in, C6thlc_in
       real, intent(in) :: C7_in, C7b_in, C7c_in
       real, intent(in) :: C8_in, C8b_in, C10_in
       real, intent(in) :: C11_in, C11b_in, C11c_in, C12_in
       real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in
       real, intent(in) :: beta_in
       real, intent(in) :: gamma_coef_in, gamma_coefb_in, gamma_coefc_in
       real, intent(in) :: c_K_in, taumin_in, taumax_in, lmin_coef_in
       real, intent(in) :: mu_in

!      Subroutine Arguments (Model Setting)
       character(len=*), intent(in) :: runfile ! Name of the nml file
       logical, intent(in)          :: stdout  ! whether to print output

!      Subroutine Arguments (Output)
       logical, intent(inout)       :: isValid    ! valid run? (T/F)

!      sample_ratio is used to determine often to sample data from the GrADS 
!      input per iteration of the main time stepping loop.
!      e.g. if dtmain is 20.0 seconds and the data was written every 60.0 sec, 
!      then sample ratio should be equal to 3.0.
       real, intent(in)             :: sample_ratio 

!      Internal Timing Variables

       integer ifinal

       integer niterlong

       integer sample_time

!      Coriolis parameter
       real :: fcor 

!      Flag for cloud water droplet sedimentation. - Brian
       logical :: cloud_sed

!      Flag for Khairoutdinov and Kogan rain microphysics. - Brian
       logical :: kk_rain

!      Flag for BUGsrad radiation scheme - dschanen
       logical :: lbugsrad

!      Definition of namelists
       namelist /model_setting/ runtype, nzmax, deltaz, 
     .                          day, month, year, 
     .                          time_initial, time_final, time_spinup,
     .                          dtmain, dtclosure,
     .                          sfctype, Tsfc, psfc, SE, LE, fcor,
     .                          cloud_sed, kk_rain, lbugsrad

!      Internal Loop Variables

       integer i, i1


!-------------------- Initialize the Model run --------------------------

!      Initialize run to valid
       isValid = .true.

!      Read namelist file
       open(10, file = runfile, status = 'old')
       read(10, nml = model_setting)
       close(10)

!      Allocate & initilialize variables,
!      setup grid, setup constants, and setup flags

       call parameterization_setup
     .      ( nzmax, deltaz, C1_in, C1b_in, C1c_in, 
     .        C2_in, C2b_in, C2c_in, C4_in, C5_in, 
     .        C6rt_in, C6rtb_in, C6rtc_in, C6thl_in, 
     .        C6thlb_in, C6thlc_in, C7_in, C7b_in, C7c_in, 
     .        C8_in, C8b_in, C10_in, C11_in, C11b_in, C11c_in, C12_in,
     .        nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in, 
     .        gamma_coef_in, gamma_coefb_in, gamma_coefc_in,
     .        mu_in, c_K_in, beta_in, lmin_coef_in, 
     .        taumin_in, taumax_in, fcor,
     .        lbugsrad, kk_rain, cloud_sed, .false. )

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

!      Intialize 

       time = time_initial

       call hoc_initialize( runfile, psfc, thlm, rtm, 
     .                      um, vm, ug, vg, wp2,
     .                      rcm, wmt, wmm, em, exner,
     .                      taut, taum, thvm, p,
     .                      rhot, rhom, Lscale, lup,
     .                      ldown, Kht, Khm
#ifdef SCALARS
     .                      ,sclrm, edsclrm
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
            call parameterization_timestep
     .          ( i1, niterlong, time_initial,
     .            time, dt, trim(runtype),
     .            p, rhom, rhot, exner, 
     .            thlm_forcing, rtm_forcing,
     .            um, vm, upwp, vpwp,
     .            thlm, rtm, wp2, wp3, 
     .            rtp2, thlp2, rtpthlp, 
     .            wprtp, wpthlp, Scm, taum, 
     .            rrm, Nrm,
     .            umt, vmt, isValid
#ifdef SCALARS
     .            ,sclrm, sclrm_forcing, 
     .            edsclrm, edsclrmt
#endif /*SCALARS*/
     .          )

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

        call parameterization_cleanup( .false. )

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
!         beta_in, gamma_coef_in, gamma_coefb_in, gamma_coefc_in, c_K_in
!         lmin_coef_in, taumin_in, taumax_in, mu_in

!         Other parameters:
!         runfile:  Namelist filename on disk to read in.
!         stdout:   Whether to print to stdout while running timesteps.

!  Output: 
!         isValid:  A false return value indicates something was unsound.
!-----------------------------------------------------------------------
       subroutine hoc_model_standard( C1_in, C1b_in, C1c_in,
     .                                C2_in, C2b_in, C2c_in,
     .                                C4_in, C5_in, 
     .                                C6rt_in, C6rtb_in, C6rtc_in,
     .                                C6thl_in, C6thlb_in, C6thlc_in,
     .                                C7_in, C7b_in, C7c_in,
     .                                C8_in, C8b_in, C10_in, 
     .                                C11_in, C11b_in, C11c_in, C12_in,
     .                                nu1_in, nu2_in, nu6_in, nu8_in,
     .                                nu_r_in, beta_in, gamma_coef_in,
     .                                gamma_coefb_in, gamma_coefc_in,
     .                                c_K_in, lmin_coef_in, 
     .                                taumin_in, taumax_in, mu_in,
     .                                runfile, isValid, stdout )
       use constants, only: lmin, beta
       use grid_class
       !use model_arrays
       use diagnostic_variables
       use prognostic_variables
       use numerical_check
       use hoc_parameterization_interface
#ifdef STATS
       use statistics
#endif

       implicit none

!      External
       intrinsic mod, real, int 

!      Subroutine Arguments (constants)
       real, intent(in) :: C1_in, C1b_in, C1c_in
       real, intent(in) :: C2_in, C2b_in, C2c_in
       real, intent(in) :: C4_in, C5_in
       real, intent(in) :: C6rt_in, C6rtb_in, C6rtc_in
       real, intent(in) :: C6thl_in, C6thlb_in, C6thlc_in
       real, intent(in) :: C7_in, C7b_in, C7c_in
       real, intent(in) :: C8_in, C8b_in, C10_in
       real, intent(in) :: C11_in, C11b_in, C11c_in
       real, intent(in) :: C12_in
       real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in
       real, intent(in) :: beta_in
       real, intent(in) :: gamma_coef_in,gamma_coefb_in, gamma_coefc_in 
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

!      Coriolis parameter
       real :: fcor 

!      Flag for cloud water droplet sedimentation. - Brian
       logical :: cloud_sed

!      Flag for Khairoutdinov and Kogan rain microphysics. - Brian
       logical :: kk_rain

!      Flag for BUGsrad radiation scheme - dschanen
       logical :: lbugsrad

!      Definition of namelists
       namelist /model_setting/ runtype, nzmax, deltaz, 
     .                          day, month, year, 
     .                          time_initial, time_final, time_spinup,
     .                          dtmain, dtclosure,
     .                          sfctype, Tsfc, psfc, SE, LE, fcor,
     .                          cloud_sed, kk_rain, lbugsrad

!      Internal Loop Variables

       integer i, i1

!-------------------- Initialize the Model run --------------------------

!      Initialize run to valid
       isValid = .true.

!      Read namelist file
       open(10, file = runfile, status = 'old')
       read(10, nml = model_setting)
       close(10)

!      Allocate & initilialize variables,
!      setup grid, setup constants, and setup flags

       call parameterization_setup
     .      ( nzmax, deltaz, C1_in, C1b_in, C1c_in, 
     .        C2_in, C2b_in, C2c_in, C4_in, C5_in, 
     .        C6rt_in, C6rtb_in, C6rtc_in, C6thl_in, 
     .        C6thlb_in, C6thlc_in, C7_in, C7b_in, C7c_in, 
     .        C8_in, C8b_in, C10_in, C11_in, C11b_in, C11c_in, C12_in,
     .        nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in, 
     .        gamma_coef_in, gamma_coefb_in, gamma_coefc_in,
     .        mu_in, c_K_in, beta_in, lmin_coef_in, 
     .        taumin_in, taumax_in, fcor,
     .        lbugsrad, kk_rain, cloud_sed, .false. )

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

!      Intialize 

       time = time_initial

       call hoc_initialize( runfile,  psfc, thlm, rtm, um, vm, 
     .                      ug, vg, wp2, rcm, wmt, wmm,
     .                      em, exner, taut, taum, thvm,
     .                      p, rhot, rhom, Lscale, lup,
     .                      ldown, Kht, Khm
#ifdef SCALARS
     .                      ,sclrm, edsclrm 
#endif /* SCALARS */
     .                    )

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
           call parameterization_timestep
     .          ( i1, niterlong, time_initial,
     .            time, dt, trim(runtype),
     .            p, rhom, rhot, exner, 
     .            thlm_forcing, rtm_forcing,
     .            um, vm, upwp, vpwp,
     .            thlm, rtm, wp2, wp3, 
     .            rtp2, thlp2, rtpthlp, 
     .            wprtp, wpthlp, Scm, taum, 
     .            rrm, Nrm,
     .            umt, vmt, isValid
#ifdef SCALARS
     .            ,sclrm, sclrm_forcing, 
     .            edsclrm, edsclrmt
#endif /*SCALARS*/
     .          )
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

        call parameterization_cleanup( .false. )

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
        subroutine hoc_initialize
     .             ( runfile, psfc, thlm, rtm, um, vm,
     .               ug, vg, wp2, rcm, wmt, wmm,
     .               em, exner, taut, taum, thvm,
     .               p, rhot, rhom, Lscale, lup,
     .               ldown, Kht, Khm, 
     .               sclrm, edsclrm )

        use constants
        use grid_class
        use sounding
        use surface
#ifdef STATS
        use statistics
#endif /*STATS*/
        implicit none

!       External functions
        real, external :: rsat
        real, external :: sat_rcm
        intrinsic min, max, trim, sqrt, size

!       Input
        character(len=*), intent(in) :: runfile  ! filename for the namelist

        real, intent(in) :: psfc

!       Output
        real, dimension(gr%nnzp), intent(out) :: thlm, rtm
        real, dimension(gr%nnzp), intent(out) :: um, vm 
        real, dimension(gr%nnzp), intent(out) :: ug, vg 
        real, dimension(gr%nnzp), intent(out) :: wp2
        real, dimension(gr%nnzp), intent(out) :: rcm, wmt, wmm
        real, dimension(gr%nnzp), intent(out) :: exner, em
        real, dimension(gr%nnzp), intent(out) :: p, rhot, rhom
        real, dimension(gr%nnzp), intent(out) :: Lscale, lup
        real, dimension(gr%nnzp), intent(out) :: ldown, Kht, Khm
        real, dimension(gr%nnzp), intent(out) :: taum, taut
        real, dimension(gr%nnzp), intent(out) :: thvm

!       Output (optional)
        real, optional, 
     .  dimension(gr%nnzp,sclrm_dimension), intent(out) :: sclrm
        real, optional, 
     .  dimension(gr%nnzp,sclrm_dimension), intent(out) :: edsclrm

!       Internal
        real, dimension(gr%nnzp) :: tmp1
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
        use diagnostic_variables
        use prognostic_variables
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

      end module hoc
