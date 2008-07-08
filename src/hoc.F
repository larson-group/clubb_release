!-----------------------------------------------------------------------
! $Id: hoc.F,v 1.204 2008-07-08 17:36:00 faschinj Exp $

        module hoc

!       Description:
!       Contains the necessary subroutines to execute individual HOC 
!       model runs, using one of the driver programs (the simplest case
!       being the hoc_standalone program).
!-----------------------------------------------------------------------

        use stats_precision, only: time_precision ! Variable(s)

        implicit none

        ! Setup hoc_model() as the sole external interface
        private :: 
     .  hoc_initialize,
     .  hoc_forcings_timestep,
     .  hoc_restart

        public  :: hoc_model

        private ! Default to private

        ! Model settings

        ! Grid definition
        integer, private :: 
     .  nzmax,    ! Vertical extent in levels              [#]
     .  grid_type ! 1 ==> evenly-spaced grid levels
                  ! 2 ==> stretched (unevenly-spaced) grid
                  !       entered on thermodynamic grid levels;
                  !       mom. levels halfway between thermo. levels
                  !       (style of SAM stretched grid).
                  ! 3 ==> stretched (unevenly-spaced) grid
                  !       entered on momentum grid levels;
                  !       thermo. levels halfway between mom. levels
                  !       (style of WRF stretched grid).

        real, private :: 
     .  deltaz, ! Change per grid level                 [m]
     .  zm_init ! Initial point on the momentum grid    [m]

!$omp   threadprivate(nzmax, grid_type, zm_init, deltaz)

        ! For grid_type 2 or 3 (stretched grid cases)
        character(len=100), private ::
     .  zt_grid_fname, ! Path and filename of thermodynamic level altitudes
     .  zm_grid_fname  ! Path and filename of momentum level altitudes

!$omp   threadprivate(zt_grid_fname, zm_grid_fname)

        integer, private :: 
     .  day, month, year ! Day of start of simulation

!$omp   threadprivate(day, month, year)

        real, private :: 
     .  rlat, ! Latitude  [Degrees North]
     .  rlon  ! Longitude [Degrees East]

!$omp threadprivate(rlat, rlon)

        character(len=50), private :: 
     .  runtype ! String identifying the model case; e.g. bomex
!$omp   threadprivate(runtype)

        ! 0: fixed sfc sensible and latent heat fluxes as
        !    given in namelist
        ! 1: bulk formula: uses given surface temperature
        !    and assumes over ocean
        integer, private :: sfctype 
!$omp   threadprivate(sfctype)

        real(kind=time_precision), private ::
     .  time_initial, ! Time of start of simulation     [s]
     .  time_final,   ! Time end of simulation          [s]
     .  time_spinup,  ! Time end of spin up period      [s]
     .  time_current  ! Current time of simulation      [s]
!$omp   threadprivate(time_initial, time_final, time_spinup,
!$omp.                time_current)

        real(kind=time_precision), private :: 
     .  dtmain,     ! Main model timestep                      [s]
     .  dtclosure,  ! Closure model timestep                   [s]
     .  dt          ! Current model timestep (based on spinup) [s]
!$omp   threadprivate(dtmain, dtclosure, dt)

        contains

!-----------------------------------------------------------------------
        subroutine hoc_model
     .            ( params, runfile, err_code, stdout, linput_fields )
!       Description:
!       Subprogram to integrate the pde equations for pdf closure.
!       This is the standard call.

!       Calls:  subroutine hoc_initialize (once)
!          subroutine hoc_forcings_timestep (ifinal times)
!          subroutine hoc_closure_timestep (ifinal*niterlong times)
!          subroutine deallocate_model_arrays (once)
!          function invalid_model_arrays
!          subroutine set_fields (passed as a parameter) 

!       Output: 
!         err_code:  An error code is returned indicating the status of the
!         model. See error_code.F.

!-----------------------------------------------------------------------

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: read_grid_heights ! Procedure(s)

        use param_index, only: nparams ! Variable(s)

        use array_index, only: iisclr_rt, iisclr_thl, iiCO2 ! Variables

        use diagnostic_variables, only: ug, vg, hydromet, em, ! Variable(s)
     .                            taut, thvm, Lscale, Kht, Khm,
     .                            um_ref, vm_ref

        use prognostic_variables, only: 
     .      Tsfc, psfc, SE, LE, thlm, rtm,    ! Variable(s)
     .      um, vm, wp2, rcm, wmt, wmm, exner,
     .      taum, p, rhom, upwp, vpwp, wpthlp,
     .      rhot, wprtp, wpthlp_sfc, wprtp_sfc,
     .      upwp_sfc, vpwp_sfc, thlm_forcing,
     .      rtm_forcing, up2, vp2, wp3, rtp2,
     .      thlp2, rtpthlp, Scm, cf

        use prognostic_variables, only: 
     .      sclrm, edsclrm, wpsclrp_sfc, ! Variables
     .      wpedsclrp_sfc, sclrm_forcing, wpsclrp

        use numerical_check, only: invalid_model_arrays ! Procedure(s)

        use inputfields, only: compute_timestep, grads_fields_reader ! Procedure(s)
     
        use inputfields, only: datafilet ! Variable(s)

        use hoc_parameterization_interface, only:
     .                                 parameterization_setup, ! Procedure(s) 
     .                                 parameterization_cleanup,
     .                                 parameterization_timestep

        use constants, only: fstdout, fstderr ! Variable(s)

        use error_code, only: clubb_var_out_of_bounds, ! Variable(s)
     .                        clubb_var_equals_NaN,
     .                        clubb_rtm_level_not_found

        use error_code, only: fatal_error, ! Procedure(s)
     .                        set_clubb_debug_level
        
        use stats_precision, only: time_precision ! Variable(s)

        use array_index, only: iisclr_rt, iisclr_thl, iiCO2 ! Variables

        use microphys_driver, only: init_microphys ! Subroutine

#ifdef STATS
        use stats_variables, only: lstats_last ! Variable(s

        use stats_subs, only: 
     .      stats_begin_timestep, stats_end_timestep, ! Procedure(s)
     .      stats_finalize, stats_init

#endif /*STATS*/

        implicit none

        ! Because Fortran I/O is not thread safe, we use this here to
        ! insure that no model uses the same file number simultaneously
        ! when doing a tuning run. -dschanen 31 Jan 2007
#ifdef _OPENMP
        integer :: omp_get_thread_num
#endif
        ! External
        intrinsic :: mod, real, int 

        ! Constant parameters
        integer, parameter :: sclr_max = 10

        ! Input Variables
        logical, intent(in) :: 
     .  stdout,       ! Whether to print output per timestep
     .  linput_fields ! Whether to set model variables from a file

        real, intent(in), dimension(nparams) :: 
     .  params  ! Model parameters, C1, nu2, etc.

        ! Subroutine Arguments (Model Setting)
        character(len=*), intent(in) :: 
     .  runfile ! Name of file containing &model_setting and &sounding

        ! Output Variables
        integer, intent(inout) :: err_code    ! valid run?

        ! Local Variables
        ! Internal Timing Variables
        integer ::
     .  ifinal,
     .  niterlong

        integer ::
     .  debug_level     ! Amount of debugging information

        real :: 
     .  fcor,           ! Coriolis parameter            [s^-1]
     .  T0,             ! Reference Temperature         [K]
     .  ts_nudge        ! Timescale for u/v nudging     [s]

        real, dimension(sclr_max) ::
     .  sclr_tol        ! Thresholds on the passive scalars     [units vary]

        real(kind=time_precision) ::
     .  time_restart    ! Time of model restart run     [s]

        logical :: 
     .  cloud_sed,    ! Flag for cloud water droplet sedimentation. - Brian
     .  kk_rain,      ! Flag for Khairoutdinov and Kogan rain microphysics. - Brian
     .  licedfs,      ! Flag for simplified ice scheme
     .  lcoamps_micro,! Flag for COAMPS microphysical scheme
     .  lbugsrad,     ! Flag for BUGsrad radiation scheme
     .  luv_nudge,    ! Whether to adjust the winds within the timestep
     .  lrestart,     ! Flag for restarting from GrADS file
     .  lKhm_aniso    ! Whether to use anisotropic Khm.  - Michael Falk 2 Feb 2007

        character(len=50) :: 
     .  restart_path_case ! GRADS file used in case of restart

        logical ::
     .  lstats ! Whether statistics are computed and output to disk

        character(len=10) ::
     .  stats_fmt  ! File format for stats; typically GrADS.

        character(len=100) ::
     .  fname_prefix ! Prefix of stats filenames, to be followed by, for example "_zt"

        real(kind=time_precision) ::
     .  stats_tsamp,  ! Stats sampling interval [s]
     .  stats_tout    ! Stats output interval   [s]

        ! Grid altitude arrays
        real, dimension(:), allocatable :: 
     .  momentum_heights, thermodynamic_heights

        ! Dummy dx and dy horizontal grid spacing.
        real :: dummy_dx, dummy_dy

        integer :: i, i1 ! Internal Loop Variables
        integer :: iinit ! initial iteration

        integer :: 
     .  iunit,          ! File unit used for I/O
     .  hydromet_dim,   ! Number of hydrometeor species
     .  sclr_dim        ! Number of passive scalars

        integer :: itime_nearest ! Used for and inputfields run [s]

        ! Definition of namelists
        namelist /model_setting/ 
     .  runtype, nzmax, grid_type, deltaz, zm_init,
     .  zt_grid_fname, zm_grid_fname, 
     .  day, month, year, rlat, rlon,
     .  time_initial, time_final, time_spinup,
     .  dtmain, dtclosure,
     .  sfctype, Tsfc, psfc, SE, LE, fcor, T0, ts_nudge,
     .  cloud_sed, kk_rain, licedfs, lcoamps_micro, 
     .  lbugsrad, lKhm_aniso, luv_nudge, lrestart, restart_path_case,
     .  time_restart, debug_level,
     .  sclr_tol,
     .  sclr_dim, iisclr_thl, iisclr_rt, iiCO2 

         namelist /stats_setting/
     .   lstats, fname_prefix, stats_tsamp, stats_tout, stats_fmt 

!-----------------------------------------------------------------------

        ! Initialize the model run 

        ! Pick some default values for model_setting
        runtype   = "generic"
        nzmax     = 100
        grid_type = 1
        deltaz    = 40.
        zm_init   = 0.
        zt_grid_fname = ''
        zm_grid_fname = ''

        day   = 1
        month = 1
        year  = 1900

        rlat = 0.
        rlon = 0.

        time_initial = 0.
        time_final   = 3600.
        time_spinup  = 0.

        dtmain    = 30.
        dtclosure = 30.

        sfctype  = 0
        tsfc     = 288.
        psfc     = 1000.e2
        SE       = 0.
        LE       = 0.
        fcor     = 1.e-4
        T0       = 300.
        ts_nudge = 86400.

        cloud_sed     = .false.
        kk_rain       = .false.
        licedfs       = .false.
        lcoamps_micro = .false.
        lbugsrad      = .false.
        lKhm_aniso    = .false.
        luv_nudge     = .false.
        lrestart      = .false.
        restart_path_case = "none"
        time_restart  = 0.
        debug_level   = 2

        sclr_dim  = 0
        iisclr_thl = -1
        iisclr_rt  = -1
        iiCO2 = -1

        sclr_tol(1:sclr_max) = 1.e-2

        ! Pick some default values for stats_setting
        lstats       = .false.
        fname_prefix = ""
        stats_fmt    = ""
        stats_tsamp  = 0.
        stats_tout   = 0.

        ! Figure out which iounit to use
#ifdef _OPENMP
        iunit = omp_get_thread_num( ) + 10
#else
        iunit = 10
#endif

        ! Read namelist file
        open(unit=iunit, file=runfile, status='old')
        read(unit=iunit, nml=model_setting)
        read(unit=iunit, nml=stats_setting)
        close(unit=iunit)

        call set_clubb_debug_level( debug_level )

        ! Allocate stretched grid altitude arrays.
        allocate( momentum_heights(1:nzmax), 
     .            thermodynamic_heights(1:nzmax) )

        ! Handle the reading of grid altitudes for 
        ! stretched (unevenly-spaced) grid options.
        ! Do some simple error checking for all grid options.
        call read_grid_heights( nzmax, grid_type,
     .                          zm_grid_fname, zt_grid_fname,
     .                          momentum_heights,
     .                          thermodynamic_heights )

        ! Dummy horizontal grid spacing variables.
        dummy_dx = 0.0
        dummy_dy = 0.0

        ! Setup microphysical fields
        call init_microphys
     .       ( kk_rain, lcoamps_micro, licedfs, hydromet_dim )

        ! Allocate & initialize variables,
        ! setup grid, setup constants, and setup flags

        call parameterization_setup
     .       ( nzmax, T0, ts_nudge, hydromet_dim, sclr_dim, 
     .         sclr_tol(1:sclr_dim), params,
     .         lbugsrad, kk_rain, licedfs, lcoamps_micro,
     .         cloud_sed, luv_nudge, lKhm_aniso,
     .         .false., grid_type, deltaz, zm_init,
     .         momentum_heights, thermodynamic_heights,
     .         dummy_dx, dummy_dy, err_code )
      

        if ( err_code == clubb_var_out_of_bounds ) return


        ! Deallocate stretched grid altitude arrays
        deallocate( momentum_heights, thermodynamic_heights )

        if ( .not. lrestart ) then

           time_current = time_initial
           iinit = 1

           call hoc_initialize( iunit, runfile, psfc, thlm, rtm, 
     .                          um, vm, ug, vg, wp2,
     .                          rcm, hydromet,
     .                          wmt, wmm, em, exner,
     .                          taut, taum, thvm, p,
     .                          rhot, rhom, Lscale,
     .                          Kht, Khm, um_ref, vm_ref,
     .                          sclrm, edsclrm )

        else  ! restart
               
           ! Joshua Fasching March 2008
           time_current = time_restart + dtmain
           ! time_current = time_restart

           ! Determining what iteration to restart at.
           ! The value is increased by 1 to sychronize with restart data.
           ! Joshua Fasching February 2008
          
           ! Ensure that iteration num, iinit, is an integer, so that model time is 
           !   incremented correctly by iteration number at end of timestep
           if ( mod( (time_restart-time_initial) , dtmain ) /= 0 ) then
              print*, "Error: (time_restart-time_initial) ", 
     .                "is not a multiple of dtmain."
              print*, "time_restart = ", time_restart
              print*, "time_initial = ", time_initial
              print*, "dtmain = ", dtmain
              stop
           end if
 
           iinit = floor( ( time_current - time_initial ) / dtmain ) + 1

           call hoc_restart( iunit, runfile, restart_path_case, 
     .                       trim( fname_prefix )//"_zt", time_restart, 
     .                       thlm, rtm, um, vm,
     .                       ug, vg, upwp, vpwp, wmt, wmm, 
     .                       um_ref, vm_ref, wpthlp, wprtp,
     .                       wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .                       sclrm, edsclrm )
        end if ! ~lrestart
#ifdef STATS

#ifdef _OPENMP
        iunit = omp_get_thread_num( ) + 50
#else
        iunit = 50
#endif
        ! Initialize statistics output
        call stats_init
     .       ( iunit, fname_prefix, 
     .         lstats, stats_fmt, stats_tsamp, stats_tout,
     .         runfile, gr%nnzp, gr%zt, gr%zm,
     .         day, month, year, rlat, rlon, time_current,
     .         dtmain )
#endif /*STATS*/

        ! Time integration
        ! Call hoc_closure_timestep once per each GrADS output time 
        ifinal  = floor(( time_final - time_initial ) / dtmain)

      
!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

        do i = iinit, ifinal, 1

#ifdef STATS
          ! When this time step is over, the time will be time + dtmain

          ! We use elapsed time for stats_begin_step
          if (.not. lrestart) then    
             call stats_begin_timestep
     .            ( time_current-time_initial+dtmain, dtmain )
          else
             ! Different elapsed time for restart
             ! Joshua Fasching March 2008     
             call stats_begin_timestep
     .            ( time_current-time_restart, dtmain )
          end if

#endif /*STATS*/

          ! If we're doing an inputfields run, get the values for our
          ! model arrays from a GrADS file
          if ( linput_fields ) then
            call compute_timestep( iunit, datafilet, .false.,
     .                             time_current, itime_nearest )
            call grads_fields_reader( max( itime_nearest, 1 ) )
          end if

          if ( invalid_model_arrays( ) ) then
              err_code = clubb_var_equals_NaN ! Check for bad values 
                                              ! in the model arrays
          end if

          if ( err_code == clubb_var_equals_NaN ) exit

          call hoc_forcings_timestep( dtmain, i, stdout, err_code )

          if ( err_code == clubb_rtm_level_not_found ) exit
         
          ! Compute number of iterations for closure loop
          if ( time_current > time_spinup ) then
            niterlong = 1
            dt        = dtmain
          else
            niterlong = floor( dtmain / dtclosure )
            dt        = dtclosure
          end if

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          do i1=1, niterlong
            call parameterization_timestep
     .           ( i, dt, fcor,
     .             thlm_forcing, rtm_forcing, wmm, wmt,
     .             wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .             p, rhom, rhot, exner,
     .             um, vm, upwp, vpwp, up2, vp2,
     .             thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .             rtp2, thlp2, rtpthlp,
     .             Scm, taum, rcm, cf,
     .             err_code, .false.,
     .             wpsclrp_sfc, wpedsclrp_sfc, 
     .             sclrm, sclrm_forcing, edsclrm,
     .             wpsclrp )
           
#ifdef STATS
            call stats_end_timestep( time_current + dtmain, dtmain )
#endif /*STATS*/

            ! Set Time
            ! Advance time here, not in parameterization_timestep,
            ! in order to facilitate use of stats.
            ! A host model, e.g. WRF, would advance time outside
            ! of hoc_closure_timestep.  Vince Larson 7 Feb 2006
            if ( i1 < niterlong ) then
              time_current = time_initial + (i-1) * dtmain 
     .                     + i1 * dtclosure
            else if ( i1 == niterlong ) then
              time_current = time_initial + i * dtmain
            end if
#ifdef STATS
            ! This was moved from above to be less confusing to the user,
            ! since before it would appear as though the last timestep
            ! was not executed. -dschanen 19 May 08
            if ( lstats_last .and. stdout ) then
              write(unit=fstdout,fmt='(a,i8,a,f10.1)') 'iteration = ', 
     .          i, '; time = ', time_current
            end if
#endif


            if ( fatal_error( err_code ) ) exit
           
          end do ! i1=1..niterlong

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
          if ( fatal_error( err_code ) ) exit
          
        end do ! i=1, ifinal

!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

! Free memory

        call parameterization_cleanup( .false. )

#ifdef STATS
        call stats_finalize( )
#endif

        return
        end subroutine hoc_model

!-----------------------------------------------------------------------
        subroutine hoc_initialize
     .             ( iunit, runfile, psfc, thlm, rtm, um, vm,
     .               ug, vg, wp2, rcm, hydromet,
     .               wmt, wmm, em, exner, taut, taum, thvm,
     .               p, rhot, rhom, Lscale,
     .               Kht, Khm, um_ref, vm_ref,
     .               sclrm, edsclrm )

!       Description:
!       Execute the necessary steps for the initialization of the 
!       HOC model run. 

!       Calls: ( all these are external except compute_length, rsat,
!       stat_rcm )
!       subroutine read_sounding
!       subroutine hydrostatic (twice)
!       subroutine compute_length

!       function rsat
!       function sat_rcm
!       function zt2zm (from grid_class)
!       function zm2zt (from grid_class)
!-----------------------------------------------------------------------

        use constants, only: 
     .      Cp, Lv, ep2, ep1, emin, wtol ! Variable(s)

        use parameters, only: 
     .      T0, taumax, taumin, c_K, sclr_dim, hydromet_dim  ! Variable(s)

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zm2zt, zt2zm ! Procedure(s)

        use sounding, only: read_sounding ! Procedure(s)

        use model_flags, only: luv_nudge ! Variable(s)

        use arm_0003, only: arm_0003_init ! Procedure(s)

        use arm_3year, only: arm_3year_init ! Procedure(s)

        use arm_97, only: arm_97_init ! Procedure(s)

        use mpace_a, only: mpace_a_init ! Procedure(s)

        use lba, only: lba_init ! Procedure(s)

        use mixing_length, only: compute_length ! Procedure(s)

        use error_code, only: clubb_no_error ! Variable(s)

        use array_index, only: iisclr_thl ! Variable(s)
        
        ! Joshua Fasching
        ! March 2008
        use saturation, only: rsat, sat_rcm ! Procedure(s)
       
        use hydrostatic_mod, only: hydrostatic ! Procedure(s)

        implicit none

        intrinsic :: min, max, trim, sqrt, size

        ! Input
        integer, intent(in) :: iunit
        character(len=*), intent(in) :: runfile  ! filename for the namelist

        real, intent(in) :: psfc ! Pressure at the surface [Pa]

        ! Output
        real, dimension(gr%nnzp), intent(out) :: 
     .  thlm,           ! Theta l mean                  [K] 
     .  rtm,            ! Total water mixing ratio      [kg/kg]
     .  um,             ! u wind                        [m/s]
     .  vm,             ! v wind                        [m/s]
     .  um_ref,         ! Initial profile of u wind     [m/s]
     .  vm_ref,         ! Initial profile of v wind     [m/s]
     .  ug,             ! u geostrophic wind            [m/s] 
     .  vg,             ! u geostrophic wind            [m/s] 
     .  wp2,            ! w'^2                          [m^2/s^2]
     .  rcm,            ! Cloud water mixing ratio      [kg/kg]
     .  wmt, wmm,       ! w wind                        [m/s]
     .  exner,          ! Exner function                [-] 
     .  em,             ! Turbulence kinetic energy     [m^2/s^2]
     .  p,              ! Pressure                      [Pa]
     .  rhot, rhom,     ! Density                       [kg/m^3]
     .  Lscale,         ! Mixing length                 [m] 
     .  Kht, Khm,       ! Eddy diffusivity              [m^2/s]
     .  taum, taut,     ! Dissipation time              [s]
     .  thvm            ! Virtual potential temperature [K]

        ! Input/Output
        real, dimension(gr%nnzp,hydromet_dim), intent(inout) :: 
     .  hydromet        ! Hydrometeor types             [units vary]

        ! Output
        real, dimension(gr%nnzp,sclr_dim), intent(out) :: 
     .  sclrm,  ! Standard passive scalar [units vary]
     .  edsclrm ! Eddy diffusivity passive scalar [units vary]

        ! Local Variables
        real, dimension(gr%nnzp) :: tmp1

        real :: cloud_top_height ! [m]
        real :: emax
        integer :: k, err_code

!-----------------------------------------------------------------------

        err_code = clubb_no_error
        
        ! Read sounding information

        call read_sounding( iunit, thlm, rtm, um, vm, ug, vg, 
     .                      runfile, runtype,
     .                      sclrm, edsclrm )


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
     .         "clex9_nov02", "clex9_oct14",
     .         "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .         "dycoms2_rf02_nd", "dycoms2_rf02_so" )
          ! thlm profile that is initially saturated at points.
          ! thlm profile remains the same as in the input sounding.
          ! use iterative method to find initial rcm.
          do k = 1, gr%nnzp, 1
             rcm(k) = sat_rcm( thlm(k), rtm(k), p(k), exner(k) )
          end do

        case default
          ! Initial profile is non-saturated thlm or any type of theta.
          thlm = thlm - Lv/(Cp*exner) * rcm

          ! Testing of passive scalars
          if ( iisclr_thl > 0 ) then
            sclrm(:,iisclr_thl) = sclrm(:,iisclr_thl) 
     .                          - Lv/(Cp*exner)*rcm
          end if

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
        case ( "generic" )
          em = 1.0
          wp2 = 2.0 / 3.0 * em

        ! GCSS BOMEX
        case ( "bomex" )
!---> Reduction of initial sounding for stability
!         do k = 1, gr%nnzp
!            em(k) = 1.0 - (gr%zm(k)/3000.0)
!            if ( em(k) < emin ) then
!               em(k) = emin
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nnzp) = em(gr%nnzp-1)

!         wp2 = 2.0 / 3.0 * em
!<--- End reduction of initial sounding for stability 24 Jan 07
          em(:) = emin
          wp2 = 2.0 / 3.0 * em

        ! GCSS ARM
        case ( "arm" )
!---> Reduction of initial sounding for stability
!         do k = 1, gr%nnzp
!            if ( gr%zm(k) < 150.0 ) then
!               em(k) = ( 0.15 * (1.0 - gr%zm(k)/150.0) ) / rhom(k)
!            else
!               em(k) = emin
!            end if
!         end do
!         em(1) = em(2)
!         em(gr%nnzp) = em(gr%nnzp-1)

!<--- End reduction of initial sounding for stability 24 Jan 07
          em(:) = emin
          wp2 = 2.0 / 3.0 * em

        ! March 2000 ARM case
        case ( "arm_0003" )
          em = 1.0
          wp2 = 2.0 / 3.0 * em
          call arm_0003_init()

        ! 3 year ARM case
        case ( "arm_3year" )
          em = 1.0
          wp2 = 2.0 / 3.0 * em
          call arm_3year_init()
         
        ! June 27 1997 ARM case
        case ( "arm_97" )
          em = 1.0
          wp2 = 2.0 / 3.0 * em
          call arm_97_init()

        ! GCSS FIRE Sc
        case ( "fire" ) 
          cloud_top_height = 700. ! 700 m is the top of the cloud in FIRE
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
               em(k) = 1.
            else
               em(k) = emin
            end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! GCSS ATEX
        case ( "atex" )
          um = max( um, -8. )

!---> Reduction of initial sounding for stability
!         do k = 1, gr%nnzp
!           em(k) = 1.0 - (gr%zm(k)/3000.0)
!           if ( em(k) < emin ) then
!             em(k) = emin
!           end if
!         end do
!         em(1) = em(2)
!         em(gr%nnzp) = em(gr%nnzp-1)

!         wp2 = 2.0 / 3.0 * em

!<--- End reduction of initial sounding for stability 24 Jan 07
          em(:) = emin
          wp2 = 2.0 / 3.0 * em

        ! GCSS DYCOMS II RF01
        case ( "dycoms2_rf01" ) 
          cloud_top_height = 800. ! 800 m is the top of the cloud in RF01
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
              em(k) = 0.5
            else
              em(k) = emin
            end if
          end do
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
          ! 4150 + 2800 m is the top of the cloud in Nov11
          cloud_top_height = 2800. + gr%zm(1)
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then

        ! Modification by Adam Smith, 08 April 2008
        ! Reducing the value of em appears to reduce error in the
        ! updated Nov.11 case
              em(k) = 0.01
        ! End of ajsmith4's modification

            else
              em(k) = emin
            end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)
        ! End Vince Larson's change.

          wp2 = 2.0 / 3.0  * em

        ! Adam Smith addition for June 25 altocumulus case.
        case ( "jun25_altocu" )

        ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
!          do k=1,gr%nnzp
!            if ( gr%zm(k) < 1400. ) then
!               em(k) = 0.1
!            else
!               em(k) = emin
!            end if
!          end do

        ! Note: emin = 1.0e-6, defined in constants.F
        ! Adam Smith, 28 June 2006
          do k = 1, gr%nnzp
            em(k) = 0.01
          end do


          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)
        ! End Vince Larson's change.

          wp2 = 2.0 / 3.0  * em

        ! End of ajsmith4's addition

        ! Adam Smith addition for CLEX-9: Nov. 02 altocumulus case.
        case ( "clex9_nov02" )

        ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
          ! 4150 + 1400 m is the top of the cloud in Nov11
          cloud_top_height = 2200. + gr%zm(1)
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
              em(k) = 0.01
            else
              em(k) = emin
            end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)
        ! End Vince Larson's change.

          wp2 = 2.0 / 3.0  * em

        ! End of ajsmith4's addition

        ! Adam Smith addition for CLEX-9: Oct. 14 altocumulus case.
        case ( "clex9_oct14" )

        ! Vince Larson reduced initial forcing.  4 Nov 2005
!          em = 1.0
!          em = 0.1
          ! 4150 + 1400 m is the top of the cloud in Nov11
          cloud_top_height = 3500. + gr%zm(1)
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
              em(k) = 0.01
            else
              em(k) = emin
            end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)
        ! End Vince Larson's change.

          wp2 = 2.0 / 3.0  * em

        ! End of ajsmith4's addition

        case ( "lba" )
          em = 0.1
          wp2 = 2./3. * em
          call lba_init()
        ! Michael Falk for mpace_a Arctic Stratus case.
        case ( "mpace_a" )

          cloud_top_height = 2000.
          emax = 1.0

          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
               em(k) = emax
            else
               em(k) = emin
            end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em
          call mpace_a_init
        ! Michael Falk for mpace_b Arctic Stratus case.
        case ( "mpace_b" )

          cloud_top_height = 1300. ! 1300 m is the cloud top in mpace_b.  Michael Falk 17 Aug 2006
          emax = 1.0

          do k=1,gr%nnzp

            if ( gr%zm(k) < cloud_top_height ) then
               em(k) = emax
            else
               em(k) = emin
            end if
          enddo
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

        ! Brian Griffin for COBRA CO2 case.
        case ( "cobra" )
          em = 0.1

          wp2 = 2.0 / 3.0 * em


        ! Michael Falk for RICO tropical cumulus case, 13 Dec 2006
        case ( "rico" )

          cloud_top_height = 1500.
          emax = 1.0
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
              em(k) = emax
            else
              em(k) = emin
            end if
          enddo

          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2(1:gr%nnzp) = 2.0 / 3.0 * em(1:gr%nnzp) ! Michael Falk reworded this.

        ! Michael Falk for GABLS2 case, 29 Dec 2006
        case ( "gabls2" )

          cloud_top_height = 800.  ! per GABLS2 specifications
          emax = 0.5
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
              em(k) = emax * (1 - (gr%zm(k)/cloud_top_height))
            else
              em(k) = emin
            end if
          end do

          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2(1:gr%nnzp) = 2.0 / 3.0 * em(1:gr%nnzp)

        end select

        ! End Initialize TKE and other fields as needed

        ! Compute mixing length

        call compute_length( thvm, thlm, rtm, rcm, 
     .               em, p, exner, Lscale, err_code )

        ! Dissipation time
        tmp1 = sqrt( max( wtol**2, zm2zt( em ) ) )
        taut = min( Lscale / tmp1, taumax )
        taum = min( ( max( zt2zm( Lscale ), 0.0 )
     .               / sqrt( max( wtol**2, em ) ) ), taumax )
!        taum = zt2zm( taut )

        ! Modification to damp noise in stable region
        do k=1,gr%nnzp
          if ( wp2(k) <= 0.005 ) then
            taut(k) = taumin
            taum(k) = taumin
          end if
        end do

        ! Eddy diffusivity coefficient
        ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)

        Kht = c_K * Lscale * tmp1
        Khm = c_K * max( zt2zm( Lscale ), 0.0 ) 
     .            * sqrt( max( em, emin ) )
!        Khm = zt2zm( Kht )

       ! Moved this to be more general -dschanen July 16 2007
       if ( luv_nudge ) then
         um_ref = um ! Michael Falk addition for nudging code.  27 Sep/1 Nov 2006
         vm_ref = vm ! ditto
       end if

       return
       end subroutine hoc_initialize
!-----------------------------------------------------------------------
       subroutine hoc_restart
     .             ( iunit, runfile, restart_path_case, 
     .               filename, time_restart, 
     .               thlm, rtm, um, vm,
     .               ug, vg, upwp, vpwp, wmt, wmm, 
     .               um_ref, vm_ref, wpthlp, wprtp,
     .               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .               sclrm, edsclrm
     .             )
!       Description:
!       Execute the necessary steps for the initialization of the 
!       HOC model to a designated point in the submitted GrADS file. 
!-----------------------------------------------------------------------
        use inputfields,only: 
     .      datafile, input_type, input_um, input_vm, ! Variable(s)
     .      input_rtm, input_thlm, input_wp2, input_ug,
     .      input_vg, input_rcm, input_wmt, input_exner,
     .      input_em, input_p, input_rhot, input_rhom,
     .      input_Lscale, input_Lup, input_Ldown, input_Kht,
     .      input_Khm, input_taum, input_taut, input_thvm, 
     .      input_rrm, input_rsnowm, input_ricem, 
     .      input_rgraupelm, input_wprtp, input_wpthlp,
     .      input_wp3, input_rtp2, input_thlp2, 
     .      input_rtpthlp, input_upwp, input_vpwp, 
     .      input_thlm_forcing, input_rtm_forcing,  
     .      input_up2, input_vp2, input_Scm, input_Ncm,
     .      input_Ncnm, input_Nim, input_cf, input_Nrm,
     .      input_Sct

        use inputfields, only: compute_timestep, grads_fields_reader ! Procedure(s)

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm ! Procedure(s)

        use constants, only: fstderr ! Variables(s)

        use parameters, only: sclr_dim ! Variables(s)

        use sounding, only: read_sounding ! Procedure(s)

        use stats_precision, only: time_precision ! Variable(s)

        use model_flags, only: luv_nudge ! Variable(s)

        use arm_0003, only: arm_0003_init ! Procedure(s)

        use arm_97, only: arm_97_init ! Procedure(s)

        use arm_3year, only: arm_3year_init ! Procedure(s)

        use mpace_a, only: mpace_a_init ! Procedure(s)

        use lba, only: lba_init ! Procedure(s)


        implicit none

        ! Input Variables
        integer, intent(in) :: iunit
        character(len=*), intent(in) :: 
     .  runfile,          ! Filename for the namelist
     .  filename,         ! GrADS file name
     .  restart_path_case ! Path to GrADS data for restart

        real(kind=time_precision), intent(in) ::
     .  time_restart

        ! Output Variables
        real, dimension(gr%nnzp), intent(inout) :: 
     .  thlm,           ! Theta l mean                 [K] 
     .  rtm,            ! Total water mixing ratio     [kg/kg]
     .  um,             ! u wind                       [m/s]
     .  vm,             ! v wind                       [m/s]
     .  um_ref,         ! Initial profile of u wind    [m/s]
     .  vm_ref,         ! Initial profile of v wind    [m/s]
     .  ug,             ! u geostrophic wind           [m/s] 
     .  vg,             ! v geostrophic wind           [m/s] 
     .  wmt, wmm,       ! w wind                       [m/s]
     .  wprtp,          ! w' r_t'                      [(kg m)(kg s)]
     .  wpthlp,         ! w' th_l'                     [(m K)/s]
     .  upwp,           ! u'w'                         [m^2/s^2]
     .  vpwp            ! v'w'                         [m^2/s^2]

        real, intent(out) ::
     .  wpthlp_sfc,     ! w'theta_l' surface flux   [(m K)/s]
     .  wprtp_sfc,      ! w'rt' surface flux        [(m kg)/(kg s)]
     .  upwp_sfc,       ! u'w' at surface           [m^2/s^2] 
     .  vpwp_sfc        ! v'w' at surface           [m^2/s^2]
        
        ! Input/Output
        real, dimension(gr%nnzp,sclr_dim), intent(inout) :: 
     .  sclrm,  ! Standard passive scalar [units vary]
     .  edsclrm ! Eddy diffusivity passive scalar [units vary]
       
        integer :: timestep

        ! Inform inputfields module
        datafile = "../"//trim( restart_path_case )
        input_type = "hoc"
        input_um   = .true.
        input_vm   = .true.
        input_rtm  = .true.
        input_thlm = .true.
        input_wp2  = .true.
        input_ug   = .true.
        input_vg   = .true.
        input_rcm  = .true.
        input_wmt  = .true.
        input_exner = .true.
        input_em = .true.
        input_p = .true.
        input_rhot = .true.
        input_rhom = .true.
        input_Lscale = .true.
        input_Lup = .true.
        input_Ldown = .true.
        input_Kht = .true.
        input_Khm = .true.
        input_taum = .true.
        input_taut = .true.
        input_thvm = .true.
        input_rrm = .true.
        input_rsnowm = .true.
        input_ricem = .true.
        input_rgraupelm = .true.
        input_wprtp = .true.
        input_wpthlp = .true.
        input_wp3 = .true.
        input_rtp2 = .true.
        input_thlp2 = .true.
        input_rtpthlp = .true.
        input_upwp = .true.
        input_vpwp = .true. 
        input_thlm_forcing = .true.
        input_rtm_forcing = .true. 
        input_up2 = .true.
        input_vp2 = .true.
        input_Scm = .true.
        input_Ncm = .true.
        input_Ncnm = .true.
        input_Nim = .true.
        input_cf  = .true.
        input_Nrm = .true.
        input_Sct = .true.
        
        ! Determine the nearest timestep in the GRADS file to the
        ! restart time.
        call compute_timestep
     .       ( iunit, "../"//trim( restart_path_case )//"_zt.ctl", 
     .         .true., time_restart, timestep )
        
        ! Sanity check for input time_restart
        if ( timestep < 0 ) then
          write(fstderr,*) "Invalid time_restart in "//
     .      "file: "//trim( runfile )
          stop
        end if

        ! Read in sounding to get appropriate nudging information for um
        ! and vm
      
        call read_sounding( iunit, thlm, rtm, um, vm, ug, vg, 
     .                      runfile, runtype,
     .                      sclrm, edsclrm )

        if ( luv_nudge ) then
          um_ref = um
          vm_ref = vm
        end if
        
        ! Read data from GrADS files
        call grads_fields_reader( timestep )
        
        ! Initialize forcing files for specific cases
        select case( trim( runtype ) )
        case( "arm_3year" )
            call arm_3year_init()
        case( "arm_97" )
            call arm_97_init()
        case( "arm_0003" )
            call arm_0003_init()
        case( "lba" )
            call lba_init()
        case( "mpace_a" )
            call mpace_a_init()
        end select
        
        wmm = zt2zm( wmt )
        
        wpthlp_sfc = wpthlp(1)
        wprtp_sfc  = wprtp(1)
        upwp_sfc   = upwp(1)
        vpwp_sfc   = vpwp(1)

        return
        end subroutine hoc_restart

!----------------------------------------------------------------------
        subroutine hoc_forcings_timestep
     .               ( dt, iteration, stdout, err_code )

!       Description:
!       Calculate tendency and surface variables

!       Calls: (* = model case)
!         subroutines *_sfclyr
!         subroutines *_tndncy
!         subroutine sfc_thermo_fluxes
!         subroutine sfc_momentum_fluxes
!----------------------------------------------------------------------

        ! Modules to be included
        use model_flags, only: kk_rain, licedfs, lcoamps_micro, ! Variable(s)
     .                   cloud_sed, lbugsrad

        use constants, only: rc_tol, fstderr ! Variable(s)

        use grid_class, only: gr ! Variable(s)

        use grid_class, only: zt2zm ! Procedure(s)

        use diagnostic_variables, only: hydromet, Ncm, radht, um_ref, ! Variable(s)
     .                            vm_ref, Frad, Ncnm, thvm, ustar,
     .                            pdf_parms, Khm, Akm_est, Akm, Nim

        use diagnostic_variables, only: wpedsclrp ! Passive scalar variables
        
        use prognostic_variables, only: rtm_forcing, thlm_forcing, ! Variable(s)
     .                            wmt, wmm, rhot, rtm, thlm, p,
     .                            exner, rcm, rhom, um, psfc, vm,
     .                            upwp_sfc, vpwp_sfc, Tsfc,
     .                            wpthlp_sfc, SE, LE, wprtp_sfc, cf


        use prognostic_variables, only: 
     .      sclrm, sclrm_forcing,  ! Passive scalar variables
     .      wpsclrp,
     .      wpsclrp_sfc, wpedsclrp_sfc
        
        use stats_precision, only: time_precision ! Variable(s)

        use numerical_check, only: isnan2d, rad_check ! Procedure(s)

        use microphys_driver, only: advance_microphys ! Procedure(s)

        use error_code, only: lapack_error, ! Procedure(s)
     .                        clubb_at_debug_level

        use array_index, only:
     .      iirsnowm, iiricem,
     .      iisclr_rt, iisclr_thl
       
        ! Case specific modules 
        use arm, only: arm_tndcy, arm_sfclyr ! Procedure(s)

        use arm_97, only: arm_97_tndcy, arm_97_sfclyr ! Procedure(s)

        use arm_0003, only: arm_0003_tndcy, arm_0003_sfclyr ! Procedure(s)

        use arm_3year, only: arm_3year_tndcy, arm_3year_sfclyr ! Procedure(s)

        use astex, only: astex_tndcy, astex_sfclyr ! Procedure(s)

        use atex, only: atex_tndcy, atex_sfclyr ! Procedure(s)

        use bomex, only: bomex_tndcy, bomex_sfclyr ! Procedure(s)

        use cobra, only: cobra_tndcy, cobra_sfclyr ! Procedure(s)

        use dycoms2_rf01, only: 
     .      dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr ! Procedure(s)

        use dycoms2_rf02, only: 
     .      dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr ! Procedure(s)

        use cloud_sed_mod, only: cloud_drop_sed ! Procedure(s)

        use fire, only: fire_tndcy, sfc_momentum_fluxes,
     .                  sfc_thermo_fluxes ! Procedure(s)

        use gabls2, only: gabls2_tndcy, gabls2_sfclyr ! Procedure(s)

        use rico, only: rico_tndcy, rico_sfclyr ! Procedure(s)

        use lba, only: lba_tndcy, lba_sfclyr ! Procedure(s)

        use mpace_a, only: mpace_a_tndcy, mpace_a_sfclyr ! Procedure(s)

        use mpace_b, only: mpace_b_tndcy, mpace_b_sfclyr ! Procedure(s)

        use nov11, only: nov11_altocu_tndcy ! Procedure(s)

        use jun25, only: jun25_altocu_tndcy ! Procedure(s)

        use clex9_nov02, only: clex9_nov02_tndcy ! Procedure(s)

        use clex9_oct14, only: clex9_oct14_tndcy ! Procedure(s)

        use wangara, only: wangara_tndcy, wangara_sfclyr ! Procedure(s)

#ifdef radoffline
        use bugsrad_hoc_mod, only: bugsrad_hoc ! Procedure(s)
#endif


        implicit none

        ! Input Variables
        real(kind=time_precision), intent(in) ::
     .    dt    ! Timestep      [s]

        integer, intent(in) ::
     .    iteration ! The iteration number

        logical, intent(in) ::
     .    stdout ! Whether to print messages to stdout

        ! Input/Output Variables
        integer, intent(inout) ::
     .    err_code

        ! Local Variables

        real, dimension(gr%nnzp) :: 
     .    rsnowm, ! Rain water mixing ratio              [kg/kg]
     .    ricem   ! Prisitine ice water mixing ratio     [kg/kg]

        integer :: lin_int_buffer

        integer :: k ! Vertical loop index variable


!-----------------------------------------------------------------------

!#######################################################################
!##############      FIND ALL DIAGNOSTIC VARIABLES        ##############
!#######################################################################

         select case ( runtype )
         case( "generic" ) ! Generic run
           ! Configure for K&K microphysics
           do k=1, gr%nnzp, 1
             if ( rcm(k) >= rc_tol ) then
              ! Ncm is in units of kg^-1.  If the coefficient is in m^-3, then
              ! it needs to be divided by rhot in order to get units of kg^-1.
              ! Brian.  Sept. 8, 2007.
!               Ncm(k) = 30.0 * (1.0 + exp( -gr%zt(k)/2000.0 )) * 1.e6
!     .                  * rhot(k) 
               Ncm(k) = 30.0 * (1.0 + exp( -gr%zt(k)/2000.0 )) * 1.e6
     .                  / rhot(k) 
             end if
           end do

         case( "arm" ) ! ARM Cu case
           call arm_tndcy
     .          ( time_current, thlm_forcing, radht, rtm_forcing,
     .            sclrm_forcing )

         case( "arm_0003" ) ! ARM March 2000 case
          call arm_0003_tndcy( time_current, wmm, wmt, thlm_forcing, 
     .                         rtm_forcing, um_ref, vm_ref,
     .                         sclrm_forcing )

         case( "arm_3year" ) ! ARM 3 year case
          call arm_3year_tndcy( time_current, wmm, wmt, thlm_forcing, 
     .                          rtm_forcing, um_ref, vm_ref,
     .                          sclrm_forcing )

         case( "arm_97" ) ! 27 June 1997 ARM case
           call arm_97_tndcy( time_current, wmm, wmt, thlm_forcing, 
     .                        rtm_forcing, um_ref, vm_ref,
     .                        sclrm_forcing )

         case( "bomex" ) ! BOMEX Cu case
           call bomex_tndcy( time_current, wmt, wmm, radht,
     .                       thlm_forcing, rtm_forcing,
     .                       sclrm_forcing )

         case( "fire" ) ! FIRE Sc case
           call fire_tndcy( time_current, rhot, rcm, exner, 
     .                      wmt, wmm, Frad, radht,
     .                      thlm_forcing, rtm_forcing,
     .                      sclrm_forcing )

         case( "wangara" ) ! Wangara dry CBL
           call wangara_tndcy( time_current, wmt, wmm, 
     .                         thlm_forcing, rtm_forcing,
     .                         sclrm_forcing )

         case( "atex" ) ! ATEX case
           call atex_tndcy( time_current, time_initial, rtm, 
     .                      rhot, rcm, exner, wmt, wmm, Frad, radht,
     .                      thlm_forcing, rtm_forcing, err_code,
     .                      sclrm_forcing )

         case( "dycoms2_rf01" ) ! DYCOMS2 RF01 case
           call dycoms2_rf01_tndcy( time_current, rhot, rhom, rtm, rcm, 
     .                              exner, wmt, wmm, Frad, radht, 
     .                              thlm_forcing, rtm_forcing, err_code,
     .                              sclrm, sclrm_forcing )

         case( "astex_a209" ) ! ASTEX Sc case for K & K
           call astex_tndcy( time_current, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing,
     .                       sclrm_forcing )

         case( "dycoms2_rf02_do", ! DYCOMS2 RF02 case with drizzle only.
     .         "dycoms2_rf02_ds", ! DYCOMS2 RF02 case with drizzle and cloud sedimentation.
     .         "dycoms2_rf02_nd", ! DYCOMS2 RF02 case with no drizzle and no cloud sedimentation.
     .         "dycoms2_rf02_so" )! DYCOMS2 RF02 case with cloud water sedimentation only.
           call dycoms2_rf02_tndcy
     .          ( time_current, time_initial, rhot, rhom, rtm, rcm, 
     .            exner, wmt, wmm, thlm_forcing, rtm_forcing,
     .            Frad, radht, Ncm, Ncnm, err_code,
     .            sclrm_forcing  )

         case( "nov11_altocu" ) ! Nov. 11 Altocumulus case.
           call nov11_altocu_tndcy
     .          ( time_current, time_initial, dt, rlat, rlon,
     .            thlm, rcm, p, exner, rhot, rtm, wmt,
     .            wmm, thlm_forcing, rtm_forcing,
     .            Frad, radht, Ncnm,
     .            sclrm_forcing )

         case( "jun25_altocu" ) ! June 25 Altocumulus case.
           call jun25_altocu_tndcy
     .          ( time_current, time_initial, rlat, dt, thlm,
     .            rcm, p, exner, rhot, rtm, wmt,
     .            wmm, thlm_forcing, rtm_forcing,
     .            Frad, radht, sclrm_forcing )

         case( "clex9_nov02" ) ! CLEX-9: Nov. 02 Altocumulus case.
           call clex9_nov02_tndcy
     .          ( time_current, time_initial, dt, rlat, rlon,
     .            thlm, rcm, p, exner, rhot, rtm, wmt,
     .            wmm, thlm_forcing, rtm_forcing,
     .            Frad, radht, Ncnm, sclrm_forcing )

         case( "clex9_oct14" ) ! CLEX-9: Oct. 14 Altocumulus case.
           call clex9_oct14_tndcy
     .          ( time_current, time_initial, dt, rlat, rlon,
     .            thlm, rcm, p, exner, rhot, rtm, wmt,
     .            wmm, thlm_forcing, rtm_forcing,
     .            Frad, radht, Ncnm ,sclrm_forcing )

         case ( "lba" )
           call lba_tndcy( time_current, wmt, wmm, radht, 
     .                     thlm_forcing, rtm_forcing,
     .                     sclrm_forcing )

         case ( "mpace_a" ) ! mpace_a arctic stratus case
           call mpace_a_tndcy
     .          ( time_current, time_initial, dtmain, rlat, thlm,
     .            exner, rhot, rtm, p, thvm, rcm,
     .            wmt, wmm, thlm_forcing, rtm_forcing,
     .            Ncnm, Ncm, Frad, radht, um_ref, vm_ref,
     .            sclrm_forcing )

         case ( "mpace_b" ) ! mpace_b arctic stratus case
           call mpace_b_tndcy
     .          ( time_current, time_initial, dtmain, rlat, thlm,
     .            exner, rhot, rtm, p, thvm, rcm,
     .            wmt, wmm, thlm_forcing, rtm_forcing,
     .            Ncnm, Ncm, Frad, radht,
     .            sclrm_forcing )

        ! Brian Griffin for COBRA CO2 case.
        case ( "cobra" )
           call cobra_tndcy( time_current, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing,
     .                       sclrm_forcing )

         case ( "rico" ) ! RICO case
           call rico_tndcy( time_current, time_initial, dtmain, exner,
     .                      rhot, rcm, kk_rain, wmt, wmm,
     .                      thlm_forcing, rtm_forcing, radht, Ncm,
     .                      sclrm_forcing )

         case ( "gabls2" ) ! GABLS 2 case
           call gabls2_tndcy( time_current, time_initial, 
     .                        rhot, rcm, kk_rain, wmt, wmm,
     .                        thlm_forcing, rtm_forcing, radht, Ncm,
     .                        sclrm_forcing )

         case default
           write(unit=fstderr,fmt=*) 
     .       "hoc_forcings_timestep: Don't know how to handle "
     .       //"LS forcing for runtype: "//trim( runtype )
           stop

         end select

         ! Bc for the second order moments

        select case ( trim( runtype ) )

        case( "generic", "fire" )  ! Normal and FIRE
          call sfc_momentum_fluxes( um(2), vm(2),
     .                              upwp_sfc, vpwp_sfc, ustar )
          ! sfctype = 0  fixed sfc sensible and latent heat fluxes 
          !                   as given in hoc.in
          ! sfctype = 1  bulk formula: uses given surface temperature 
          !                   and assumes over ocean
          if ( sfctype == 0 ) then
            wpthlp_sfc = SE
            wprtp_sfc  = LE
            if ( iisclr_thl > 0 ) wpsclrp(:,iisclr_thl) = SE
            if ( iisclr_rt > 0 ) wpsclrp(:,iisclr_rt)   = LE
            if ( iisclr_thl > 0 ) wpedsclrp(:,iisclr_thl) = SE
            if ( iisclr_rt > 0 ) wpedsclrp(:,iisclr_rt)   = LE
          else if ( sfctype == 1 ) then
            call sfc_thermo_fluxes( um(2), vm(2),
     .                              Tsfc, psfc, 
     .                              thlm(2), rtm(2),
     .                              wpthlp_sfc, wprtp_sfc,
     .                              sclrm(2,:), wpsclrp(1,:)
     .                            )

          else
            write(unit=fstderr,fmt=*) 
     .        "Invalid value of sfctype = ", sfctype
            stop

          end if

        case( "arm" )
          call arm_sfclyr( time_current, gr%zt(2), 1.1, 
     .                     thlm(2), um(2), vm(2),
     .                     upwp_sfc, vpwp_sfc, 
     .                     wpthlp_sfc, wprtp_sfc, ustar,
     .                     wpsclrp_sfc, wpedsclrp_sfc )

        case( "arm_0003" )
           call arm_0003_sfclyr( time_current, gr%zt(2), rhom(1),
     .                         thlm(2), um(2), vm(2),
     .                         upwp_sfc, vpwp_sfc, 
     .                         wpthlp_sfc, wprtp_sfc, ustar,
     .                         wpsclrp_sfc, wpedsclrp_sfc )

        case( "arm_3year" )
           call arm_3year_sfclyr( time_current, gr%zt(2), rhom(1),
     .                            thlm(2), um(2), vm(2),
     .                            upwp_sfc, vpwp_sfc, 
     .                            wpthlp_sfc, wprtp_sfc, ustar,
     .                            wpsclrp_sfc, wpedsclrp_sfc )


         case ( "arm_97" )
           call arm_97_sfclyr( time_current, gr%zt(2), rhom(1),
     .                         thlm(2), um(2), vm(2),
     .                         upwp_sfc, vpwp_sfc, 
     .                         wpthlp_sfc, wprtp_sfc, ustar,
     .                         wpsclrp_sfc, wpedsclrp_sfc ) 
        case( "bomex" ) 
          call bomex_sfclyr( um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc,
     .                       wpthlp_sfc, wprtp_sfc, ustar,
     .                       wpsclrp_sfc, wpedsclrp_sfc )

        case( "wangara" )
          call wangara_sfclyr( time_current, um(2), vm(2),
     .                         upwp_sfc, vpwp_sfc,
     .                         wpthlp_sfc, wprtp_sfc, ustar,
     .                         wpsclrp_sfc, wpedsclrp_sfc )

        case( "atex" )
          call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                      upwp_sfc, vpwp_sfc,
     .                      wpthlp_sfc, wprtp_sfc, ustar,
     .                      sclrm(2,:), wpsclrp_sfc, wpedsclrp_sfc )

        case( "dycoms2_rf01" ) 
          call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, 
     .                              exner(1), um(2), vm(2), 
     .                              thlm(2), rtm(2), 
     .                              rhom(1), upwp_sfc, vpwp_sfc, 
     .                              wpthlp_sfc, wprtp_sfc, ustar,
     .                              sclrm(2,:), wpsclrp_sfc, 
     .                              wpedsclrp_sfc )
        case( "astex_a209" )
          call astex_sfclyr( rhom(1), um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc, wpthlp_sfc, 
     .                       wprtp_sfc ,wpsclrp_sfc, wpedsclrp_sfc )

        case( "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .        "dycoms2_rf02_nd", "dycoms2_rf02_so" )
          call dycoms2_rf02_sfclyr( um(2), vm(2),
     .                              upwp_sfc, vpwp_sfc,
     .                              wpthlp_sfc, wprtp_sfc, ustar,
     .                              wpsclrp_sfc, wpedsclrp_sfc )

        case( "nov11_altocu" )
          ! There are no surface momentum or heat fluxes
          ! for the Nov. 11 Altocumulus case.

        case( "jun25_altocu" )
          ! There are no surface momentum or heat fluxes
          ! for the Jun. 25 Altocumulus case.

        case( "clex9_nov02" )
          ! There are no surface momentum or heat fluxes
          ! for the CLEX-9: Nov. 02 Altocumulus case.
       
        case( "clex9_oct14" )
          ! There are no surface momentum or heat fluxes
          ! for the CLEX-9: Oct. 14 Altocumulus case.

         case ( "lba" )
           call lba_sfclyr( time_current, gr%zt(2), rhom(1),
     .                      thlm(2), um(2), vm(2),
     .                      upwp_sfc, vpwp_sfc, 
     .                      wpthlp_sfc, wprtp_sfc, ustar,
     .                      wpsclrp_sfc, wpedsclrp_sfc )

         case ( "mpace_a" )
           call mpace_a_sfclyr
     .          ( time_current, rhom(1), um(2), vm(2), upwp_sfc,
     .            vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar,
     .            wpsclrp_sfc, wpedsclrp_sfc )

         case ( "mpace_b" )
           call mpace_b_sfclyr( rhom(1), um(2), vm(2), upwp_sfc,
     .                          vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar,
     .                          wpsclrp_sfc, wpedsclrp_sfc )

        ! Brian Griffin for COBRA CO2 case.
        case ( "cobra" )
          call cobra_sfclyr( time_current, gr%zt(2), rhom(1),
     .                       thlm(2), um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc,
     .                       wpthlp_sfc, wprtp_sfc, ustar,
     .                       wpsclrp_sfc, wpedsclrp_sfc )

         case ( "rico" )
           call rico_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                        gr%zt(2), 299.8, 101540., ! 299.8 K is the RICO SST; 101540 Pa is the sfc pressure.
!     .                        gr%zt(2), Tsfc, psfc,
     .                        upwp_sfc, vpwp_sfc, wpthlp_sfc,
     .                        wprtp_sfc, ustar,
     .                        sclrm(2,:), wpsclrp_sfc, wpedsclrp_sfc )

         case ( "gabls2" )
          call gabls2_sfclyr
     .         ( time_current, time_initial, gr%zt(2), 97200.,
     .           um(2), vm(2), thlm(2), rtm(2),
     .           upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, ustar,
     .           sclrm(2,:), wpsclrp_sfc, wpedsclrp_sfc )

        case default
          write(unit=fstderr,fmt=*) 
     .      "Invalid value of runtype = ", runtype
          stop

        end select ! runtype

        !----------------------------------------------------------------
        ! Compute Microphysics
        !----------------------------------------------------------------
        ! Call Khairoutdinov and Kogan (2000) scheme, or COAMPS
        ! for rain microphysics.
        
        if ( kk_rain .or. lcoamps_micro .or. licedfs ) then
          call advance_microphys
     .         ( runtype, dt, time_current, time_initial,
     .           thlm, p, exner, rhot, rhom, rtm, rcm, Ncm, 
     .           pdf_parms, wmt, wmm, Khm, AKm_est, Akm, 
     .           Ncnm, Nim,
     .           hydromet,
     .           rtm_forcing, thlm_forcing, err_code )

          if ( lapack_error(err_code) ) return

        end if

        if ( cloud_sed ) then
          call cloud_drop_sed( rcm, Ncm, rhom, rhot, exner,
     .                         rtm_forcing, thlm_forcing ) 
        end if


        if ( lbugsrad ) then
#ifdef radoffline /*This directive is needed for BUGSrad to work with HOC.*/

          ! Assign pointers to snow and ice
          if ( iirsnowm > 0 ) then
            rsnowm = hydromet(1:gr%nnzp,iirsnowm)
          else
            rsnowm = 0.0
          end if

          if ( iiricem > 0 ) then
            ricem = hydromet(1:gr%nnzp,iiricem)
          else
            ricem = 0.0
          end if

          ! NaN checks added to detect possible errors
          ! with BUGSrad
          ! Joshua Fasching November 2007

          if ( isnan2d( thlm ) ) then
            print *, "thlm before BUGSrad is NaN" 
          endif

          if ( isnan2d( rcm ) ) then
            print *, "rcm before BUGSrad is NaN" 
          endif

          if ( isnan2d( rtm ) ) then
            print *, "rtm before BUGSrad is NaN" 
          endif

          if ( isnan2d( rsnowm ) ) then
            print *, "rsnowm before BUGSrad is NaN" 
          endif

          if ( isnan2d( ricem ) ) then
            print *, "ricem before BUGSrad is NaN" 
          endif

          if ( isnan2d( cf ) ) then
            print *, "cf before BUGSrad is NaN" 
          end if

          if ( isnan2d( p ) ) then
            print *, "p before BUGSrad is NaN" 
          end if

          if ( isnan2d( exner ) ) then
            print *, "exner before BUGSrad is NaN" 
          end if

          if ( isnan2d( rhom ) ) then
            print *, "rhom before BUGSrad is NaN" 
          end if

          if ( isnan2d( thlm_forcing ) ) then
            print *, "thlm_forcing before BUGSrad is NaN" 
          end if

          ! Check for impossible negative values
          if ( clubb_at_debug_level( 2 ) ) then
            call rad_check( thlm, rcm, rtm, ricem,
     .                      cf, p, exner, rhom )
          end if

          ! Initially we will set this to a constant for testing purposes
          ! lin_int_buffer = 20

          ! Use a a new formula that creates and evenly spaced grid
          ! between the model domain top and the standard atmosphere
          ! table.  e.g. if the HOC model top is 3200m, and the spacing
          ! between gr%nnzp-1 and gr%nnzp is 40m, then lin_int_buffer is
          ! 19 and each layer of the buffer is 40m deep. -dschanen 14 May 08
          lin_int_buffer = 
     .    max( int( ( 1000.-mod( gr%zm(gr%nnzp), 1000. ) )
     .           * gr%dzm(gr%nnzp) ) - 1, 0 )

          ! print *, "buffer = ", lin_int_buffer !%% debug

          call bugsrad_hoc( gr%zm, gr%nnzp, lin_int_buffer, ! In
     .                      rlat, rlon,                     ! In
     .                      day, month, year, time_current, ! In
     .                      thlm, rcm, rtm, rsnowm, ricem,  ! In
     .                      cf, p, zt2zm( p ), exner, rhom, ! In
     .                      radht, Frad,                    ! Out
     .                      thlm_forcing )                  ! In/Out
          
          if ( isnan2d( thlm_forcing ) ) then
            print *, "thlm_forcing after BUGSrad is NaN" 
!           print *,thlm_forcing

          end if

          if ( isnan2d( Frad ) ) then
            print *, "Frad after BUGSrad is NaN" 
!	    print *, Frad
          end if

          if ( isnan2d( radht ) ) then
            print *, "radht after BUGSrad is NaN"
!           print *,radht                                             
          end if

#else
          stop "Cannot call BUGSrad with these compile options."
#endif /*radoffline*/

        end if ! lbugsrad

        return
        end subroutine hoc_forcings_timestep

        end module hoc
