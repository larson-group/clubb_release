!-----------------------------------------------------------------------
! $Id: hoc.F,v 1.126 2007-07-06 23:01:39 dschanen Exp $

#ifdef SCALARS /* Test scalars for theta and rt */
#define SCLR_THETA 1
#define SCLR_RT 2
#endif /*SCALARS*/

        module hoc

!       Description:
!       Contains the necessary subroutines to execute individual HOC 
!       model runs, using one of the driver programs (the simplest case
!       being the hoc_standalone program).
!-----------------------------------------------------------------------

        implicit none

        ! Setup hoc_model() as the sole external interface
        private :: hoc_initialize
        private :: hoc_forcings_timestep

        private :: hoc_model_standard, hoc_model_field_input

        public  :: hoc_model

        ! Overloaded interface for inputfields;  There may be a better
        ! way to do this, but I don't know what it would be.
        ! -dschanen 21 August 2006
        interface hoc_model
          module procedure hoc_model_field_input, hoc_model_standard
        end interface

        private ! Default to private

        ! Model settings

        ! Grid definition
        integer, private :: 
     .  nzmax  ! Vertical extent in levels              [#]

        real, private :: 
     .  deltaz, ! Change per grid level                 [m]
     .  zm_init ! Initial point on the momentum grid    [m]

!$omp   threadprivate(nzmax, zm_init, deltaz)

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

        real, private :: 
     .  time_initial, ! Time of start of simulation     [s]
     .  time_final,   ! Time end of simulation          [s]
     .  time_spinup,  ! Time end of spin up period      [s]
     .  time          ! Current time of simulation      [s]
!$omp   threadprivate(time_initial, time_final, time_spinup, time)

        real, private :: 
     .  dtmain,     ! Main model timestep                      [s]
     .  dtclosure,  ! Closure model timestep                   [s]
     .  dt          ! Current model timestep (based on spinup) [s]
!$omp   threadprivate(dtmain, dtclosure, dt)

        contains
!-----------------------------------------------------------------------

        subroutine hoc_model_field_input
     .             ( params, runfile, isValid, stdout,
     .               set_fields, sample_ratio )
!       Description:
!       Subprogram to integrate the pde equations for pdf closure
!       This version is modified to set fields from pre-calculated data

!       Calls:  subroutine hoc_initialize (once)
!          subroutine hoc_forcings_timestep (ifinal times)
!          subroutine hoc_closure_timestep (ifinal*niterlong times)
!          subroutine deallocate_model_arrays (once)
!          function invalid_model_arrays
!          subroutine set_fields (passed as a parameter) 


!-----------------------------------------------------------------------

        use diagnostic_variables
        use prognostic_variables
        use grid_class
        use param_index, only: nparams
        use numerical_check
        use hoc_parameterization_interface
        use constants, only: fstdout, fstderr
#ifdef STATS
        use hoc_stats
#endif

        implicit none

        interface
          subroutine set_fields
     .               ( timestep,
     .                um, vm, rtm, thlm, 
     .                wp2, wprtp, wpthlp, wp3,
     .                rtp2, thlp2, rtpthlp, upwp, vpwp )
          implicit none

          integer, intent(in) :: 
     .    timestep

          real, optional, dimension(:), intent(out) :: 
     .    um,     ! u wind                              [m/s]
     .    vm,     ! v wind                              [m/s]
     .    rtm,    ! total water mixing ratio            [kg/kg]
     .    thlm,   ! liquid potential temperature        [K]
     .    wp2,    ! w'^2                                [m^2/s^2]
     .    wprtp,  ! w'r_t'                              [(m kg)/(s kg)]
     .    wpthlp, ! w'th_l'                             [(m K)/ s] 
     .    wp3,    ! w'^3                                [m^3/s^3]
     .    rtp2,   ! r_t'^2                              [kg/kg/s]
     .    thlp2,  ! th_l'^2                             [K^2/s]
     .    rtpthlp,! r_t' th_l'                          [(kg K)/(kg s)]
     .    upwp,   ! u'w'                                [m^2/s^2]        
     .    vpwp    ! v'w'                                [m^2/s^2]

          end subroutine set_fields
        end interface
     
        ! External
        intrinsic :: mod, real, int 

        ! Input Variables
        real, intent(in), dimension(nparams) :: 
     .  params  ! Model parameters, C1, nu2, etc.

        ! Subroutine Arguments (Model Setting)
        character(len=*), intent(in) :: 
     .  runfile ! Name of file containing &model_setting and &sounding
        logical, intent(in) :: 
     .  stdout  ! Whether to print output per timestep

        ! Sample_ratio is used to determine often to sample data from the GrADS 
        ! input per iteration of the main time stepping loop.
        ! E.g. if dtmain is 20.0 seconds and the data was written every 60.0 sec, 
        ! then sample ratio should be equal to 3.0.
        real, intent(in) :: sample_ratio 

        ! Output Variables
        logical, intent(inout) :: isValid    ! valid run? (T/F)

        ! Local Variables
        ! Internal Timing Variables
        integer ::
     .  ifinal,
     .  niterlong,
     .  sample_time

        real :: 
     .  fcor, ! Coriolis parameter       [s^-1]
     .  T0    ! Reference Temperature    [K]

        logical :: 
     .  cloud_sed,    ! Flag for cloud water droplet sedimentation. - Brian
     .  kk_rain,      ! Flag for Khairoutdinov and Kogan rain microphysics. - Brian
     .  licedfs,      ! Flag for simplified ice scheme
     .  lcoamps_micro,! Flag for COAMPS microphysical scheme
     .  lbugsrad,     ! Flag for BUGSrad interactive radiation scheme
     .  lmpace_nudge, ! Whether to adjust the winds within the timestep
     .  lKhm_aniso    ! Whether to use anisotropic Khm.  - Michael Falk 2 Feb 2007

        ! Dummy grid height arrays
        real, dimension(:), allocatable :: 
     .  dummy_array_m, dummy_array_t

       integer :: i, i1 ! Internal Loop Variables

       ! Definition of namelists
       namelist /model_setting/ 
     .   runtype, nzmax, deltaz, zm_init,
     .   day, month, year, rlat, rlon,
     .   time_initial, time_final, time_spinup,
     .   dtmain, dtclosure,
     .   sfctype, Tsfc, psfc, SE, LE, fcor, T0,
     .   cloud_sed, kk_rain, licedfs, lcoamps_micro, 
     .   lbugsrad, lKhm_aniso

       !  Initialize the model run 

       ! Initialize run to valid
       !isValid = .true.

       ! Read namelist file
       open(unit=10, file=runfile, status='old')
       read(unit=10, nml=model_setting)
       close(unit=10)

       ! Dummy grid height arrays
       allocate( dummy_array_m(1:nzmax), dummy_array_t(1:nzmax) )

       dummy_array_m = 0.0
       dummy_array_t = 0.0

       if ( trim( runtype ) == "mpace" ) then 
         lmpace_nudge = .true.
       else
         lmpace_nudge = .false.
       end if

       ! Allocate & initialize variables,
       ! setup grid, setup constants, and setup flags

       call parameterization_setup
     .      ( nzmax, deltaz, zm_init, T0, params,
     .        lbugsrad, kk_rain, licedfs, lcoamps_micro,
     .        cloud_sed,
     .        lmpace_nudge, lKhm_aniso, .false., 
     .        dummy_array_m, dummy_array_t, isValid )

       if ( .not. isValid ) return


       ! Dummy grid height arrays
       deallocate( dummy_array_m, dummy_array_t )

       time = time_initial

       call hoc_initialize( 10, runfile, psfc, thlm, rtm, 
     .                      um, vm, ug, vg, wp2,
     .                      rcm, hydromet,
     .                      wmt, wmm, em, exner,
     .                      taut, taum, thvm, p,
     .                      rhot, rhom, Lscale, lup,
     .                      ldown, Kht, Khm
#ifdef SCALARS
     .                      ,sclrm, edsclrm
#endif /*SCALARS*/
     .                     )

#ifdef STATS
!      Initialize statistics output
       call stats_init( 50, runfile, gr%nnzp, gr%zt, gr%zm,
     .                  day, month, year, rlat, rlon, time, dtmain )
#endif /*STATS*/

       ! Time integration
       ! Call hoc_closure_timestep once per each GrADS output time 
       ifinal  = ( time_final - time_initial ) / dtmain

       sample_time = 1

!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

       ! Run from iteration 1 to iteration ifinal + 1 in order to output
       ! information to the GrADS file for the final model output time.
       ! Brian.
       do i = 1, ifinal-1, 1

#ifdef STATS
         ! When this time step is over, the time will be time + dtmain

         ! We use elapsed time for stats_begin_step
         call stats_begin_timestep( time-time_initial+dtmain, dtmain )

         if ( lstats_last .and. stdout ) then
           write(unit=fstdout,fmt='(a,i7,a,f8.1)') 'iteration = ',
     .       i, '; time = ', time
         end if
#endif /*STATS*/

         ! Sample from the data file
         
         if ( mod( real( i ), sample_ratio) == 0.0 ) then
           sample_time = int( i / sample_ratio ) + 1
         end if
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
         end if

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         do i1=1, niterlong
           call parameterization_timestep
     .          ( i, time_initial, time, dt, fcor,
     .            thlm_forcing, rtm_forcing, wmm, wmt,
     .            wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .            p, rhom, rhot, exner,
     .            um, vm, upwp, vpwp, up2, vp2,
     .            thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .            rtp2, thlp2, rtpthlp,
     .            Scm, taum, rcm, Ncm, Ncnm, Nim,
     .            hydromet, cf, isValid, .false.
#ifdef SCALARS
     .            ,wpsclrp_sfc, wpedsclrp_sfc,
     .            sclrm, sclrm_forcing, edsclrm
#endif /*SCALARS*/
     .          )

#ifdef STATS
           call stats_end_timestep( time + dtmain, dtmain )
#endif /*STATS*/

           ! Set Time
           ! Advance time here, not in parameterization_timestep,
           ! in order to facilitate use of stats.
           ! A host model, e.g. WRF, would advance time outside
           ! of hoc_closure_timestep.  Vince Larson 7 Feb 2006
           if ( i1 < niterlong ) then
             time = time_initial + (i-1) * dtmain + i1 * dtclosure
           else if ( i1 == niterlong ) then
             time = time_initial + i * dtmain
           end if

           if ( .not. isValid ) exit
         end do ! i1=1..niterlong

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         if ( .not. isValid ) exit
       end do ! i=1, ifinal

!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

! Free memory

        call parameterization_cleanup( .false. )

#ifdef STATS
        call stats_finalize( )
#endif

        return
        end subroutine hoc_model_field_input

!-----------------------------------------------------------------------
        subroutine hoc_model_standard
     .            ( params, runfile, isValid, stdout )
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
!         isValid:  A false return value indicates something was unsound.

!-----------------------------------------------------------------------

        use grid_class
        use param_index, only: nparams
        use diagnostic_variables
        use prognostic_variables
        use numerical_check
        use hoc_parameterization_interface
        use constants, only: fstdout, fstderr
#ifdef STATS
        use hoc_stats
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

        ! Input Variables
        real, intent(in), dimension(nparams) :: 
     .  params  ! Model parameters, C1, nu2, etc.

        ! Subroutine Arguments (Model Setting)
        character(len=*), intent(in) :: 
     .  runfile ! Name of file containing &model_setting and &sounding

        logical, intent(in) :: 
     .  stdout  ! Whether to print output per timestep

        ! Output Variables
        logical, intent(inout) :: isValid    ! valid run? (T/F)

        ! Local Variables
        ! Internal Timing Variables
        integer ::
     .  ifinal,
     .  niterlong

        real :: 
     .  fcor, ! Coriolis parameter [s^-1]
     .  T0    ! Reference Temperature    [K]

        logical :: 
     .  cloud_sed,    ! Flag for cloud water droplet sedimentation. - Brian
     .  kk_rain,      ! Flag for Khairoutdinov and Kogan rain microphysics. - Brian
     .  licedfs,      ! Flag for simplified ice scheme
     .  lcoamps_micro,! Flag for COAMPS microphysical scheme
     .  lbugsrad,     ! Flag for BUGsrad radiation scheme
     .  lmpace_nudge, ! Whether to adjust the winds within the timestep
     .  lKhm_aniso    ! Whether to use anisotropic Khm.  - Michael Falk 2 Feb 2007

        ! Dummy grid height arrays
        real, dimension(:), allocatable :: 
     .  dummy_array_m, dummy_array_t

        integer :: i, i1 ! Internal Loop Variables

        integer :: iunit

        ! Definition of namelists
        namelist /model_setting/ 
     .  runtype, nzmax, deltaz, zm_init,
     .  day, month, year, rlat, rlon,
     .  time_initial, time_final, time_spinup,
     .  dtmain, dtclosure,
     .  sfctype, Tsfc, psfc, SE, LE, fcor, T0,
     .  cloud_sed, kk_rain, licedfs, lcoamps_micro, 
     .  lbugsrad, lKhm_aniso

!-----------------------------------------------------------------------

        ! Initialize the model run 

        ! Figure out which iounit to use
#ifdef _OPENMP
        iunit = omp_get_thread_num( ) + 10
#else
        iunit = 10
#endif

        ! Read namelist file
        open(unit=iunit, file=runfile, status='old')
        read(unit=iunit, nml=model_setting)
        close(unit=iunit)

       ! Dummy grid height arrays
       allocate( dummy_array_m(1:nzmax), dummy_array_t(1:nzmax) )

       dummy_array_m = 0.0
       dummy_array_t = 0.0

       if ( trim( runtype ) == "mpace" ) then 
         lmpace_nudge = .true.
       else
         lmpace_nudge = .false.
       end if

       ! Allocate & initialize variables,
       ! setup grid, setup constants, and setup flags

       call parameterization_setup
     .      ( nzmax, deltaz, zm_init, T0, params,
     .        lbugsrad, kk_rain, licedfs, lcoamps_micro,
     .        cloud_sed,
     .        lmpace_nudge, lKhm_aniso, .false.,
     .        dummy_array_m, dummy_array_t, isValid )

       if ( .not. isValid ) return


       ! Dummy grid height arrays
       deallocate( dummy_array_m, dummy_array_t )

       time = time_initial

       call hoc_initialize( iunit, runfile, psfc, thlm, rtm, 
     .                      um, vm, ug, vg, wp2,
     .                      rcm, hydromet,
     .                      wmt, wmm, em, exner,
     .                      taut, taum, thvm, p,
     .                      rhot, rhom, Lscale, lup,
     .                      ldown, Kht, Khm
#ifdef SCALARS
     .                      ,sclrm, edsclrm
#endif /*SCALARS*/
     .                     )

       um_init = um ! Michael Falk addition for nudging code.  27 Sep/1 Nov 2006
       vm_init = vm ! ditto

#ifdef STATS
#ifdef _OPENMP
       iunit = omp_get_thread_num( ) + 50
#else
       iunit = 50
#endif
!      Initialize statistics output
       call stats_init( iunit, runfile, gr%nnzp, gr%zt, gr%zm,
     .                  day, month, year, rlat, rlon, time, dtmain )
#endif /*STATS*/

       ! Time integration
       ! Call hoc_closure_timestep once per each GrADS output time 
       ifinal  = ( time_final - time_initial ) / dtmain


!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

       do i = 1, ifinal, 1

#ifdef STATS
         ! When this time step is over, the time will be time + dtmain

         ! We use elapsed time for stats_begin_step
         call stats_begin_timestep( time-time_initial+dtmain, dtmain )

         if ( lstats_last .and. stdout ) then
           write(unit=fstdout,fmt='(a,i7,a,f8.1)') 'iteration = ', 
     .       i, '; time = ', time
         end if
#endif /*STATS*/

         isValid = .not. invalid_model_arrays( ) ! Check for bad values 
                                                 ! in the model arrays
         if ( .not. isValid ) exit

         call hoc_forcings_timestep( i, stdout, isValid )

         if ( .not. isValid ) exit

         ! Compute number of iterations for closure loop
         if ( time > time_spinup ) then
           niterlong = 1
           dt        = dtmain
         else
           niterlong = int( dtmain / dtclosure )
           dt        = dtclosure
         end if

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         do i1=1, niterlong
           call parameterization_timestep
     .          ( i, time_initial, time, dt, fcor,
     .            thlm_forcing, rtm_forcing, wmm, wmt,
     .            wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .            p, rhom, rhot, exner,
     .            um, vm, upwp, vpwp, up2, vp2,
     .            thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .            rtp2, thlp2, rtpthlp,
     .            Scm, taum, rcm, Ncm, Ncnm, Nim,
     .            hydromet, cf, isValid, .false.
#ifdef SCALARS
     .            ,wpsclrp_sfc, wpedsclrp_sfc, 
     .            sclrm, sclrm_forcing, edsclrm
#endif /*SCALARS*/
     .          )

#ifdef STATS
           call stats_end_timestep( time + dtmain, dtmain )
#endif /*STATS*/

           ! Set Time
           ! Advance time here, not in parameterization_timestep,
           ! in order to facilitate use of stats.
           ! A host model, e.g. WRF, would advance time outside
           ! of hoc_closure_timestep.  Vince Larson 7 Feb 2006
           if ( i1 < niterlong ) then
             time = time_initial + (i-1) * dtmain + i1 * dtclosure
           else if ( i1 == niterlong ) then
             time = time_initial + i * dtmain
           end if

           if ( .not. isValid ) exit
         end do ! i1=1..niterlong

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         if ( .not. isValid ) exit
       end do ! i=1, ifinal

!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

! Free memory

        call parameterization_cleanup( .false. )

#ifdef STATS
        call stats_finalize( )
#endif

        return
        end subroutine hoc_model_standard


!-----------------------------------------------------------------------
        subroutine hoc_initialize
     .             ( iunit, runfile, psfc, thlm, rtm, um, vm,
     .               ug, vg, wp2, rcm, hydromet,
     .               wmt, wmm, em, exner, taut, taum, thvm,
     .               p, rhot, rhom, Lscale, lup,
     .               ldown, Kht, Khm
#ifdef SCALARS
     .               ,sclrm, edsclrm
#endif
     .             )

!       Description:
!       Execute the necessary steps for the initialization of the 
!       HOC model run. 

!       Calls: (all these are external)
!       subroutine read_sounding
!       subroutine hydrostatic (twice)
!       subroutine length

!       function rsat
!       function sat_rcm
!       function zt2zm (from grid_class)
!-----------------------------------------------------------------------

        use constants
        use parameters
        use grid_class
        use sounding
#ifdef STATS
        use hoc_stats
#endif /*STATS*/

        implicit none

        ! External
        real, external :: rsat
        real, external :: sat_rcm
        external :: length

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
     .  lup,            ! Upward mixing length          [m] 
     .  ldown,          ! Downward mixing length        [m] 
     .  Kht, Khm,       ! Eddy diffusivity              [m^2/s]
     .  taum, taut,     ! Dissipation time              [s]
     .  thvm            ! Virtual potential temperature [K]

        real, dimension(gr%nnzp,hydromet_dim), target, intent(out) :: 
     .  hydromet        ! Hydrometeor types             [units vary]

#ifdef SCALARS
!       Output (optional)
        real, dimension(gr%nnzp,sclr_dim), intent(out) :: 
     .  sclrm,  ! Standard passive scalar [units vary]
     .  edsclrm ! Eddy diffusivity passive scalar [units vary]
#endif

        ! Local Variables
        real, dimension(gr%nnzp) :: tmp1

        real, dimension(:), pointer ::
     .  rrm,        ! Rain water mixing ratio           [kg/kg]
     .  ricem,      ! Ice water mixing ratio            [kg/kg]
     .  rgraupelm,  ! Graupel water mixing ratio        [kg/kg] 
     .  rsnowm      ! Snow water mixing ratio           [kg/kg]

        real :: cloud_top_height ! [m]
        real :: emax

        integer :: k

!-----------------------------------------------------------------------

        ! Assign pointers to each hydrometeor field
        rrm       => hydromet(1:gr%nnzp,1)
        rsnowm    => hydromet(1:gr%nnzp,3)
        ricem     => hydromet(1:gr%nnzp,4)
        rgraupelm => hydromet(1:gr%nnzp,5)

        ! Read sounding information

        call read_sounding( iunit, thlm, rtm, um, vm, ug, vg, 
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
          end do

        case default
          ! Initial profile is non-saturated thlm or any type of theta.
          thlm = thlm - Lv/(Cp*exner) * rcm
#ifdef SCALARS
          ! The following case used the formulas for theta in order to test
          ! the scalar code.  However, it needs to be commented out because
          ! it does not apply to such things as CO2.
!          sclrm(:,SCLR_THETA) = sclrm(:,SCLR_THETA) - Lv/(Cp*exner)*rcm
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

!         wp2 = 2.0 / 3.0 * em
!<--- End reduction of initial sounding for stability 24 Jan 07

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
          ! 4150 + 1400 m is the top of the cloud in Nov11
          cloud_top_height = 1400. + gr%zm(1)
          do k=1,gr%nnzp
            if ( gr%zm(k) < cloud_top_height ) then
              em(k) = 0.1
            else
              em(k) = emin
            end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)
        ! End Vince Larson's change.

          wp2 = 2.0 / 3.0  * em

        ! Adam for June 25 altocumulus case.
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

        case ( "lba" )
          em = 0.1
          wp2 = 2./3. * em

        ! Michael Falk for MPACE Arctic Stratus case.
        case ( "mpace" )

          cloud_top_height = 1300. ! 1300 m is the cloud top in MPACE.  Michael Falk 17 Aug 2006
          emax = 1.0
          do k=1,gr%nnzp

            rrm(k) = 0.
            ricem(k) = 0.
            rgraupelm(k) = 0.
            rsnowm(k) = 0.

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

        call length( thvm, thlm, rtm, rcm, 
     .               em, p, exner, Lscale, lup, ldown )

        ! Dissipation time
        tmp1 = sqrt( max( emin, zm2zt( em ) ) )
        taut = min( Lscale / tmp1, taumax )
        taum = min( ( zt2zm( Lscale ) / sqrt( em ) ), taumax )
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
        Khm = c_K * zt2zm( Lscale ) * sqrt( em )
!        Khm = zt2zm( Kht )

       return
       end subroutine hoc_initialize

!----------------------------------------------------------------------
        subroutine hoc_forcings_timestep( iteration, stdout, isValid )

!       Description:
!       Calculate tendency and surface variables

!       Calls: (* = model case)
!         subroutines *_sfclyr
!         subroutines *_tndncy
!         subroutine sfc_thermo_fluxes
!         subroutine sfc_momentum_fluxes
!----------------------------------------------------------------------

        ! Modules to be included
        use model_flags
        use constants
        use parameters
        use grid_class
        use diagnostic_variables
        use prognostic_variables
        use arm, only: arm_tndcy, arm_sfclyr
        use astex, only: astex_tndcy, astex_sfclyr
        use atex, only: atex_tndcy, atex_sfclyr
        use bomex, only: bomex_tndcy, bomex_sfclyr
        use cobra, only: cobra_tndcy, cobra_sfclyr
        use dycoms2_rf01, only: dycoms2_rf01_tndcy, dycoms2_rf01_sfclyr
        use dycoms2_rf02, only: dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr
        use fire, only: fire_tndcy, sfc_momentum_fluxes,
     .                  sfc_thermo_fluxes
        use gabls2, only: gabls2_tndcy, gabls2_sfclyr
        use rico, only: rico_tndcy, rico_sfclyr
        use lba, only: lba_tndcy, lba_sfclyr
        use mpace, only: mpace_tndcy, mpace_sfclyr
        use nov11, only: nov11_altocu_tndcy ! Above surface
        use jun25, only: jun25_altocu_tndcy ! Above surface
        use wangara, only: wangara_tndcy, wangara_sfclyr

        implicit none

        ! External
#ifdef radoffline
        external :: bugsrad_hoc
#endif

        ! Input Variables
        integer, intent(in) ::
     .    iteration

        logical, intent(in) ::
     .    stdout

        ! Input/Output Variables
        logical, intent(inout) ::
     .    isValid

        ! Local Variables

        real, dimension(:), pointer :: 
     .    rrm,   ! Rain water mixing ratio              [kg/kg]
     .    ricem  ! Prisitine ice water mixing ratio     [kg/kg]

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
               Ncm(k) = 30.0 * (1.0 + exp( -gr%zt(k)/2000.0 )) * 1.e6
     .                  * rhot(k) 
             end if
           end do

         case( "arm" ) ! ARM Cu case
           call arm_tndcy( time, thlm_forcing, radht, rtm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing 
#endif /* SCALARS */
     .                   )

         case( "bomex" ) ! BOMEX Cu case
           call bomex_tndcy( time, wmt, wmm, radht,
     .                       thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                       ,sclrm_forcing 
#endif /* SCALARS */
     .                     )

         case( "fire" ) ! FIRE Sc case
           call fire_tndcy( time, rhot, rcm, exner, 
     .                      wmt, wmm, Frad, radht,
     .                      thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing 
#endif /* SCALARS */
     .                    )

         case( "wangara" ) ! Wangara dry CBL
           call wangara_tndcy( time, wmt, wmm, 
     .                         thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                        ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                       )

         case( "atex" ) ! ATEX case
           call atex_tndcy( time, time_initial, rtm, 
     .                      rhot, rcm, exner, wmt, wmm, Frad, radht,
     .                      thlm_forcing, rtm_forcing, isValid
#ifdef SCALARS
     .                     ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                    )

         case( "dycoms2_rf01" ) ! DYCOMS2 RF01 case
           call dycoms2_rf01_tndcy( time, rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, Frad, radht, thlm_forcing,
     .                              rtm_forcing, isValid  
#ifdef SCALARS
     .                             ,sclrm, sclrm_forcing  ! optional variables
#endif /* SCALARS */
     .                            )

         case( "astex_a209" ) ! ASTEX Sc case for K & K
           call astex_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                      ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                     )

         case( "dycoms2_rf02_do", ! DYCOMS2 RF02 case with drizzle only.
     .         "dycoms2_rf02_ds", ! DYCOMS2 RF02 case with drizzle and cloud sedimentation.
     .         "dycoms2_rf02_nd", ! DYCOMS2 RF02 case with no drizzle and no cloud sedimentation.
     .         "dycoms2_rf02_so" )! DYCOMS2 RF02 case with cloud water sedimentation only.
           ! Added time, time_initial, and Ncnm to arguments for testing
           ! the COAMPS microphysics with RF02 -dschanen 30 Nov 2006
           call dycoms2_rf02_tndcy
     .          ( time, time_initial, rhot, rhom, rtm, rcm, exner,
     .            wmt, wmm, thlm_forcing, rtm_forcing,
     .            Frad, radht, Ncm, Ncnm, isValid
#ifdef SCALARS
     .           ,sclrm_forcing  ! Optional variable
#endif /* SCALARS */
     .                            )

         case( "nov11_altocu" ) ! Nov. 11 Altocumulus case.
           call nov11_altocu_tndcy( time, time_initial, dt, rlat, rlon,
     .                              thlm, rcm, p, exner, rhot, rtm, wmt,
     .                              wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, Ncnm
#ifdef SCALARS
     .                             ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                            )

         case( "jun25_altocu" ) ! June 25 Altocumulus case.
           call jun25_altocu_tndcy( time, time_initial, rlat, dt, thlm,
     .                              rcm, p, exner, rhot, rtm, wmt,
     .                              wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht
#ifdef SCALARS
     .                             ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                            )

         case ( "lba" )
           call lba_tndcy( time, wmt, wmm, radht, 
     .                     thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                    ,sclrm_forcing
#endif /* SCALARS */
     .                    )

         case ( "mpace" ) ! MPACE arctic stratus case
           call mpace_tndcy( time, time_initial, dtmain, rlat, thlm,
     .                       exner, rhot, rtm, p, thvm, rcm,
     .                       wmt, wmm, thlm_forcing, rtm_forcing,
     .                       Ncnm, Ncm, Frad, radht
#ifdef SCALARS
     .                      ,sclrm_forcing
#endif /* SCALARS */
     .                    )

        ! Brian Griffin for COBRA CO2 case.
        case ( "cobra" )
           call cobra_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                      ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                     )

         case ( "rico" ) ! RICO case
           call rico_tndcy( time, time_initial, dtmain, exner,
     .                      rhot, rcm, kk_rain, wmt, wmm,
     .                      thlm_forcing, rtm_forcing, radht, Ncm
#ifdef SCALARS
     .                     ,sclrm_forcing
#endif /* SCALARS */
     .                    )

         case ( "gabls2" ) ! GABLS 2 case
           call gabls2_tndcy( time, time_initial, 
     .                      rhot, rcm, kk_rain, wmt, wmm,
     .                      thlm_forcing, rtm_forcing, radht, Ncm
#ifdef SCALARS
     .                     ,sclrm_forcing
#endif /* SCALARS */
     .                    )

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
     .                              upwp_sfc, vpwp_sfc )
          ! sfctype = 0  fixed sfc sensible and latent heat fluxes 
          !                   as given in hoc.in
          ! sfctype = 1  bulk formula: uses given surface temperature 
          !                   and assumes over ocean
          if ( sfctype == 0 ) then
            wpthlp_sfc = SE
            wprtp_sfc  = LE
#ifdef SCALARS
            wpsclrp(:,SCLR_THETA) = SE
            wpsclrp(:,SCLR_RT)    = LE
            wpedsclrp(:,SCLR_THETA) = SE
            wpedsclrp(:,SCLR_RT)    = LE
#endif /* SCALARS */
          else if ( sfctype == 1 ) then
            call sfc_thermo_fluxes( um(2), vm(2),
     .                              Tsfc, psfc, 
     .                              thlm(2), rtm(2),
     .                              wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                              ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                            )

          else
            write(unit=fstderr,fmt=*) 
     .        "Invalid value of sfctype = ", sfctype
            stop

          end if

        case( "arm" )
          call arm_sfclyr( time, gr%zt(2), 1.1, 
     .                     thlm(2), um(2), vm(2),
     .                     upwp_sfc, vpwp_sfc, 
     .                     wpthlp_sfc, wprtp_sfc, ustar
#ifdef SCALARS
     .                     ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                    )

        case( "bomex" ) 
          call bomex_sfclyr( um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc,
     .                       wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                       ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                     )

        case( "wangara" )
          call wangara_sfclyr( time, um(2), vm(2),
     .                         upwp_sfc, vpwp_sfc,
     .                         wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                         ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                       )

        case( "atex" )
          call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                      upwp_sfc, vpwp_sfc,
     .                      wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                     ,sclrm(2,:), wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                    )

        case( "dycoms2_rf01" ) 
          call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, 
     .                              exner(1), um(2), vm(2), 
     .                              thlm(2), rtm(2), 
     .                              rhom(1), upwp_sfc, vpwp_sfc, 
     .                              wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                             ,sclrm(2,:), wpsclrp_sfc, 
     .                              wpedsclrp_sfc
#endif /* SCALARS */
     .                             )
        case( "astex_a209" )
          call astex_sfclyr( rhom(1), um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc, wpthlp_sfc, 
     .                       wprtp_sfc
#ifdef SCALARS
     .                       ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                     )
        case( "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .        "dycoms2_rf02_nd", "dycoms2_rf02_so" )
          call dycoms2_rf02_sfclyr( um(2), vm(2),
     .                              upwp_sfc, vpwp_sfc,
     .                              wpthlp_sfc, wprtp_sfc, ustar
#ifdef SCALARS
     .                              ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                             )

        case( "nov11_altocu" )
          ! There are no surface momentum or heat fluxes
          ! for the Nov. 11 Altocumulus case.

        case( "jun25_altocu" )
          ! There are no surface momentum or heat fluxes
          ! for the Jun. 25 Altocumulus case.

         case ( "lba" )
           call lba_sfclyr( time, gr%zt(2), rhom(1),
     .                      thlm(2), um(2), vm(2),
     .                      upwp_sfc, vpwp_sfc, 
     .                      wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                     ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                     )

         case ( "mpace" )
           call mpace_sfclyr( rhom(1), um(2), vm(2), upwp_sfc, vpwp_sfc,
     .                        wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                        ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                      )

        ! Brian Griffin for COBRA CO2 case.
        case ( "cobra" )
          call cobra_sfclyr( time, gr%zt(2), rhom(1),
     .                       thlm(2), um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc,
     .                       wpthlp_sfc, wprtp_sfc, ustar
#ifdef SCALARS
     .                      ,wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                     )

         case ( "rico" )
           call rico_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                        gr%zt(2), 299.8, 101540., ! 299.8 K is the RICO SST; 101540 Pa is the sfc pressure.
!     .                        gr%zt(2), Tsfc, psfc,
     .                        upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                       ,sclrm(2,:), wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                      )

         case ( "gabls2" )
          call gabls2_sfclyr(time, time_initial, gr%zt(2), 97200.,
     .                      um(2), vm(2), thlm(2), rtm(2),
     .                      upwp_sfc, vpwp_sfc,
     .                      wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                     ,sclrm(2,:), wpsclrp_sfc, wpedsclrp_sfc
#endif /* SCALARS */
     .                    )

        case default
          write(unit=fstderr,fmt=*) 
     .      "Invalid value of runtype = ", runtype
          stop

        end select ! runtype

        if ( lbugsrad ) then
#ifdef radoffline /*This directive is needed for this to work at all.*/

          ! Assign pointers to rain and ice
          rrm   => hydromet(1:gr%nnzp,1)
          ricem => hydromet(1:gr%nnzp,4)

          call bugsrad_hoc( gr%zm(gr%nnzp), gr%nnzp, rlat, rlon,
     .                      day, month, year, time,
     .                      thlm, rcm, rtm, rrm, ricem, 
     .                      cf, p, exner, rhom, Tsfc,
     .                      radht, Frad, thlm_forcing )
#else
          stop "Cannot call BUGSrad with these compile options."
#endif
        end if ! lbugsrad

        return
        end subroutine hoc_forcings_timestep

        end module hoc
