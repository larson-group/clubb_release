!-----------------------------------------------------------------------
! $Id: hoc.F,v 1.95 2006-11-02 01:44:43 dschanen Exp $

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
     .  nzmax  ! Vertical extent [-]
        real, private :: 
     .  deltaz ! Change per grid level [m]

        integer, private :: 
     .  day, month, year ! Day of start of simulation

        real, private :: 
     .  rlat, ! Latitude  [Degrees North]
     .  rlon  ! Longitude [Degrees East]

        character(len=50), private :: 
     .  runtype ! String identifying the model case; e.g. bomex

        ! 0: fixed sfc sensible and latent heat fluxes as
        !    given in namelist
        ! 1: bulk formula: uses given surface temperature
        !    and assumes over ocean
        integer, private :: sfctype 

        real, private :: 
     .  time_initial, ! Time of start of simulation     [s]
     .  time_final,   ! Time end of simulation          [s]
     .  time_spinup,  ! Time end of spin up period      [s]
     .  time          ! Current time of simulation      [s]

        real, private :: 
     .  dtmain,     ! Main model timestep                      [s]
     .  dtclosure,  ! Closure model timestep                   [s]
     .  dt          ! Current model timestep (based on spinup) [s]

        contains
!-----------------------------------------------------------------------

        subroutine hoc_model_field_input
     .             ( C1_in, C1b_in, C1c_in, 
     .               C2_in, C2b_in, C2c_in, C4_in, C5_in, 
     .               C6rt_in, C6rtb_in, C6rtc_in,
     .               C6thl_in, C6thlb_in, C6thlc_in,
     .               C7_in, C7b_in, C7c_in,
     .               C8_in, C8b_in, C10_in, 
     .               C11_in, C11b_in, C11c_in, C12_in,
     .               nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .               beta_in, 
     .               gamma_coef_in, gamma_coefb_in, gamma_coefc_in, 
     .               c_K_in, lmin_coef_in,
     .               taumin_in, taumax_in, mu_in,
     .               runfile, isValid, stdout,
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
        use numerical_check
        use hoc_parameterization_interface
        use constants, only: fstdout, fstderr
#ifdef STATS
        use statistics
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
        real, intent(in) :: 
     .  C1_in, C1b_in, C1c_in,
     .  C2_in, C2b_in, C2c_in, C4_in, C5_in,
     .  C6rt_in, C6rtb_in, C6rtc_in,
     .  C6thl_in, C6thlb_in, C6thlc_in,
     .  C7_in, C7b_in, C7c_in,
     .  C8_in, C8b_in, C10_in,
     .  C11_in, C11b_in, C11c_in, C12_in,
     .  nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .  beta_in,
     .  gamma_coef_in, gamma_coefb_in, gamma_coefc_in,
     .  c_K_in, taumin_in, taumax_in, lmin_coef_in, mu_in

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
        logical, intent(out) :: isValid    ! valid run? (T/F)

        ! Local Variables
        ! Internal Timing Variables
        integer ::
     .  ifinal,
     .  niterlong,
     .  sample_time

        real :: fcor ! Coriolis parameter [s^-1]

        logical :: 
     .  cloud_sed,    ! Flag for cloud water droplet sedimentation. - Brian
     .  kk_rain,      ! Flag for Khairoutdinov and Kogan rain microphysics. - Brian
     .  lbugsrad,     ! Flag for BUGsrad radiation scheme - dschanen
     .  lmpace_nudge  ! Whether to adjust the winds within the timestep

        ! Flag for the lowest model level being at the surface (true)
        ! or not at the surface (false).  Example of case where it is not
        ! at the surface:  Nov. 11 altocumulus case.
        logical ::
     .  bottom_at_sfc 

        ! Dummy grid height arrays
        real, dimension(:), allocatable :: 
     .  dummy_array_m, dummy_array_t

       integer :: i, i1 ! Internal Loop Variables

       ! Definition of namelists
       namelist /model_setting/ runtype, nzmax, deltaz, 
     .   day, month, year, rlat, rlon,
     .   time_initial, time_final, time_spinup,
     .   dtmain, dtclosure,
     .   sfctype, Tsfc, psfc, SE, LE, fcor,
     .   cloud_sed, kk_rain, lbugsrad, bottom_at_sfc

       !  Initialize the model run 

       ! Initialize run to valid
       isValid = .true.

       ! Read namelist file
       open(unit=10, file=runfile, status='old')
       read(unit=10, nml=model_setting)
       close(unit=10)

       ! Model sanity Check

       ! When tuning budgets, we will allow these invalid results
!       if (( beta_in < 0.0 ) .or. ( beta_in > 3.0 )) then
!         print *, 'beta < 0 OR beta > 3'
!         isValid = .false. 
!       endif 

       if ( lmin_coef_in * deltaz < 4.0 ) then
         print *, 'lmin < 4.0'
         isValid = .false.
       end if
       if ( .not. isValid ) return

       ! Dummy grid height arrays
       allocate( dummy_array_m(1:nzmax), dummy_array_t(1:nzmax) )

       dummy_array_m = 0.0
       dummy_array_t = 0.0

       if ( trim( runtype ) == "mpace" ) then 
         lmpace_nudge = .true.
       else
         lmpace_nudge = .false.
       end if

       ! Allocate & initilialize variables,
       ! setup grid, setup constants, and setup flags
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
     .        lbugsrad, kk_rain, cloud_sed, bottom_at_sfc,
     .        lmpace_nudge, .false., dummy_array_m, dummy_array_t )


       ! Dummy grid height arrays
       deallocate( dummy_array_m, dummy_array_t )

       time = time_initial

       call hoc_initialize( runfile, psfc, thlm, rtm, 
     .                      um, vm, ug, vg, wp2,
     .                      rcm, rrm, rim, rgm, rsnowm,
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

         call stats_begin_timestep( time + dtmain, dtmain )

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
     .          ( i1, niterlong, time_initial, time, dt,
     .            thlm_forcing, rtm_forcing, wmm, wmt,
     .            wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .            p, rhom, rhot, exner,
     .            um, vm, upwp, vpwp,
     .            thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .            rtp2, thlp2, rtpthlp,
     .            Scm, taum, rcm, Ncm, rrm, Nrm, cf,
     .            isValid, .false.
#ifdef SCALARS
     .            ,sclrm, sclrm_forcing,
     .            edsclrm, edsclrmt
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
     .            ( C1_in, C1b_in, C1c_in,
     .              C2_in, C2b_in, C2c_in, C4_in, C5_in, 
     .              C6rt_in, C6rtb_in, C6rtc_in,
     .              C6thl_in, C6thlb_in, C6thlc_in,
     .              C7_in, C7b_in, C7c_in,
     .              C8_in, C8b_in, C10_in, 
     .              C11_in, C11b_in, C11c_in, C12_in,
     .              nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in, 
     .              beta_in, 
     .              gamma_coef_in, gamma_coefb_in, gamma_coefc_in,
     .              c_K_in, lmin_coef_in, 
     .              taumin_in, taumax_in, mu_in,
     .              runfile, isValid, stdout )
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
        use diagnostic_variables
        use prognostic_variables
        use numerical_check
        use hoc_parameterization_interface
        use constants, only: fstdout, fstderr
#ifdef STATS
        use statistics
#endif /*STATS*/

        implicit none

        ! External
        intrinsic :: mod, real, int 

        ! Input Variables
        real, intent(in) :: 
     .  C1_in, C1b_in, C1c_in,
     .  C2_in, C2b_in, C2c_in, C4_in, C5_in,
     .  C6rt_in, C6rtb_in, C6rtc_in,
     .  C6thl_in, C6thlb_in, C6thlc_in,
     .  C7_in, C7b_in, C7c_in,
     .  C8_in, C8b_in, C10_in,
     .  C11_in, C11b_in, C11c_in, C12_in,
     .  nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .  beta_in,
     .  gamma_coef_in, gamma_coefb_in, gamma_coefc_in,
     .  c_K_in, taumin_in, taumax_in, lmin_coef_in, mu_in

        ! Subroutine Arguments (Model Setting)
        character(len=*), intent(in) :: 
     .  runfile ! Name of file containing &model_setting and &sounding

        logical, intent(in) :: 
     .  stdout  ! Whether to print output per timestep

        ! Output Variables
        logical, intent(out) :: isValid    ! valid run? (T/F)

        ! Local Variables
        ! Internal Timing Variables
        integer ::
     .  ifinal,
     .  niterlong

        real :: fcor ! Coriolis parameter [s^-1]

        logical :: 
     .  cloud_sed,    ! Flag for cloud water droplet sedimentation. - Brian
     .  kk_rain,      ! Flag for Khairoutdinov and Kogan rain microphysics. - Brian
     .  lbugsrad,     ! Flag for BUGsrad radiation scheme - dschanen
     .  lmpace_nudge  ! Whether to adjust the winds within the timestep

        ! Flag for the lowest model level being at the surface (true)
        ! or not at the surface (false).  Example of case where it is not
        ! at the surface:  Nov. 11 altocumulus case.
        logical ::
     .  bottom_at_sfc 

        ! Dummy grid height arrays
        real, dimension(:), allocatable :: 
     .  dummy_array_m, dummy_array_t

        integer :: i, i1 ! Internal Loop Variables

        ! Definition of namelists
        namelist /model_setting/ runtype, nzmax, deltaz, 
     .  day, month, year, rlat, rlon,
     .  time_initial, time_final, time_spinup,
     .  dtmain, dtclosure,
     .  sfctype, Tsfc, psfc, SE, LE, fcor,
     .  cloud_sed, kk_rain, lbugsrad, bottom_at_sfc

!-----------------------------------------------------------------------

        !  Initialize the model run 

        ! Initialize run to valid
        isValid = .true.

        ! Read namelist file
        open(unit=10, file=runfile, status='old')
        read(unit=10, nml=model_setting)
        close(unit=10)

       ! Model sanity Check

        if (( beta_in < 0.0 ) .or. ( beta_in > 3.0 )) then
          print *, 'beta < 0 OR beta > 3'
          isValid = .false. 
        endif 

        if ( lmin_coef_in * deltaz < 4.0 ) then
          print *, 'lmin < 4.0'
          isValid = .false.
        end if
        if ( .not. isValid ) return

       ! Dummy grid height arrays
       allocate( dummy_array_m(1:nzmax), dummy_array_t(1:nzmax) )

       dummy_array_m = 0.0
       dummy_array_t = 0.0

       ! Allocate & initilialize variables,
       ! setup grid, setup constants, and setup flags

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
     .        lbugsrad, kk_rain, cloud_sed, bottom_at_sfc,
     .        lmpace_nudge, .false., dummy_array_m, dummy_array_t )


       ! Dummy grid height arrays
       deallocate( dummy_array_m, dummy_array_t )

       time = time_initial

       call hoc_initialize( runfile, psfc, thlm, rtm, 
     .                      um, vm, ug, vg, wp2,
     .                      rcm, rrm, rim, rgm, rsnowm,
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
!      Initialize statistics output
       call stats_init( 50, runfile, gr%nnzp, gr%zt, gr%zm,
     .                  day, month, year, rlat, rlon, time, dtmain )
#endif /*STATS*/

       ! Time integration
       ! Call hoc_closure_timestep once per each GrADS output time 
       ifinal  = ( time_final - time_initial ) / dtmain


!<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<

       do i = 1, ifinal, 1

#ifdef STATS
         ! When this time step is over, the time will be time + dtmain

         call stats_begin_timestep( time + dtmain, dtmain )

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
     .          ( i1, niterlong, time_initial, time, dt,
     .            thlm_forcing, rtm_forcing, wmm, wmt,
     .            wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc,
     .            p, rhom, rhot, exner,
     .            um, vm, upwp, vpwp,
     .            thlm, rtm, wprtp, wpthlp, wp2, wp3,
     .            rtp2, thlp2, rtpthlp,
     .            Scm, taum, rcm, Ncm, rrm, Nrm, cf,
     .            isValid, .false.
#ifdef SCALARS
     .            ,sclrm, sclrm_forcing,
     .            edsclrm, edsclrmt
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
     .             ( runfile, psfc, thlm, rtm, um, vm,
     .               ug, vg, wp2, rcm, rrm, rim,
     .               rgm, rsnowm, wmt, wmm,
     .               em, exner, taut, taum, thvm,
     .               p, rhot, rhom, Lscale, lup,
     .               ldown, Kht, Khm, 
     .               sclrm, edsclrm )

!       Description:
!       Execute the necessary steps for the initialization of the 
!       HOC model run. 

!       Calls: (all these are external)
!       subroutine read_sounding
!       subroutine hydrostatic (twice)
!       subroutine sfc_var
!       subroutine length

!       function rsat
!       function sat_rcm
!       function zt2zm (from grid_class)
!-----------------------------------------------------------------------

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

        intrinsic :: min, max, trim, sqrt, size

        ! Input
        character(len=*), intent(in) :: runfile  ! filename for the namelist

        real, intent(in) :: psfc ! Pressure at the surface [Pa]

!       ! Output
        real, dimension(gr%nnzp), intent(out) :: 
     .  thlm,           ! Theta l mean                  [K] 
     .  rtm,            ! Total water mixing ratio      [kg/kg]
     .  um,             ! u wind                        [m/s]
     .  vm,             ! v wind                        [m/s]
     .  ug,             ! u geostrophic wind            [m/s] 
     .  vg,             ! u geostrophic wind            [m/s] 
     .  wp2,            ! w'^2                          [m^2/s^2]
     .  rcm,            ! Cloud water mixing ratio      [kg/kg]
     .  rrm,            ! Rain water mixing ratio       [kg/kg] Michael Falk, 1 Nov 2006
     .  rim,            ! Ice water mixing ratio        [kg/kg]           "
     .  rgm,            ! Graupel water mixing ratio    [kg/kg]           "
     .  rsnowm,         ! Snow water mixing ratio       [kg/kg]           "
     .  wmt, wmm,       ! w wind                        [m/s]
     .  exner,          ! Exner function                [-] 
     .  em,             ! Em                            [m^2/s^2]
     .  p,              ! Pressure                      [Pa]
     .  rhot, rhom,     ! Density                       [kg/m^3]
     .  Lscale,         ! Mixing length                 [m] 
     .  lup,            ! Upward mixing length          [m] 
     .  ldown,          ! Downward mixing length        [m] 
     .  Kht, Khm,       ! Eddy diffusivity              [m^2/s]
     .  taum, taut,     ! Dissipation time              [s]
     .  thvm            ! Virtual potential temperature [K]

!       Output (optional)
        real, optional, dimension(gr%nnzp,sclr_dim), intent(out) :: 
     .  sclrm,  ! Standard passive scalar [units vary]
     .  edsclrm ! Eddy diffusivity passive scalar [units vary]

        ! Local Variables
        real, dimension(gr%nnzp) :: tmp1

        real :: cloud_top_height ! [m]
        real :: emax

        integer :: k

!-----------------------------------------------------------------------

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
          end do

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
             end if
          end do
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
             end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

          wp2 = 2.0 / 3.0 * em

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

          do k = 1, gr%nnzp
             em(k) = 1.0 - (gr%zm(k)/3000.0)
             if ( em(k) < emin ) then
                em(k) = emin
             end if
          end do
          em(1) = em(2)
          em(gr%nnzp) = em(gr%nnzp-1)

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
          cloud_top_height = 1400. ! 1400 m is the top of the cloud in Nov11
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

        ! Michael Falk for MPACE Arctic Stratus case.
        case ( "mpace" )

          cloud_top_height = 1300. ! 1300 m is the cloud top in MPACE.  Michael Falk 17 Aug 2006
          emax = 1.0
          do k=1,gr%nnzp

            rrm(k) = 0.
            rim(k) = 0.
            rgm(k) = 0.
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
!         subroutine sfc_var
!----------------------------------------------------------------------

        ! Modules to be included
        use model_flags
        use constants
        use grid_class
        use diagnostic_variables
        use prognostic_variables
        use gcss
        use surface

        implicit none

        ! Input Variables
        integer, intent(in) :: iteration
        logical, intent(in) :: stdout

        ! Input/Output Variables
        logical, intent(inout) :: isValid

        ! Michael Falk, 1 November 2006, thinks these should be renamed
        ! and moved to diag_variables.F .
        real snowv(1,1,gr%nnzp)
        real falli(1,1,gr%nnzp)
        real fallg(1,1,gr%nnzp)
        real fallr(1,1,gr%nnzp)
        ! eMFc

!-----------------------------------------------------------------------

!#######################################################################
!##############      FIND ALL DIAGNOSTIC VARIABLES        ##############
!#######################################################################

         select case ( runtype )
         case( "normal" ) ! Normal run

         case( "arm" ) ! ARM Cu case
           call arm_tndcy( time, thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing 
#endif /* SCALARS */
     .                   )

         case( "bomex" ) ! BOMEX Cu case
           call bomex_tndcy( time, wmt, wmm, 
     .                       thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                       ,sclrm_forcing 
#endif /* SCALARS */
     .                     )

         case( "fire" ) ! FIRE Sc case
           call fire_tndcy( time, wmt, wmm, thlm_forcing, 
     .                      rtm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing 
#endif /* SCALARS */
     .                    )
           thlm_forcing = 0.0
#ifdef SCALARS
           sclrm_forcing(:, SCLR_THETA) = 0.0
#endif /* SCALARS */

           call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing(:,SCLR_THETA)
#endif /* SCALARS */
     .                   )

         case( "wangara" ) ! Wangara dry CBL
           call wangara_tndcy( time, wmt, wmm, 
     .                         thlm_forcing, rtm_forcing
#ifdef SCALARS
     .                        ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                       )

         case( "atex" ) ! ATEX case
           call atex_tndcy( time, time_initial, rtm, 
     .                      wmt, wmm, thlm_forcing,
     .                      rtm_forcing, isValid
#ifdef SCALARS
     .                     ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                    )

           call cloud_rad( rhot, rcm, exner, Frad, radht, thlm_forcing
#ifdef SCALARS
     .                     ,sclrm_forcing(:, SCLR_THETA)
#endif /* SCALARS */
     .                   )

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
     .                             ,sclrm_forcing   ! optional variable
#endif /* SCALARS */
     .                     )

         case( "dycoms2_rf02_do" ) ! DYCOMS2 RF02 case with drizzle only.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, Ncm, isValid
#ifdef SCALARS
     .                              ,sclrm_forcing  ! optional variable
#endif /* SCALARS */
     .                            )

         case( "dycoms2_rf02_ds" ) ! DYCOMS2 RF02 case with drizzle and cloud sedimentation.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, Ncm, isValid
#ifdef SCALARS
     .                              ,sclrm_forcing  ! optional variable
#endif /* SCALARS */
     .                            )

         case( "dycoms2_rf02_nd" ) ! DYCOMS2 RF02 case with no drizzle and no cloud sedimentation.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, Ncm, isValid
#ifdef SCALARS
     .                              ,sclrm_forcing  ! optional variable
#endif /*SCALARS*/
     .                            )

         case( "dycoms2_rf02_so" ) ! DYCOMS2 RF02 case with cloud water sedimentation only.
           call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, exner,
     .                              wmt, wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht, Ncm, isValid
#ifdef SCALARS
     .                             ,sclrm_forcing   ! optional variable
#endif /*SCALARS*/
     .                            )

         case( "nov11_altocu" ) ! Nov. 11 Altocumulus case.
           call nov11_altocu_tndcy( time, time_initial, dt, thlm,
     .                              rcm, p, exner, rhot, rtm, wmt,
     .                              wmm, thlm_forcing, rtm_forcing,
     .                              Frad, radht
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

         case ( "mpace" ) ! MPACE case
           call mpace_tndcy( time, time_initial, dtmain, thlm,
     .                       rim, rrm, rgm, rsnowm,
     .                       exner, rhot, rtm, p, thvm, rcm,
     .                       wmt, wmm, thlm_forcing, rtm_forcing,
     .                       Frad, radht, 
     .                       nc3, ni3, nr3, ncn3,
     .                       snowv, falli, fallg, fallr
#ifdef SCALARS
     .                     ,sclrm_forcing
#endif /* SCALARS */
     .                    )

#ifdef SCALARS
           sclrm_forcing(:, SCLR_THETA) = 0.
#endif /* SCALARS */

         case default
           write(unit=fstderr,fmt=*) 
     .       "hoc_forcings_timestep: Don't know how to handle "
     .       //"LS forcing for runtype: "//trim( runtype )
           stop

         end select

         ! Bc for the second order moments

        select case ( trim( runtype ) )

        case( "normal", "fire" )  ! Normal and FIRE
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
     .                     ,wpsclrp(1,:)
#endif /* SCALARS */
     .                    )

        case( "bomex" ) 
          call bomex_sfclyr( um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc,
     .                       wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                       ,wpsclrp(1,:)
#endif /* SCALARS */
     .                     )

        case( "wangara" )
          call wangara_sfclyr( time, um(2), vm(2),
     .                         upwp_sfc, vpwp_sfc,
     .                         wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                         ,wpsclrp(1,:)
#endif /* SCALARS */
     .                       )

        case( "atex" )
          call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                      upwp_sfc, vpwp_sfc,
     .                      wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                     ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                    )

        case( "dycoms2_rf01" ) 
          call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, 
     .                              exner(1), um(2), vm(2), 
     .                              thlm(2), rtm(2), 
     .                              rhom(1), upwp_sfc, vpwp_sfc, 
     .                              wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                              ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                             )
        case( "astex_a209" )
          call astex_sfclyr( rhom(1), um(2), vm(2),
     .                       upwp_sfc, vpwp_sfc, wpthlp_sfc, 
     .                       wprtp_sfc
#ifdef SCALARS
     .                       ,wpsclrp(1,:)
#endif /* SCALARS */
     .                     )
        case( "dycoms2_rf02_do", "dycoms2_rf02_ds",
     .        "dycoms2_rf02_nd", "dycoms2_rf02_so" )
          call dycoms2_rf02_sfclyr( um(2), vm(2),
     .                              upwp_sfc, vpwp_sfc,
     .                              wpthlp_sfc, wprtp_sfc, ustar
#ifdef SCALARS
     .                              ,sclrm(2,:), wpsclrp(1,:)
#endif /* SCALARS */
     .                             )

        case( "nov11_altocu" )
          ! There are no surface momentum or heat fluxes
          ! for the Nov. 11 Altocumulus case.

        case( "jun25_altocu" )
          ! There are no surface momentum or heat fluxes
          ! for the Jun. 25 Altocumulus case.

         case ( "mpace" )
           call mpace_sfclyr( rhom(1), um(2), vm(2), upwp_sfc, vpwp_sfc,
     .                        wpthlp_sfc, wprtp_sfc
#ifdef SCALARS
     .                        ,wpsclrp(1,:)
#endif /* SCALARS */
     .                      )

        case default
          write(unit=fstderr,fmt=*) 
     .      "Invalid value of runtype = ", runtype
          stop

        end select ! runtype

        if ( lbugsrad ) then
#ifdef radoffline /*This directive is needed for this to work at all.*/
          call bugsrad_hoc( deltaz*gr%nnzp, gr%nnzp, rlat, rlon,
     .                      day, month, year, time,
     .                      thlm, rcm, rtm,
     .                      rrm, cf, p, exner, rhom, Tsfc,
     .                      radht, Frad, thlm_forcing )
#else
          stop "Cannot call BUGSrad with these compile options."
#endif
        end if ! lbugsrad

        return
        end subroutine hoc_forcings_timestep

        end module hoc
