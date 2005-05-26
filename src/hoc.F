       module hoc

       contains
! Program to integrate the pde equations for pdf closure

       subroutine hoc_sub (C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                     C8_in, C10_in, C11_in, 
     .                     nu1_in, nu2_in, nu6_in, nu8_in, beta_in,
     .                     gamma_coef_in, c_K_in, 
     .                     lmin_coef_in, taumin_in, taumax_in, mu_in,
     .                     runfile, run_name, isValid, stdout)

! Modules to be included

        use constants
        use grid_class
        use model_arrays
        use outputfile_class
#ifdef STATS
        use statistics
#endif STATS

        implicit none

! External functions
        
        real rsat
        external rsat

! Model setting
        character(len=*), intent(in) ::  runfile       ! Name of the nml file

        integer runtype                 ! = 0 normal run
                                        ! = 1 GCSS ARM
                                        ! = 2 GCSS BOMEX
                                        ! = 3 GCSS FIRE
                                        ! = 4 Wangara day 33
                                        ! = 6 GCSS ATEX
                                        ! = 7 GCSS DYCOMS II (RF01)
                                        ! = 8 GCSS ASTEX with K & K changes 
                                        ! = 9 GCSS DYCOMS II (RF02)
                                        ! = 10 GCSS DYCOMS II (RF02) no drizzle

        integer sfctype                 ! = 0  fixed sfc sensible and
                                        !      latent heat fluxes as
                                        !      given in hoc.in
                                        ! = 1  bulk formula: uses
                                        !      given surface temperature
                                        !      and assumes over ocean

        integer day, month, year        ! day of start ot simulation
        real time_initial               ! time (s) of start of simulation
        real time_final                 ! time (s) end of simulation
        real time_spinup                ! time (s) end of spin up period
        real time                       ! actual time of simulation

        real dtmain                     ! host model timestep (s)
        real dtclosure                  ! closure model timestep (s)
        real dt
        real dtoutput                   ! interval of file output (s)


! Grid definition

        integer nzmax
        real deltaz

! Subroutine Arguments (constants)
        real, intent(in) :: C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                      C8_in, C10_in, C11_in
        real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in
        real, intent(in) :: beta_in
        real, intent(in) :: gamma_coef_in
        real, intent(in) :: c_K_in, taumin_in, taumax_in, lmin_coef_in
        real, intent(in) :: mu_in

! Variables

        integer i, i1, i2, k
        integer ifinal
        integer ioutput

        integer niterlong
        logical spunup                  ! Set but never used anywhere
                                        ! dschanen 4/4/2005

! Surface data

        real Tsfc                       ! surface temperature (K)
        real psfc                       ! surface pressure (Pa)
        real SE                         ! sensible heat flux (K/s)
        real LE                         ! latent heat flux (1/s)
        real ustar                      ! - Brian DYCOMS_2

! Flag for Khairoutdinov and Kogan rain microphysics. - Brian

        logical                :: kk_rain

! Misc
        real n2

        logical, intent(inout) :: isValid ! " " valid run ? (T/F)
        logical, intent(in)    :: stdout  ! whether to print to the terminal

! Definition of namelists
        namelist /model_setting/ runtype, nzmax, deltaz, 
     .                           day, month, year, 
     .                           time_initial, time_final, time_spinup,
     .                           dtmain, dtclosure, dtoutput, 
     .                           sfctype, Tsfc, psfc, SE, LE, fcor

! GrADS file variables

        character(len=200) fdir,fname
        type (outputfile) fzt, fzm, fpdf, fsfc
        character(len=*), intent(in) :: run_name

! Local variables for Latin Hypercube sampling
!        logical, parameter ::  LH_on = .true. ! moved to constants.F

        real crt1, crt2, cthl1, cthl2   ! coeffs of s from pdf_closure_new

        integer i_rmd 

! Number of variables to sample
        integer, parameter :: d = 5

! n = number of calls to microphysics per timestep (normally=2)
        integer, parameter :: n = 12

! sequence_length = nt/n = number of timesteps before sequence repeats.
        integer, parameter :: sequence_length = 1

! nt = number of random samples before sequence of repeats (normally=10)
        integer nt

! Sample drawn from uniform distribution
        double precision :: X_u(1:n,1:(d+1))

! Sample that is transformed ultimately to normal-lognormal
        double precision :: X_nl(1:n,1:d)

! A true/false flag that determines whether
!     the PDF allows us to construct a sample
        logical :: sample_flag

        nt = n*sequence_length        
! End Latin hypercube section


! Read namelist file
        open(10, file = runfile, status = 'old' )
        read(10, nml = model_setting )
        close(10)

! Define constants
      call setup_constants(deltaz, C1_in, C2_in, C4_in, C5_in, C6_in, 
     .                     C7_in, C8_in, C10_in, C11_in,  
     .                     nu1_in, nu2_in, nu6_in, nu8_in,
     .                     beta_in, lmin_coef_in, taumin_in, taumax_in)

! Allocate memory for dynamic arrays and initialize them
! Vince Larson fed in LH vars.
!        call allocate_model_arrays( nzmax )
        call allocate_model_arrays( nzmax, n, nt, d ) 
! End Vince Larson's change

        call initialize_model_arrays()

! Define Sc
!
!        Sc = 0.4  ! Commented out for DYCOMS_2

! Set grid up

!        call gridsetup( nzmax, deltaz )
        call gridsetup( nzmax, deltaz, runtype )

! Read sounding information

        call read_sounding(thlm, rtm, um, vm, ug, vg, runfile, runtype)

! At this point, thlm actually contains theta (except for DYCOMS). 
! We need to compute liquid water content, and initilialize thlm properly

! First, compute approximate pressure using theta

        call hydrostatic(thlm, psfc, p, exner, rhot, rhom)

! Second, use this pressure to compute liquid water
! from excess saturation

        do k=1,gr%nnzp
           rcm(k) = max( rtm(k) - rsat(p(k),thlm(k)*exner(k)), 0. )
        end do

! Compute initial theta-l

        select case (runtype)
          case (7, 8, 9, 10 ) 
             ! Do nothing

          case default
            thlm = thlm - Lv/(Cp*exner) * rcm

        end select

! Now, compute initial thetav

        thvm = thlm + ep1 * T0 * rtm 
     .              + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

! Recompute more accurate initial exner function and pressure using thvm

        call hydrostatic(thvm, psfc, p, exner, rhot, rhom)

! Initialize imposed w

        do k=1,gr%nnzp
           wmt(k) = 0.0
           wmm(k) = 0.0
        end do

! Initialize TKE and other fields as needed
        select case ( runtype )
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

            wp2 = 2./3. * em

          case ( 6 )
            um = max( um, -8. )

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

            wp2 = 2./3. * em

      ! Brian for DYCOMS II RF02
          case ( 9, 10 ) 
            em = 1.0
            wp2 = (2.0/3.0) * em

        end select

! End Initialize TKE and other fields as needed

!       compute mixing length

        call length(thvm, thlm, rtm, rcm, 
     .              em, p, exner, Lscale, lup, ldown, mu_in)

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
!       c_K_in is 0.548 usually (Duynkerke and Driedonks 1987)

        Kht = c_K_in * Lscale * tmp1
        Khm = zt2zm( Kht )

!       surface values of second order moments

        call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),
     .                wp2(1), thlp2(1), rtp2(1), rtpthlp(1) )
       
! Time integration

        ifinal = ( time_final - time_initial )  / dtmain
        ioutput = dtoutput / dtmain

        time = time_initial
        spunup = .false.

! Initialize run to valid
        isValid = .true.
#ifdef STATS
! Initialize statistics output
        call stats_init(50, runfile, gr%nnzp, gr%zt, gr%zm,
     .                  day, month, year, time, dtmain)

#endif STATS

!<<<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<<<

!        do i=1, ifinal
        do i = 0, ifinal   ! Brian

          time = time_initial + i * dtmain

!           Model sanity Check
          if (( beta .lt. 0. ) .or. ( beta .gt. 3.0 )) then
            print *, 'beta < 0 OR beta > 3'
            isValid = .false. 
          endif 

          if (lmin < 4.0) then
            print *, 'lmin < 4.0'
            isValid = .false.
          endif

!         Check for a Not-A-Number or Infinity in any model array 
          if ( isValid ) isValid = Test_Variables( )

          if ( .not. isValid ) then
             call deallocate_model_arrays( ) 
!             deallocate( fzt%var, fzt%z, fzm%var, fzm%z, 
!     .                   fpdf%var, fpdf%z, fsfc%var, fsfc%z )
             return
          endif

!          print*,'hoc.F: just before iteration printout'

!          output to file every so often if called for
#ifdef STATS
           call stats_begin_timestep( time, dtmain )
           if ( lstats_last .and. stdout ) then
             write(*,'(a,i7,a,f8.1)') 'iteration = ',i,'; time = ',time
           endif
!           print*, 'hoc.F: just after iteration printout'
#endif STATS
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
     .                           thlm_forcing, rtm_forcing )
               kk_rain = .false.

             case ( 3 ) ! FIRE Sc case
               call fire_tndcy( time, wmt, wmm, thlm_forcing, 
     .                          rtm_forcing)
               thlm_forcing = 0.

               call cloud_rad( rhot, rcm, Frad, thlm_forcing )
               kk_rain = .false.

             case ( 4 ) ! Wangara dry CBL
               call wangara_tndcy( time, wmt, wmm, 
     .                             thlm_forcing, rtm_forcing )
               kk_rain = .false.

             case ( 6 ) ! ATEX case
               call atex_tndcy( time, rtm, wmt, wmm, thlm_forcing,
     .                          rtm_forcing, isValid )
               call cloud_rad( rhot, rcm, Frad, thlm_forcing )
               kk_rain = .false.

             case ( 7 ) ! DYCOMS2 RF01 case
               call dycoms2_rf01_tndcy( time, rhot, rhom, rtm, rcm,
     .                            wmt, wmm, Frad,
     .                            thlm_forcing, rtm_forcing, 
     .                            isValid )
               kk_rain = .false.

             case ( 8 ) ! ASTEX Sc case for K & K
               call astex_tndcy( time, wmt, wmm, 
     .                           thlm_forcing, rtm_forcing )
               kk_rain = .true.

             case ( 9 ) ! DYCOMS2 RF02 case with drizzle
               call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, wmt, wmm,
     .                                  thlm_forcing, rtm_forcing, 
     .                                  isValid ) 
               kk_rain = .true.

             case ( 10 ) ! DYCOMS2 RF02 case with no drizzle
               call dycoms2_rf02_tndcy( rhot, rhom, rtm, rcm, wmt, wmm,
     .                                   thlm_forcing, rtm_forcing, 
     .                                   isValid ) 
               kk_rain = .false.

             case default
               stop "HOC_SUB: Don't know how to do this runtype"

           end select

!           print*, 'hoc.F: got past runtype selection'

!           Added to FIRE case above -dschanen 4/08/2005
!           if ( runtype == 3 ) then
!             thlm_forcing = 0.
!             call cloud_rad(rhot,rcm,Frad,thlm_forcing)
!          end if

! Compute mixing length

           call length(thvm, thlm, rtm, rcm, 
     .                 em, p, exner, Lscale, lup, ldown, mu_in)

! Compute Sc with new formula from Vince

!  Linear formula
!         Sc = 0.6 * ( 1.0 - 
!     .    min (
!     .    max ( abs(wpthlp/(sqrt(wp2)*sqrt(thlp2)+0.01*wtol*thltol)) , 
!     .          abs(wprtp/(sqrt(wp2)*sqrt(rtp2)+0.01*wtol*rttol)) ) ,
!     .    1.0      ) )

!  Quadratic formula
!  note that gamma_coef_in is typically 0.45 
!  Divide Sc into Scm and Sct.
        Scm = gamma_coef_in * ( 1.0 - 
     .    min (
     .    max ( (wpthlp/(sqrt(wp2)*sqrt(thlp2)+0.01*wtol*thltol))**2 , 
     .          (wprtp/(sqrt(wp2)*sqrt(rtp2)+0.01*wtol*rttol))**2 ),
     .    1.0      ) )

!        Scm = 0.45 * ( 1.0 - 
!     .    min (
!     .    max ( (wpthlp/(sqrt(wp2)*sqrt(thlp2)+0.01*wtol*thltol))**2 , 
!     .          (wprtp/(sqrt(wp2)*sqrt(rtp2)+0.01*wtol*rttol))**2 ),
!     .    1.0      ) )

        Sct = zm2zt(Scm)


! Compute number of iterations for closure loop

           if ( time > time_spinup ) then
              spunup = .true.
              niterlong = 1
              dt = dtmain
           else
              spunup = .false.
              niterlong = int( dtmain / dtclosure )
              dt = dtclosure
           endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<


           do i1=1, niterlong
!             Diagnose variances

              call diag_var( taum, 
     .                       wmm, rtm, wprtp, thlm, wpthlp,
     .                       wp2, wp3, Scm, 
     .                       rtp2, thlp2, rtpthlp,
     .                       .false., dtclosure, isValid)

!             Check for degeneratation 
              if ( .not. isValid ) then
                call deallocate_model_arrays( )
!                deallocate( fzt%var, fzt%z, fzm%var, fzm%z, 
!     .                      fpdf%var, fpdf%z, fsfc%var, fsfc%z )
                return
              endif

!              print*, 'hoc.F: got past deallocate_model_arrays'

!             Dissipation time

              tmp1 = sqrt( zm2zt( em ) )
              taut = min( Lscale / tmp1, taumax )

!             Modification to damp noise in stable region
              do k=1,gr%nnzp
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

!             Eddy diffusivity coefficient
!             c_K_in is 0.548 usually (Duynkerke and Driedonks 1987)

              Kht = c_K_in * Lscale * tmp1
              Khm = zt2zm( Kht )

!             Advance wp2/wp3 one timestep

              call timestep_wp23
     .        (dtclosure, Lscale, taum, taut, Sct, pdf_parms(:,13), Khm, ! Brian
     .         wmt, wmm,
     .         um, upwp, vm, vpwp, wpthvp, wp2thvp,
     .         wp2, wp3, em, isValid)

! Latin hypercube sample generation 
! Generate p_height_time, an nnzp x nt x d array of random integers
              if ( LH_on ) then
                 i_rmd = mod(i1-1,sequence_length)
                 if ( i_rmd .eq. 0) then
                  call permute_height_time(gr%nnzp,nt,d+1,p_height_time)
                 endif
              endif
! End Latin hypercube generation

!              print*, 'hoc.F: i_rmd=', i_rmd

!             Call closure scheme and boundary conditions

              call sfc_var( upwp(1), vpwp(1), wpthlp(1), wprtp(1),
     .                      wp2(1), thlp2(1), rtp2(1), rtpthlp(1) )

              wpthlp(gr%nnzp)  = 0.
              wprtp(gr%nnzp)   = 0.
              wp2(gr%nnzp)     = 0.
              thlp2(gr%nnzp)   = 0.
              rtp2(gr%nnzp)    = 0.
              rtpthlp(gr%nnzp) = 0.

              wp3(1)           = 0.
              wp3(2)           = 0.
              wp3(gr%nnzp)     = 0.

              do k=2,gr%nnzp-1

                 call pdf_closure_new
     .                ( p(k), exner(k),
     .                  wmt(k), zm2zt(wp2,k), wp3(k), Sct(k),
     .                  rtm(k), zm2zt(rtp2,k), zm2zt(wprtp,k),
     .                  thlm(k), zm2zt(thlp2,k), zm2zt(wpthlp,k),
     .                  wp4(k),
     .                  wprtp2(k), wp2rtp(k),
     .                  wpthlp2(k), wp2thlp(k),
     .                  zm2zt(rtpthlp,k),                     
     .                  wprtpthlp(k),
     .                  cf(k), rcm(k),
     .                  wpthvp(k), wp2thvp(k), 
     .                  rtpthvp(k), thlpthvp(k),
     .                  wprcp(k), wp2rcp(k),
     .                  rtprcp(k), thlprcp(k),
     .                  pdf_parms(k,:),
     .                  crt1, crt2, cthl1, cthl2
     .                )

! Latin hypercube sampling.
                 if ( LH_on ) then

! Choose which rows of LH sample to feed into closure.
                    p_matrix(1:n,1:(d+1)) =
     .                 p_height_time(k,n*i_rmd+1:n*i_rmd+n,1:(d+1))

!                    print*, 'hoc.F: got past p_matrix'

! Generate LH sample, represented by X_u and X_nl, for level k
                    call lh_sampler(n, nt, d, p_matrix,
     .                       cf(k),
     .                       pdf_parms(k,:),
     .                       crt1, crt2, cthl1, cthl2,
     .                       rrm(k),
     .                       X_u, X_nl, sample_flag )

!                      print*, 'hoc.F: got past lh_sampler'

! Perform LH and analytic microphysical calculations
                    call micro_calcs(n, d, X_u, X_nl, sample_flag,
     .                       pdf_parms(k,:),
     .                       AKm_est(k), AKm(k), AKstd(k), AKstd_cld(k),
     .                       AKm_rcm(k), AKm_rcc(k), rcm_est(k) )

!                    print*, 'k, AKm_est=', k, AKm_est(k)
!                    print*, 'k, AKm=', k, AKm(k)

                endif ! End LH sampling part.

              enddo

!              print*, 'hoc.F: AKm=', AKm
!              print*, 'hoc.F: AKm_est=', AKm_est
              

!             Interpolate variables back to m grid
              wp4         = zt2zm( wp4 )
              wpthvp      = zt2zm( wpthvp )
              wpthvp(1)   = wpthlp( 1 )
              thlpthvp    = zt2zm( thlpthvp )
              thlpthvp(1) = thlp2( 1 ) 
              rtpthvp     = zt2zm( rtpthvp )   
              wprcp       = zt2zm( wprcp )       
              rtprcp      = zt2zm( rtprcp )     
              thlprcp     = zt2zm( thlprcp )   

 

!             Alternate approach for lower wp4 boundary condition:
!             call pdf scheme at the 1st momentum level, zm(1)

              call pdf_closure_new
     .             ( p(1), exner(1),
     .               wmm(1), wp2(1), zt2zm(wp3,1), Scm,
     .               zt2zm(rtm,1), rtp2(1), wprtp(1),
     .               zt2zm(thlm,1), thlp2(1), wpthlp(1),
     .               wp4(1),
     .               wprtp2(1), wp2rtp(1),
     .               wpthlp2(1), wp2thlp(1),
     .               rtpthlp(1),
     .               wprtpthlp(1),
     .               cf(1), rcm(1),
     .               wpthvp(1), wp2thvp(1), 
     .               rtpthvp(1), thlpthvp(1),
     .               wprcp(1), wp2rcp(1),
     .               rtprcp(1), thlprcp(1),
     .               pdf_parms(1,:),
     .               crt1, crt2, cthl1, cthl2
     .             )

              wp4(gr%nnzp)      = 0.
              wpthvp(gr%nnzp)   = 0.
              wp2thvp(gr%nnzp)  = 0.
              rtpthvp(gr%nnzp)  = 0.
              thlpthvp(gr%nnzp) = 0.                        

!             call Khairoutdinov and Kogan (2000) scheme for rain microphysics.

              if ( kk_rain ) then
                call rain( dt, thlm, p, exner, rhot, wmt, pdf_parms,
     .                     rcm, rtm, rrm, Nrm, cond_rrm, auto_rrm,
     .                     accr_rrm, cond_Nrm, auto_Nrm, sed_rrm,
     .                     diff_rrm, subs_rrm, Vr, mean_vol_rad,
! LH Kessler autoconversion.  Vince Larson 22 May 2005:
     .                     AKm_est, AKm,
     .                     runtype )

                rain_rate = ( rrm * Vr ) * ( 86400.0 * 1000.0 )  ! rate: mm/day

                F_precip = ( rrm * Vr ) * rhot * Lv              ! Precip. flux (W/m^2)

              endif

!             Bc for the second order moments

              select case ( runtype )
                case ( 0, 3)                           ! Normal and FIRE
                  call sfc_momentum_fluxes( um(2), vm(2),
     .                                      upwp(1), vpwp(1) )
                  if ( sfctype == 0 ) then
                    wpthlp(1) = SE
                    wprtp(1) = LE
                  else if ( sfctype == 1 ) then
                    call sfc_thermo_fluxes( um(2), vm(2),
     .                                      Tsfc, psfc, 
     .                                      thlm(2), rtm(2),
     .                                      wpthlp(1), wprtp(1) )
                  else
                     write(*,*) 'Invalid value of sfctype = ', sfctype
                     stop
                 endif

!               For GCSS-ARM case
                case ( 1 )
                  call arm_sfclyr( time, gr%zt(2), 1.1, 
     .                             thlm(2), um(2), vm(2),
     .                             upwp(1), vpwp(1), 
     .                             wpthlp(1), wprtp(1), ustar )

                case ( 2 ) 
                  call bomex_sfclyr( um(2), vm(2),
     .                               upwp(1), vpwp(1),
     .                               wpthlp(1), wprtp(1) )

                case ( 4 )
                  call wangara_sfclyr( time, um(2), vm(2),
     .                                 upwp(1), vpwp(1),
     .                                 wpthlp(1), wprtp(1) )

                case ( 6 )
                  call atex_sfclyr( um(2), vm(2), thlm(2), rtm(2),
     .                              upwp(1), vpwp(1),
     .                              wpthlp(1), wprtp(1) )

                case ( 7 ) 
                  call dycoms2_rf01_sfclyr( sfctype, Tsfc, psfc, 
     .                                      exner(1), um(2), vm(2), 
     .                                      thlm(2), rtm(2), 
     .                                      rhom(1), upwp(1), vpwp(1), 
     .                                      wpthlp(1), wprtp(1) )
                case ( 8 )
                  call astex_sfclyr( rhom(1), um(2), vm(2),
     .                               upwp(1), vpwp(1), wpthlp(1), 
     .                               wprtp(1) )
                case ( 9, 10 )
                  call dycoms2_rf02_sfclyr( um(2), vm(2), 
     .                                      thlm(2), rtm(2),
     .                                      psfc, upwp(1), vpwp(1), 
     .                                      wpthlp(1), wprtp(1) )

                case default
                 write(*,*) 'Invalid value of runtype = ', runtype
                 stop

              end select


!             Compute thvm
              thvm            = thlm + ep1 * T0 * rtm 
     .                          + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

!             Advance rtm/wprtp and thlm/wpthlp one time step
              call timestep_mixing
     .             ( dt, taum, wmm, wmt, wp2, wp3, Scm, exner, 
     .               rtm, wprtp, rtpthvp, rtm_forcing,     
     .               cond_rrm, auto_rrm, accr_rrm,             
     .               thlm, wpthlp, thlpthvp, thlm_forcing, isValid )

           enddo
#ifdef STATS
!           Store values of surface fluxes for statistics
            if (lstats_samp) then
              if ( ish > 0 ) then
                sfc%x(1,ish) = sfc%x(1,ish) + wpthlp(1) * rhom(1) * cp
                sfc%n(1,ish) = sfc%n(1,ish) + 1
              endif
              if ( ilh > 0 ) then
                sfc%x(1,ilh) = sfc%x(1,ilh) + wprtp(1)  * rhom(1) * Lv
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

!<<<<<<<<<<<<<<<<<<<<<<<<<<<< Closure loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!          Update winds

           umt
     .     = - wmt * ddzm( zt2zm( um ) )
     .       + fcor * ( vm - vg )
           call update_wind( um, umt, upwp(1), Khm, dtmain, upwp, 
     .                       isValid)
           um(1)       = 2. * um(2) - um(3)
           um(gr%nnzp) = 2. * um(gr%nnzp-1) - um(gr%nnzp-2 )

!          Check for a degenerate state
           if ( .not. isValid ) then
              call deallocate_model_arrays( )
!              deallocate( fzt%var, fzt%z, fzm%var, fzm%z, 
!     .                    fpdf%var, fpdf%z, fsfc%var, fsfc%z )
              return
           endif

           vmt 
     .     = - wmt * ddzm( zt2zm( vm ) )
     .       - fcor * ( um - ug )
           call update_wind( vm, vmt, vpwp(1), Khm,dtmain,vpwp,isValid)
           vm(1)       = 2. * vm(2) - vm(3)
           vm(gr%nnzp) = 2. * vm(gr%nnzp-1) - vm(gr%nnzp-2 )
           
           if ( .not. isValid ) then
              call deallocate_model_arrays( )
!              deallocate( fzt%var, fzt%z, fzm%var, fzm%z, 
!     .                    fpdf%var, fpdf%z, fsfc%var, fsfc%z )
              return
           endif

           ! Shear production.  Brian
           DO k = 1, gr%nnzp-1, 1
              shear(k) = -upwp(k) * ( um(k+1) - um(k) ) * gr%dzm(k)
     .                  - vpwp(k) * ( vm(k+1) - vm(k) ) * gr%dzm(k)
           ENDDO
           shear(gr%nnzp) = 0.0

           ! Store the saturation mixing ratio for output purposes.  Brian
           DO k = 1, gr%nnzp, 1
              rsm(k) =
     .         rsat(p(k), (thlm(k)+(Lv/(Cp*exner(k)))*rcm(k))*exner(k))
           ENDDO

! Added to allow tuning without using the zm file
#ifdef STATS
       wp2zt     = zm2zt( wp2 )
       thlp2zt   = zm2zt( thlp2 )
       wpthlpzt  = zm2zt( wpthlp )
       wprtpzt   = zm2zt( wprtp )
       rtp2zt    = zm2zt( rtp2 )
#endif STATS


!          Altered by Brian

#ifdef STATS
           call stats_accumulate( time, dtmain )
!           call stats_end_timestep( time, dtmain )
!           print*, 'hoc.F: got past stats_accumulate'
           call stats_end_timestep( time, dtmain, dtclosure )   ! Brian
!           print*, 'hoc.F: got past stats_end_timestep'
#endif STATS


        enddo

!<<<<<<<<<<<<<<<<<<<<<<<<<< Main time stepping loop <<<<<<<<<<<<<<<<<<<<<<<<

! Free memory

        call deallocate_model_arrays( )

! Needed when hoc became a called subroutine rather than a program
!        deallocate( fzt%var, fzt%z, fzm%var, fzm%z, 
!     .              fpdf%var, fpdf%z, fsfc%var, fsfc%z )

        contains

!       This was taken out of the main time stepping loop for readability

        logical function Test_Variables(  ) 

        implicit none

!       Internal Variables
        logical             :: isNaN

        Test_Variables = .true.
           
!       Check whether any variable array contains a NaN for
!       um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp, 
!       wp2, & wp3.       
        call TestNaN(um, isNaN, nzmax)   
        if ( isNaN ) then
!             print *, "NaN! "
!             print *, "um= ", um
          Test_Variables = .false.
        endif
        call TestNaN(vm, isNaN, nzmax)    
        if ( isNaN ) then
!          print *, "NaN! " 
!          print *, "vm= ", vm
          Test_Variables = .false.
        endif
        call TestNaN(thlm, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN! " 
!          print *, "thlm= ", thlm
          Test_Variables = .false.
        endif
        call TestNaN(rtm, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN! " 
!          print *, "rtm= ", rtm
          Test_Variables = .false.
        endif
        call TestNaN(rtp2, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN! " 
!          print *, "rtp2= ", rtp2
           Test_Variables = .false.
        endif
        call TestNaN(thlp2, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN! " 
!          print *, "thlp2= ", thlp2
          Test_Variables = .false.
        endif
        call TestNaN(wprtp, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN!" 
!          print *, "wprtp= ", wprtp
          Test_Variables = .false.
        endif
        call TestNaN(wpthlp, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN!"
!          print *, "wpthlp= ", wpthlp
          Test_Variables = .false.
        endif
        call TestNaN(rtpthlp, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN!" 
!          print *, "rtpthlp= ", rtpthlp
          Test_Variables = .false.
        endif
        call TestNaN(wp2, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN!" 
!          print *, "wp2= ", wp2
          Test_Variables = .false.
        endif
        call TestNaN(wp3, isNaN, nzmax)   
        if ( isNaN ) then
!          print *, "NaN!" 
!          print *, "wp3= ", wp3
          Test_Variables = .false.
        endif
      
        if ( .not. Test_Variables) print *, "NaN in an array" 

        return

        end function Test_Variables
!------------------------------------------------------------------------

        end subroutine hoc_sub

!------------------------------------------------------------------------
! subroutine to check if a real is NaN, +inf or -inf
! This would be better as a logical function, but I had 
! technical issues doing that. 
! Paremeter NaN does not seem to work w/ '.eq.'
! Got a tip from Andy Vaught to use transfer( ) for portability
! Ideally this would be done with intrinsic module ieee in fortran 2003
!------------------------------------------------------------------------

        subroutine TestNaN( tVar, isNaN, tDim )

        implicit none
      
!----- Arguments ------------------------------------------------------
        integer, intent(in)               :: tDim
        real, intent(in), dimension(tDim) :: tVar

        logical, intent(out)              :: isNaN

!---- Internal Variables ----------------------------------------------

        integer                           :: i          ! loop variable
        integer                           :: bits_pinf, bits_ninf 

        real                              :: PosInf, NegInf

        bits_pinf = Z"7F800000"
        bits_ninf = Z"FF800000"

        PosInf    = transfer(bits_pinf, 1.0)
        NegInf    = transfer(bits_ninf, 1.0)

        isNaN     = .false.

        do i = 1, tDim
          if (tVar(i) .eq. PosInf) then
             isNaN = .true. 
             exit
          endif
          if (tVar(i) .eq. NegInf) then
             isNaN = .true. 
             exit
          endif
        enddo
        end subroutine TestNaN
        
        end module hoc
!------------------------------------------------------------------------
