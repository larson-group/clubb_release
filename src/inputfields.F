!-----------------------------------------------------------------------
! $Id: inputfields.F,v 1.8 2008-02-22 19:30:11 faschinj Exp $

! Module inputfields

!  This exists because I wanted to keep the grads_reader code 
!  generalized and bypass having to pass the datafile as a parameter,
!  since there may be situations where the fields are be calculated 
!  analytically or by calling a different model without reading a datafile.
!  Using a module also saves the trouble of writing an interface definition
!  within the hoc_inputfields code.
!-----------------------------------------------------------------------
      module inputfields
        implicit none

!----- Run information--------------------------------------------------
        character(len=80), public :: datafile
        character(3), public      :: input_type

        logical, public :: input_um, input_vm, input_rtm, input_thlm,
     .                     input_wp2, input_wprtp, input_wpthlp, 
     .                     input_wp3, input_rtp2, input_thlp2, 
     .                     input_rtpthlp, input_upwp, input_vpwp,
     .                     input_ug, input_vg, input_rcm, 
     .                     input_wmt, input_exner, input_em,
     .                     input_p, input_rhot, input_rhom,
     .                     input_Lscale, input_lup, input_ldown,
     .                     input_Kht, input_Khm, input_taum, input_taut,
     .                     input_thvm, input_rrm,input_Nrm, 
     .                     input_rsnowm, input_ricem, input_rgraupelm, 
     .                     input_thlm_forcing, input_rtm_forcing,
     .                     input_rtm_mc, input_thlm_mc, 
     .                     input_up2, input_vp2, input_Scm, input_Ncm, 
     .                     input_Ncnm, input_Nim, input_cf, input_Sct 


        public grads_fields_reader, compute_timestep
        private interpolate_m

!-----------------------------------------------------------------------
! SUBROUTINE grads_fields_reader( ) 

! Thanks to the interface feature in Fortran 90, these could really be
! be calculated from anything.  This subroutine obtains them by reading
! GrADS data from a specified file.

! Calls:
! subroutine open_grads_read
! subroutine get_var
! subroutine close_grads_read

!-----------------------------------------------------------------------
        contains
        subroutine grads_fields_reader( timestep,
     .                                  um, vm, rtm, thlm,
     .                                  wp2, wprtp, wpthlp, wp3,
     .                                  rtp2, thlp2, rtpthlp, 
     .                                  upwp, vpwp, ug, vg, rcm, 
     .                                  wmt, exner, em, p, rhot, 
     .                                  rhom, Lscale, lup, ldown,
     .                                  Kht, Khm, taum, taut, thvm,
     .                                  rrm, Nrm, rsnowm, ricem, 
     .                                  rgraupelm, thlm_forcing, 
     .                                  rtm_forcing, rtm_mc, thlm_mc, 
     .                                  up2, vp2, Scm, Ncm, Ncnm, Nim,
     .                                  cf, Sct )
        use grid_class
        use inputfile_class

        implicit none

!       Arguments
        integer, intent(in)      :: timestep

        real, optional, intent(out), 
     .        dimension(:) :: um, vm, rtm, thlm, wp2,
     .                        wprtp, wpthlp, wp3, rtp2,
     .                        thlp2, rtpthlp, upwp, vpwp,
     .                        ug, vg, rcm, 
     .                        wmt, exner, em, p, rhot, 
     .                        rhom, Lscale, lup, ldown,
     .                        Kht, Khm, taum, taut, thvm,
     .                        rrm, Nrm, rsnowm, ricem, rgraupelm, 
     .                        thlm_forcing, rtm_forcing,
     .                        rtm_mc, thlm_mc, 
     .                        up2, vp2, Scm, Ncm, Ncnm, Nim,
     .                        cf, Sct  

!       Internal
        logical error
        type (inputgrads) fread_var
        real, dimension(gr%nnzp + 1) :: tmp1
        integer k

        select case( input_type )
        case( "hoc" )
          call open_grads_read( 15, trim(datafile)//"_zt.ctl", 
     .                          fread_var )

!         zt file 

          if (present( um ) .and. input_um) then
            call get_var( fread_var, "um", timestep,
     .                    um(1:gr%nnzp), error )
          endif
      
          if (present( vm ) .and. input_vm) then 
            call get_var( fread_var, "vm", timestep,
     .                    vm(1:gr%nnzp), error )
          endif

          if (present( rtm ) .and. input_rtm) then 
            call get_var( fread_var, "rtm", timestep,
     .                    rtm(1:gr%nnzp), error )
          endif

          if (present( thlm ) .and. input_thlm) then 
            call get_var( fread_var, "thlm", 
     .                    timestep,
     .                    thlm(1:gr%nnzp), error )
          endif

          if (present( wp3 ) .and. input_wp3) then 
            call get_var( fread_var, "wp3", timestep,
     .                    wp3(1:gr%nnzp), error )
          endif
          if (present( taut ) .and. input_taut) then 
            call get_var( fread_var, "taut", timestep,
     .                    taut(1:gr%nnzp), error )
          endif
          if (present( rrm ) .and. input_rrm) then 
            call get_var( fread_var, "rrm", timestep,
     .                    rrm(1:gr%nnzp), error )
          endif
          if (present( rsnowm ) .and. input_rsnowm) then 
            call get_var( fread_var, "rsnowm", timestep,
     .                    rsnowm(1:gr%nnzp), error )
          endif
          if (present( ricem ) .and. input_ricem) then 
            call get_var( fread_var, "ricem", timestep,
     .                    ricem(1:gr%nnzp), error )
          endif
          if (present( rgraupelm ) .and. input_rgraupelm) then 
            call get_var( fread_var, "rgraupelm", timestep,
     .                    rgraupelm(1:gr%nnzp), error )
          endif

!--------------------------------------------------------
! Added variables for hoc_restart
          if (present( p ) .and. input_p) then
            call get_var( fread_var, "p", timestep,
     .                    p(1:gr%nnzp), error )
          endif
          if (present( exner ) .and. input_exner) then
            call get_var( fread_var , "exner", timestep,
     .                    exner(1:gr%nnzp), error)
          endif
          if (present( ug ) .and. input_ug) then
            call get_var( fread_var , "ug", timestep,
     .                    ug(1:gr%nnzp), error)
          endif
          if (present( vg ) .and. input_vg) then
            call get_var( fread_var , "vg", timestep,
     .                    vg(1:gr%nnzp), error)
          endif
          if (present( rcm ) .and. input_rcm) then
            call get_var( fread_var , "rcm", timestep,
     .                    rcm(1:gr%nnzp), error)
          endif
          if (present( wmt ) .and. input_wmt) then
            call get_var( fread_var , "wm", timestep,
     .                    wmt(1:gr%nnzp), error)
          endif
          if (present( rhot ) .and. input_rhot) then
            call get_var( fread_var , "rhot", timestep,
     .                    rhot(1:gr%nnzp),error)
          endif
          if (present( Lscale ) .and. input_Lscale) then
            call get_var( fread_var , "lscale", timestep,
     .                    Lscale(1:gr%nnzp),error)
          endif
          if (present( lup ) .and. input_lup) then
            call get_var( fread_var , "lup", timestep,
     .                    lup(1:gr%nnzp),error)
          endif
          if (present( ldown ) .and. input_ldown) then
            call get_var( fread_var , "ldown", timestep,
     .                    ldown(1:gr%nnzp),error)
          endif
          if (present( Kht ) .and. input_Kht) then
            call get_var( fread_var , "kht", timestep,
     .                    Kht(1:gr%nnzp),error)
          endif
          if (present( thvm ) .and. input_thvm) then
            call get_var( fread_var , "thvm", timestep,
     .                    thvm(1:gr%nnzp),error)
          endif
          if (present( thlm_forcing ) .and. input_thlm_forcing ) then
            call get_var( fread_var , "thlm_f", timestep,
     .                    thlm_forcing(1:gr%nnzp),error)
          endif
          if (present( rtm_forcing ) .and. input_rtm_forcing ) then
            call get_var( fread_var , "rtm_f", timestep,
     .                    rtm_forcing(1:gr%nnzp),error)
          endif
          if (present( rtm_mc ) .and. input_rtm_mc ) then
            call get_var( fread_var , "rtm_mc", timestep,
     .                    rtm_mc(1:gr%nnzp),error)
          endif
          if (present( thlm_mc ) .and. input_thlm_mc ) then
            call get_var( fread_var , "thlm_mc", timestep,
     .                    thlm_mc(1:gr%nnzp),error)
          endif
          if (present( Ncm ) .and. input_Ncm) then
            call get_var( fread_var , "Ncm", timestep,
     .                    Ncm(1:gr%nnzp),error)
          endif
          if (present( Ncnm ) .and. input_Ncnm) then
            call get_var( fread_var , "Ncnm", timestep,
     .                    Ncnm(1:gr%nnzp),error)
          endif
          if (present( Nim ) .and. input_Nim) then
            call get_var( fread_var , "Nim", timestep,
     .                    Nim(1:gr%nnzp),error)
          endif
          if (present( cf ) .and. input_cf) then
            call get_var( fread_var , "cf", timestep,
     .                    cf(1:gr%nnzp),error)
          endif
          if (present( Nrm ) .and. input_Nrm ) then
            call get_var( fread_var , "Nrm", timestep,
     .                    Nrm(1:gr%nnzp),error)
          endif
          if (present( Sct ) .and. input_Sct ) then
            call get_var( fread_var , "sc", timestep,
     .                    Sct(1:gr%nnzp),error)
          endif

!--------------------------------------------------------
          call close_grads_read( fread_var )

!         zm file

          call open_grads_read( 15, trim(datafile)//"_zm.ctl",
     .                          fread_var )

          if (present( wp2 ) .and. input_wp2) then 
            call get_var( fread_var, "wp2", timestep,
     .                    wp2(1:gr%nnzp), error )
          endif

          if (present( wprtp ) .and. input_wprtp) then 
            call get_var( fread_var, "wprtp", 
     .                    timestep, wprtp(1:gr%nnzp),
     .                    error )
          endif

          if (present( wpthlp ) .and. input_wpthlp) then 
            call get_var( fread_var, "wpthlp", 
     .                    timestep, 
     .                    wpthlp(1:gr%nnzp), 
     .                    error )
          endif

          if (present( rtp2 ) .and. input_rtp2) then 
             call get_var( fread_var, "rtp2", 
     .                     timestep,
     .                     rtp2(1:gr%nnzp), error )
          endif

          if (present( thlp2 ) .and. input_thlp2) then 
             call get_var( fread_var, "thlp2", 
     .                     timestep,
     .                     thlp2(1:gr%nnzp), error )
          endif

          if (present( rtpthlp ) .and. input_rtpthlp) then 
             call get_var( fread_var, "rtpthlp", 
     .                     timestep, 
     .                     rtpthlp(1:gr%nnzp),
     .                     error )
          endif

          if (present( upwp ) .and. input_upwp) then 
             call get_var( fread_var, "upwp", 
     .                     timestep,
     .                     upwp(1:gr%nnzp), error )
          endif

          if (present( vpwp ) .and. input_vpwp) then 
             call get_var( fread_var, "vpwp", 
     .                     timestep,
     .                     vpwp(1:gr%nnzp), error )
          endif
!-----------------------------------------------------------
         if (present( em ) .and. input_em) then
            call get_var( fread_var, "em",
     .                    timestep,
     .                    em(1:gr%nnzp), error )
         endif
         if (present( rhom ) .and. input_rhom) then
            call get_var( fread_var, "rhom",
     .                    timestep,
     .                    rhom(1:gr%nnzp), error )
         endif
         if (present( Khm ) .and. input_Khm) then
            call get_var( fread_var, "khm",
     .                    timestep,
     .                    Khm(1:gr%nnzp), error )
         endif
         if (present( taum ) .and. input_taum) then
            call get_var( fread_var, "taum",
     .                    timestep,
     .                    taum(1:gr%nnzp), error )
         endif
         if (present( up2 ) .and. input_up2) then
            call get_var( fread_var, "up2",
     .                    timestep,
     .                    up2(1:gr%nnzp), error )
         endif
         if (present( vp2 ) .and. input_vp2) then
            call get_var( fread_var, "vp2",
     .                    timestep,
     .                    vp2(1:gr%nnzp), error )
         endif
         if (present( Scm ) .and. input_Scm ) then
            call get_var( fread_var, "scm",
     .                    timestep,
     .                    Scm(1:gr%nnzp), error )
         endif






!-----------------------------------------------------------


          if ( error ) stop "oops, get_var failed in field_reader"

          call close_grads_read( fread_var )

        case( "rf1" )   ! special case for DYCOMS II RF01
          call open_grads_read( 15, trim(datafile)//"_coamps_sm.ctl", 
     .                          fread_var )
!         stats_sm

          if (present( um ) .and. input_um) then
            call get_var( fread_var, "um", timestep,
     .                    tmp1(1:gr%nnzp-2), error )
            um(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
            um(2) = 2. * um(3) - um(4) ! Double extrapolation
            um(1) = 2. * um(2) - um(3)
          endif

          if (present( vm ) .and. input_vm) then 
            call get_var( fread_var, "vm", timestep,
     .                    tmp1(1:gr%nnzp-2), error )
            vm(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
            vm(2) = 2. * vm(3) - vm(4) ! Double extrapolation
            vm(1) = 2. * vm(2) - vm(3)
          endif

          if (present( rtm ) .and. input_rtm) then 
            call get_var( fread_var, "qtm", timestep,
     .                    tmp1(1:gr%nnzp-2), error )
            rtm(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
            rtm(1:2) = rtm(3)
          endif

          if (present( thlm ) .and. input_thlm) then
            call get_var( fread_var, "thlm", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp-2), error )
            thlm(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
            thlm(1:2) = thlm(3)
          endif

          if (present( wp3 ) .and. input_wp3) then 
            call get_var( fread_var, "wp3", timestep,
     .                    tmp1(1:gr%nnzp-2), error )
            wp3(3:gr%nnzp) = tmp1(1:gr%nnzp-2) 
            wp3(1:2) = 0.
          endif

          if (present( wprtp ) .and. input_wprtp) then
            call get_var( fread_var, "wpqtp", 
     .                    timestep, tmp1(1:gr%nnzp-1),
     .                    error )
            call interpolate_m( tmp1(1:gr%nnzp-1), wprtp(3:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            wprtp(2) = 2. * wprtp(3) - wprtp(4) 
            wprtp(1) = 2. * wprtp(2) - wprtp(3) 
          endif

          if (present( wpthlp ) .and. input_wpthlp) then 
            call get_var( fread_var, "wpthlp", 
     .                    timestep, tmp1(1:gr%nnzp-1), 
     .                    error )
            call interpolate_m( tmp1(1:gr%nnzp-1), wpthlp(3:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            wpthlp(2) = 2. * wpthlp(3) - wpthlp(4)
            wpthlp(1) = 2. * wpthlp(2) - wpthlp(3)
          endif

          if (present( rtp2 ) .and. input_rtp2) then 
            call get_var( fread_var, "qtp2", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp-1), error )
            call interpolate_m( tmp1(2:gr%nnzp), rtp2(3:gr%nnzp) )
            ! Using a linear interpolation here resulted in negatives.
            rtp2(1:2) =  rtp2(3)
            if ( any (rtp2(1:gr%nnzp) < 0.0 ) ) then
! %% debug
!              print *, "Some values of rtp2 are negative, compensating."
! %% debug
              do k=1, gr%nnzp
                rtp2(k) = max(rtp2(k), 0.0)
              enddo
            endif
          endif

          if (present( thlp2 ) .and. input_thlp2) then 
            call get_var( fread_var, "thlp2", 
     .                    timestep, tmp1(1:gr%nnzp), error )
            ! Using a linear interpolation here resulted in negatives.
            thlp2(1:2) = thlp2(3)
          endif

          if (present( rtpthlp ) .and. input_rtpthlp) then 
            call get_var( fread_var, "qtpthlp", 
     .                    timestep, tmp1(1:gr%nnzp-1),
     .                    error )
            call interpolate_m( tmp1(1:gr%nnzp-1), rtpthlp(3:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            rtpthlp(2) = 2. * rtpthlp(3) - rtpthlp(4)
            rtpthlp(1) = 2. * rtpthlp(2) - rtpthlp(3)
          endif

          if ( error ) stop "oops, get_var failed in field_reader"

          call close_grads_read( fread_var )

        case( "les" )
          call open_grads_read( 15, trim(datafile)//"_coamps_sm.ctl", 
     .                          fread_var )

!         stats_sm
          if (present( um ) .and. input_um) then
            call get_var( fread_var, "um", timestep,
     .                    tmp1(1:gr%nnzp), error )
            um(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
            um(1) = 2. * um(2) - um(3) ! Extrapolated as in hoc.F
          endif

          if (present( vm ) .and. input_vm) then 
            call get_var( fread_var, "vm", timestep,
     .                    tmp1(1:gr%nnzp), error )
            vm(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
            vm(1) = 2. * vm(2) - vm(3)  ! Extrapolated as in hoc.F
          endif

          if (present( rtm ) .and. input_rtm) then 
            call get_var( fread_var, "qtm", timestep,
     .                    tmp1(1:gr%nnzp), error )
            rtm(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
            rtm(1) = rtm(2)  ! Computed as in mixing.F
          endif

          if (present( thlm ) .and. input_thlm) then
            call get_var( fread_var, "thlm", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp), error )
            thlm(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
            thlm(1) = thlm(2) ! Computed as in mixing.F
          endif

          if (present( wp3 ) .and. input_wp3) then 
            call get_var( fread_var, "wp3", timestep,
     .                    tmp1(1:gr%nnzp), error )
            wp3(2:gr%nnzp) = tmp1(1:gr%nnzp-1) 
            wp3(1) = 0.  ! Computed as in hoc.F
          endif

!          if (present( wp2 )) then 
!            call get_var( fread_var, "wp2", timestep,
!     .                    tmp1(1:gr%nnzp), error )
!            call interpolate_m( tmp1(1:gr%nnzp), wp2(2:gr%nnzp) )
!            ! Use linear extrapolation, assuming equal grid spacing,
!            ! although the boundary flux is usually solved in LES
!            ! or hoc via a subroutine like sfc_var.
!            wp2(1) = 2. * wp2(2) - wp2(3)
!          endif

          if (present( wprtp ) .and. input_wprtp) then
            call get_var( fread_var, "wpqtp", 
     .                    timestep, tmp1(1:gr%nnzp),
     .                    error )
            call interpolate_m( tmp1(1:gr%nnzp), wprtp(2:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            wprtp(1) = 2. * wprtp(2) - wprtp(3) 
          endif

          if (present( wpthlp ) .and. input_wpthlp) then 
            call get_var( fread_var, "wpthlp", 
     .                    timestep, tmp1(1:gr%nnzp), 
     .                    error )
            call interpolate_m( tmp1(1:gr%nnzp), wpthlp(2:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            wpthlp(1) = 2. * wpthlp(2) - wpthlp(3)
          endif

          if (present( rtp2 ) .and. input_rtp2) then 
            call get_var( fread_var, "qtp2", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp), error )
            call interpolate_m( tmp1(1:gr%nnzp), rtp2(2:gr%nnzp) )
            ! Using a linear interpolation here resulted in negatives.
            rtp2(1) =  rtp2(2)
            if ( any (rtp2(1:gr%nnzp) < 0.0 ) ) then
! %% debug
!              print *, "Some values of rtp2 are negative, compensating."
! %% debug
              do k=1, gr%nnzp
                rtp2(k) = max(rtp2(k), 0.0)
              enddo
            endif
          endif

          if (present( thlp2 ) .and. input_thlp2) then 
            call get_var( fread_var, "thlp2", 
     .                    timestep, tmp1(1:gr%nnzp), error )
            ! Using a linear interpolation here resulted in negatives.
            thlp2(1) = thlp2(2)
          endif

          if (present( rtpthlp ) .and. input_rtpthlp) then 
            call get_var( fread_var, "qtpthlp", 
     .                    timestep, 
     .                    tmp1(1:gr%nnzp),
     .                    error )
            call interpolate_m( tmp1(1:gr%nnzp), rtpthlp(2:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            rtpthlp(1) = 2. * rtpthlp(2) - rtpthlp(3)
          endif

          if ( error ) stop "oops, get_var failed in field_reader"

          call close_grads_read( fread_var )
        end select


        select case( input_type )
        case( "les", "rf1" )
!         stats_sw

          call open_grads_read( 15, trim(datafile)//"_coamps_sw.ctl", 
     .                          fread_var )
         ! no interpolation is required, however, the stats_sw files have
         ! an extra top z-level, and wpup_sgs must be added to make the
         ! u'w' and v'w' terms as they are in Hoc.

          if (present( upwp ) .and. input_upwp) then 
            call get_var( fread_var, "wpup", 
     .                    timestep, tmp1(1:gr%nnzp+1), error )
            upwp(1:gr%nnzp) = tmp1(1:gr%nnzp) 

            call get_var( fread_var, "wpup_sgs", 
     .                    timestep, tmp1(1:gr%nnzp+1), error )
            upwp(1:gr%nnzp) = tmp1(1:gr%nnzp) + upwp(1:gr%nnzp)
          endif
          if ( error ) stop "get_var failed for upwp in field_reader"

          if (present( vpwp ) .and. input_vpwp) then
            call get_var( fread_var, "wpvp", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp+1), error )
            vpwp(1:gr%nnzp) = tmp1(1:gr%nnzp) 
            call get_var( fread_var, "wpvp_sgs", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp+1), error )
            vpwp(1:gr%nnzp) = tmp1(1:gr%nnzp) + vpwp(1:gr%nnzp)
          endif
          if ( error ) stop "get_var failed for vpwp in field_reader"

          if (present( wp2 ) .and. input_wp2) then 
            call get_var( fread_var, "wp2", 
     .                    timestep,
     .                    tmp1(1:gr%nnzp+1), error )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            wp2(1:gr%nnzp) = tmp1(1:gr%nnzp)
            if ( any (wp2(1:gr%nnzp) < 0.0 ) ) then
! %% debug
!              print *, "Some values of wp2 are negative, compensating."
! %% debug
              do k=1, gr%nnzp
                wp2(k) = max(wp2(k), 0.0)
              enddo
            endif
          endif
          if ( error ) stop "get_var failed for wp2 in field_reader"

          call close_grads_read( fread_var )
        end select

        return
        end subroutine grads_fields_reader
!-----------------------------------------------------------------------
! Works like zt2zm with a variation to take into account the higher LES
!-----------------------------------------------------------------------
        subroutine interpolate_m( x_in_t, x_in_m )

        implicit none
        
        real, intent(in), dimension(:)  :: x_in_t ! nnzp
        real, intent(out), dimension(:) :: x_in_m ! nnzp-1

        integer k

        do k=1, size(x_in_t)-1, 1
          x_in_m(k) = ( x_in_t(k) + x_in_t(k+1) ) / 2
        enddo

        return
        end subroutine interpolate_m
!-------------------------------------------------------------------------
        subroutine compute_timestep( current_time, nearest_timestep )
!
!       Description: Determines the nearest timestep stored in the GRADS
!       file for the current_time in seconds recieved.
!
!-------------------------------------------------------------------------
          use inputfile_class
          implicit none
          
          ! Input Variable(s)
          real, intent(in) :: 
     .       current_time ! Time since start day
          
          ! Output Variable(s)
          integer, intent(out) :: 
     .       nearest_timestep ! Nearest timestep to current_time
          
          ! Local Variables
          type (inputgrads) fread_var         
          real delta_time
          
          call open_grads_read( 15, trim(datafile)//"_zt.ctl", 
     .                          fread_var )

          ! Difference between GrADS file start time and the restart
          ! time
          delta_time = current_time - fread_var%time
          
          ! Reporting
          print *, "Time Initial: ", fread_var%time
          print *, "Restart Time: ", current_time
          print *, "Delta Time: ", delta_time
          print *, "dt: ", fread_var%dtwrite
      
          ! Determines the closest recorded timestep to the restart
          ! time.
          nearest_timestep = nint(delta_time / 60.0)
          
          print *,"Nearest Timestep: ", nearest_timestep
          
          call close_grads_read( fread_var )
          
        end subroutine compute_timestep
        
      end module inputfields
!-----------------------------------------------------------------------
      
