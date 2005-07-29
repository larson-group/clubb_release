!-----------------------------------------------------------------------
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
     .                     input_rtpthlp, input_upwp, input_vpwp

        public grads_fields_reader
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
     .                                  upwp, vpwp )
        use grid_class
        use inputfile_class

        implicit none

!       Arguments
        integer, intent(in)      :: timestep

        real, optional, intent(out), 
     .        dimension(gr%nnzp) :: um, vm, rtm, thlm, wp2,
     .                              wprtp, wpthlp, wp3, rtp2,
     .                              thlp2, rtpthlp, upwp, vpwp

!       Internal
        logical error
        type (inputgrads) fread_var
        real, dimension(gr%nnzp + 1) :: tmp1

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
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            rtp2(1) = 2. * rtp2(2) - rtp2(3)
          endif

          if (present( thlp2 ) .and. input_thlp2) then 
            call get_var( fread_var, "thlp2", 
     .                    timestep, tmp1(1:gr%nnzp), error )
            call interpolate_m( tmp1(1:gr%nnzp), thlp2(2:gr%nnzp) )
            ! Use linear extrapolation, assuming equal grid spacing,
            ! although the boundary flux is usually solved in LES
            ! or hoc via a subroutine like sfc_var.
            thlp2(1) = 2. * thlp2(2) - thlp2(3)
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
            wp2(1:gr%nnzp) = tmp1(1:gr%nnzp)
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

      end module inputfields
!-----------------------------------------------------------------------
