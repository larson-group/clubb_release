!------------------------------------------------------------------------
! $Id: sounding.F,v 1.4 2005-06-23 20:07:53 dschanen Exp $

! subroutine read_sounding( )

! Subroutine to initialize model variables from namelist file
! Notes: for whatever reason, this subroutine is very flakey, modify it
!        with caution.
!------------------------------------------------------------------------
        subroutine read_sounding( thlm, rtm, um, vm, ugm, vgm, runfile, 
     .                            runtype )
!     .                            sclrm ) ! optional

        use grid_class
        use constants, only: sclrm_dimension

        implicit none

! External

        real, external :: linint

! Input

        character(len=*), intent(in) :: runfile
        integer, intent(in)          :: runtype

! Parameter

        integer, parameter :: nmaxsnd = 200

! Output

        real, intent(out), dimension(gr%nnzp) :: thlm
        real, intent(out), dimension(gr%nnzp) :: rtm
        real, intent(out), dimension(gr%nnzp) :: um
        real, intent(out), dimension(gr%nnzp) :: vm
        real, intent(out), dimension(gr%nnzp) :: ugm
        real, intent(out), dimension(gr%nnzp) :: vgm

!        real, optional, intent(out), 
!     .        dimension(sclrm_dimension, gr%nnzp) :: sclrm

! Input variables from namelist

        integer nlevels
        real z(nmaxsnd),theta(nmaxsnd),rt(nmaxsnd)
        real u(nmaxsnd),v(nmaxsnd),ug(nmaxsnd),vg(nmaxsnd)
        real, dimension(nmaxsnd) :: sclr

        namelist /sounding/ nlevels, z, theta, rt, u, v, ug, vg
        namelist /scalers/ sclr    ! unimplemented
    
! Internal

        integer i,k

! Read sounding namelist
        open(10, file = runfile, status = 'old' )
        read(10, nml = sounding )
        close(10)
        if ( nlevels > nmaxsnd ) then
           write(*,*) 'Error in sounding: nlevels > nmaxsnd'
           write(*,*) 'nlevels = ',nlevels
           write(*,*) 'nmaxsnd = ',nmaxsnd
           stop 'STOP in sounding'
        endif

! Error check: if lowest themodynamic grid height is lower than the
! lowest value from the input sounding, then the linear interpolation
! scheme will fail

        if ( gr%zt(2) < z(1) ) then
           write(*,*) 'First level of input sounding must be',
     .     ' below first thermodynamic level'
           write(*,*) ' first sounding level z(1) = ',z(1)
           write(*,*) ' first thermodynamic level gr%zt(2) = ',gr%zt(2)
           stop 'STOP in sounding'
        endif

! First sounding level should be near ground value

        if ( abs(z(1)) > 1.e-8 ) then
           write(*,*) 'First level of input sounding must be z=0'
           stop 'STOP in sounding'
        else
            um(1)   = u(1)
            vm(1)   = v(1)
            ugm(1)  = ug(1)
            vgm(1)  = vg(1)
            thlm(1) = theta(1)
            rtm(1)  = rt(1)
        endif

! Use linear interpolation from two nearest prescribed grid points
! (one above and one below) to initialize mean quantities in the model
! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

        do i=2, gr%nnzp
          k=1
!  40      if ( z(k) >= gr%zt(i) ) then
          do while ( z(k) < gr%zt(i) )
            k=k+1
            if ( k > nlevels ) then
              write(*,*) 'STOP Not enough sounding data to ',
     .                   'initialize grid:'
              write(*,'(a,f7.1,/a,f7.1)') '  highest sounding level'
     .        ,z(nlevels)
     .        ,'  should be higher than highest thermodynamic point'
     .        ,gr%zt(gr%nnzp)

              stop 'STOP in sounding'
            endif  ! k > nlevels
          enddo
            IF ( runtype /= 9 ) THEN  ! Regular situation w/ linear int.

            um(i)   = linint( gr%zt(i), z(k), z(k-1), u(k), u(k-1) )
            vm(i)   = linint( gr%zt(i), z(k), z(k-1), v(k), v(k-1) )
            ugm(i)  = linint( gr%zt(i), z(k), z(k-1), ug(k), ug(k-1) )
            vgm(i)  = linint( gr%zt(i), z(k), z(k-1), vg(k), vg(k-1) )
            thlm(i) = linint( gr%zt(i), z(k), z(k-1), 
     .                        theta(k), theta(k-1) )
            rtm(i)  = linint( gr%zt(i), z(k), z(k-1), rt(k), rt(k-1) )

            ELSE  ! DYCOMS II RF02 case

               IF ( gr%zt(i) < 795.0 ) THEN
                  um(i)   =  6.5 + (2.0*gr%zt(i))/1000.0
                  vm(i)   = -8.5 + (5.0*gr%zt(i))/1000.0
                  ugm(i)  = um(i)
                  vgm(i)  = vm(i)
                  thlm(i) = 288.3
                  rtm(i)  = (9.45)/1000.0
               ELSE
                  um(i)   =  6.5 + (2.0*gr%zt(i))/1000.0
                  vm(i)   = -8.5 + (5.0*gr%zt(i))/1000.0
                  ugm(i)  = um(i)
                  vgm(i)  = vm(i)
                  thlm(i) = 295.0 + ( (gr%zt(i) - 795.0)**(1.0/3.0) )
                  rtm(i)  = (  5.0 - 3.0 
     .            * ( 1.0 - EXP( (795.0 - gr%zt(i))/500.0 ) )  )/1000.0
               ENDIF

            ENDIF ! runtype /= 9
!          else
!            goto 40
!          endif      !  z(k) >= gr%zt(i) 
        enddo        !  i=2, gr%nnzp

        return
        end subroutine read_sounding

!------------------------------------------------------------------------
! Function to linearly interpolate quantities to a height zmid
! given the values vtop at a height ztop and vbot at a height
! zbot
!------------------------------------------------------------------------
        real function linint( zmid, ztop, zbot, vtop, vbot )
        implicit none

! Input 

        real zmid, ztop, zbot, vtop, vbot

! Compute linear interpolation

        linint = ( (zmid-zbot)/(ztop-zbot) ) * (vtop-vbot) + vbot

        return
        end function linint

!------------------------------------------------------------------------
! Subroutine to initialize one generic model variable from file
!------------------------------------------------------------------------
        subroutine read_profile( fname, x )

        use grid_class

        implicit none

! External

        real linint
        external linint

! Parameter

        integer, parameter :: nmaxsnd = 200

! Input

        character*(*) fname

! Output

        real x(1:gr%nnzp)

! Input variables from namelist

        integer nlevels
        real z(nmaxsnd),var(nmaxsnd)

        namelist /profile/ nlevels, z, var

! Internal

        integer i,k

! Read sounding namelist

        open(10, file = trim(fname), status = 'old' )
        read(10, nml = profile )
        close(10)

        if ( nlevels > nmaxsnd ) then
           write(*,*) 'Error in sounding: nlevels > nmaxsnd'
           write(*,*) 'nlevels = ',nlevels
           write(*,*) 'nmaxsnd = ',nmaxsnd
           stop 'STOP in sounding'
        endif

! Error check: if lowest themodynamic grid height is lower than the
! lowest value from the input sounding, then the linear interpolation
! scheme will fail

        if ( gr%zt(2) < z(1) ) then
           write(*,*) 'First level of input sounding must be',
     .     ' below first thermodynamic level'
           write(*,*) ' first sounding level z(1) = ',z(1)
           write(*,*) ' first thermodynamic level gr%zt(2) = ',gr%zt(2)
           stop 'STOP in sounding'
        endif

! Use linear interpolation from two nearest prescribed grid points
! (one above and one below) to initialize mean quantities in the model
! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

        do i = 2, gr%nnzp
          k = 1
!  40      if ( z(k) >= gr%zt(i) ) then
          do while ( z(k) < gr%zt(i) )
            k = k + 1
              if ( k > nlevels ) then
                write(*,*) 'STOP Not enough sounding data to ',
     .                     'initialize grid:'
                write(*,'(a,f7.1,/a,f7.1)') '  highest sounding level'
     .          ,z(nlevels)
     .          ,'  should be higher than highest thermodynamic point'
     .          ,gr%zt(gr%nnzp)

                write(*,*) ' filename: ',fname
                stop 'STOP in read_profile'
              endif
            x(i) = linint( gr%zt(i), z(k), z(k-1), var(k), var(k-1) )
          enddo ! while
!          else
!            goto 40
!          endif
        end do ! i=2, gr%nzzp

        return
        end subroutine read_profile
