! $Id: sounding.F,v 1.9 2005-07-25 22:28:43 dschanen Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
      module sounding

      implicit none

      public read_sounding
      private linint
      public read_profile ! not currently used at all

      contains
!------------------------------------------------------------------------
! subroutine read_sounding( )

! Subroutine to initialize model variables from namelist file
!------------------------------------------------------------------------
        subroutine read_sounding( thlm, rtm, um, vm, ugm, vgm, runfile, 
     .                            runtype, 
     .                            sclrm ) ! optional

        use grid_class
        use constants, only: sclrm_dimension

        implicit none

! Input

        character(len=*), intent(in) :: runfile
        character(len=*), intent(in) :: runtype

! Parameter

        integer, parameter :: nmaxsnd = 200

! Output

        real, intent(out), dimension(gr%nnzp) :: thlm
        real, intent(out), dimension(gr%nnzp) :: rtm
        real, intent(out), dimension(gr%nnzp) :: um
        real, intent(out), dimension(gr%nnzp) :: vm
        real, intent(out), dimension(gr%nnzp) :: ugm
        real, intent(out), dimension(gr%nnzp) :: vgm

! Optional Output
        real, optional, intent(out), 
     .        dimension(gr%nnzp, sclrm_dimension) :: sclrm

! Input variables from namelist

        integer nlevels
        real z(nmaxsnd),theta(nmaxsnd),rt(nmaxsnd)
        real u(nmaxsnd),v(nmaxsnd),ug(nmaxsnd),vg(nmaxsnd)

        real, dimension(nmaxsnd, sclrm_dimension)  :: sclr

        namelist /sounding/ nlevels, z, theta, rt, u, v, ug, vg
        namelist /scalar_sounding/ sclr
    
! Internal

        integer i, j, k
        logical  :: sclr_sounding

! Read sounding namelist
        if ( present( sclrm ) ) then 
          sclr_sounding = .true.
        else
          sclr_sounding = .false.
        endif
        
        open(10, file = runfile, status = 'old')
        read(10, nml = sounding)

        if (sclr_sounding) read(10, nml = scalar_sounding)

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
          if ( sclr_sounding ) sclrm(1,:) = sclr(1,:)  
        endif

! Use linear interpolation from two nearest prescribed grid points
! (one above and one below) to initialize mean quantities in the model
! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

        do i=2, gr%nnzp
          k=1
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

            ! Regular situation w/ linear int.
            IF (      trim(runtype) /= "dycoms2_rf02_d" 
     .          .AND. trim(runtype) /= "dycoms2_rf02_nd" ) THEN  
              um(i)   = linint( gr%zt(i), z(k), z(k-1), u(k), u(k-1) )
              vm(i)   = linint( gr%zt(i), z(k), z(k-1), v(k), v(k-1) )
              ugm(i)  = linint( gr%zt(i), z(k), z(k-1), ug(k), ug(k-1) )
              vgm(i)  = linint( gr%zt(i), z(k), z(k-1), vg(k), vg(k-1) )
              thlm(i) = linint( gr%zt(i), z(k), z(k-1), 
     .                      theta(k), theta(k-1) )
              rtm(i)  = linint( gr%zt(i), z(k), z(k-1), rt(k), rt(k-1) )
              if ( sclr_sounding ) then
                do j = 1, sclrm_dimension 
                  sclrm(i,j) = linint( gr%zt(i), z(k), z(k-1), 
     .                                 sclr(k,j), sclr(k-1,j) )
                enddo
              endif
            ELSE  ! DYCOMS II RF02 case

              IF ( gr%zt(i) < 795.0 ) THEN
                um(i)   =  6.5 + (2.0*gr%zt(i))/1000.0
                vm(i)   = -8.5 + (5.0*gr%zt(i))/1000.0
                ugm(i)  = um(i)
                vgm(i)  = vm(i)
                thlm(i) = 288.3
                rtm(i)  = (9.45)/1000.0
!               Scalars; change if not equal to th and rt
                if ( sclr_sounding ) then
                  sclrm(i, SCLR_THETA) = thlm(i)
                  sclrm(i, SCLR_RT)    = rtm(i)
                endif
              ELSE
                um(i)   =  6.5 + (2.0*gr%zt(i))/1000.0
                vm(i)   = -8.5 + (5.0*gr%zt(i))/1000.0
                ugm(i)  = um(i)
                vgm(i)  = vm(i)
                thlm(i) = 295.0 + ( (gr%zt(i) - 795.0)**(1.0/3.0) )
                rtm(i)  = (  5.0 - 3.0 
     .          * ( 1.0 - EXP( (795.0 - gr%zt(i))/500.0 ) )  )/1000.0
                ! Same as above
                if ( sclr_sounding ) then
                  sclrm(i, SCLR_THETA) = thlm(i)
                  sclrm(i, SCLR_RT)    = rtm(i)
                endif
              ENDIF

            ENDIF ! runtype
          enddo ! do while 
        enddo   ! i=2, gr%nnzp

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

!        real linint
!        external linint

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
        end do ! i=2, gr%nzzp

        return
        end subroutine read_profile

      end module sounding
!-----------------------------------------------------------------------
