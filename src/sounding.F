! $Id: sounding.F,v 1.15 2008-03-20 19:59:54 dschanen Exp $
#define SCLR_THETA 1
#define SCLR_RT 2
        module sounding

        implicit none

        public :: 
     .  read_sounding,
     .  read_profile ! Not currently used in HOC

        private :: linint

        contains
!------------------------------------------------------------------------
        subroutine read_sounding( iunit, 
     .                            thlm, rtm, um, vm, ugm, vgm, runfile, 
     .                            runtype, 
     .                            sclrm, edsclrm ) ! optional

!       Description:
!       Subroutine to initialize model variables from a namelist file
!------------------------------------------------------------------------

        use grid_class
        use constants, only: sclr_dim, fstderr

        implicit none

        ! Constant parameter
        integer, parameter :: nmaxsnd = 600


        ! Input variables
        integer, intent(in) :: iunit ! File unit to use for namelist

        character(len=*), intent(in) :: 
     .  runfile,  ! Contains namelists
     .  runtype   ! String for DYCOMS II RF02

        ! Output variables
        real, intent(out), dimension(gr%nnzp) :: 
     .  thlm, ! Liquid potential temperature    [K]
     .  rtm,  ! Total water mixing ratio        [kg/kg]
     .  um,   ! u wind                          [m/s]
     .  vm,   ! v wind                          [m/s]
     .  ugm,  ! u geostrophic wind              [m/s]
     .  vgm   ! v geostrophic wind              [m/s]

        ! Optional output variables
        real, optional, intent(out), dimension(gr%nnzp, sclr_dim) :: 
     .  sclrm, edsclrm ! Passive scalar input   [units vary] 

        ! Local variables

        ! Input variables from namelist
        integer :: nlevels  ! Levels in the input sounding

        real, dimension(nmaxsnd) ::
     .  z,     ! Altitude                               [m]
     .  theta, ! Liquid potential temperature sounding  [K]
     .  rt,    ! Total water mixing ratio sounding      [kg/kg]
     .  u,     ! u wind sounding                        [m/s] 
     .  v,     ! v wind sounding                        [m/s]
     .  ug,    ! u geostrophic wind sounding            [m/s]
     .  vg     ! v geostrophic wind sounding            [m/s]

        real, dimension(nmaxsnd, sclr_dim) :: 
     .  sclr, edsclr ! Passive scalar input sounding    [units vary]

        ! Namelists
        namelist /sounding/ nlevels, z, theta, rt, u, v, ug, vg

        namelist /scalar_sounding/ sclr, edsclr

        integer :: i, j, k  ! Loop indices

        ! Whether this routine needs to read in data for the 
        ! passive scalars variables
        logical :: sclr_sounding 

        ! Read sounding namelist
        if ( present( sclrm ) .and. present( edsclrm ) ) then 
          sclr_sounding = .true.
        else
          sclr_sounding = .false.
        end if

        
        open(unit = iunit, file = runfile, status = 'old')
        read(unit = iunit, nml = sounding)

        if ( sclr_sounding ) read(unit = iunit, nml = scalar_sounding)

        close(unit=iunit)

        if ( nlevels > nmaxsnd ) then
           write(fstderr,*) 'Error in sounding: nlevels > nmaxsnd'
           write(fstderr,*) 'nlevels = ',nlevels
           write(fstderr,*) 'nmaxsnd = ',nmaxsnd
           stop 'STOP in sounding'
        end if

        ! Error check: if lowest themodynamic grid height is lower than the
        ! lowest value from the input sounding, then the linear interpolation
        ! scheme will fail

        if ( gr%zt(2) < z(1) ) then
           write(fstderr,*) 'First level of input sounding must be',
     .     ' below first thermodynamic level'
           write(fstderr,*) ' first sounding level z(1) = ',z(1)
           write(fstderr,*) ' first thermodynamic level gr%zt(2) = ',
     .       gr%zt(2)
           stop 'STOP in sounding'
        end if

        ! First sounding level should be near ground value

        ! dschanen 1 May 2007
        ! We have changed this for Nov. 11 and June 25, both of which
        ! begin above the ground.
!       if ( abs(z(1)) > 1.e-8 ) then
        if ( .false. ) then
           write(fstderr,*) 'First level of input sounding must be z=0'
           stop 'STOP in sounding'
        else
          um(1)   = u(1)
          vm(1)   = v(1)
          ugm(1)  = ug(1)
          vgm(1)  = vg(1)
          thlm(1) = theta(1)
          rtm(1)  = rt(1)
          if ( sclr_sounding ) sclrm(1,:)   = sclr(1,:)
          if ( sclr_sounding ) edsclrm(1,:) = edsclr(1,:)
        end if

! Use linear interpolation from two nearest prescribed grid points
! (one above and one below) to initialize mean quantities in the model
! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

        do i=2, gr%nnzp
          k=1
          do while ( z(k) < gr%zt(i) )
            k=k+1
            if ( k > nlevels ) then
              write(fstderr,*) 'STOP Not enough sounding data to ',
     .          'initialize grid:'
              write(fstderr,'(a,f7.1,/a,f7.1)') 
     .          '  highest sounding level', z(nlevels),
     .          '  should be higher than highest thermodynamic point',
     .          gr%zt(gr%nnzp)
              stop 'STOP in sounding'
            end if  ! k > nlevels

            ! Regular situation w/ linear int.
            IF (      trim( runtype ) /= "dycoms2_rf02_do" 
     .          .AND. trim( runtype ) /= "dycoms2_rf02_ds" 
     .          .AND. trim( runtype ) /= "dycoms2_rf02_nd" 
     .          .AND. trim( runtype ) /= "dycoms2_rf02_so" ) THEN  

              um(i)   = linint( gr%zt(i), z(k), z(k-1), u(k), u(k-1) )
              vm(i)   = linint( gr%zt(i), z(k), z(k-1), v(k), v(k-1) )
              ugm(i)  = linint( gr%zt(i), z(k), z(k-1), ug(k), ug(k-1) )
              vgm(i)  = linint( gr%zt(i), z(k), z(k-1), vg(k), vg(k-1) )
              thlm(i) = linint( gr%zt(i), z(k), z(k-1), 
     .                      theta(k), theta(k-1) )
              rtm(i)  = linint( gr%zt(i), z(k), z(k-1), rt(k), rt(k-1) )
              if ( sclr_sounding ) then
                do j = 1, sclr_dim 
                  sclrm(i,j) = linint( gr%zt(i), z(k), z(k-1), 
     .                                 sclr(k,j), sclr(k-1,j) )
                  edsclrm(i,j) = linint( gr%zt(i), z(k), z(k-1), 
     .                                   edsclr(k,j), edsclr(k-1,j) )
                end do

              end if

            ELSE  ! DYCOMS II RF02 case

              IF ( gr%zt(i) < 795.0 ) THEN
                um(i)   =  3.0 + (4.3*gr%zt(i))/1000.0
                vm(i)   = -9.0 + (5.6*gr%zt(i))/1000.0
                ugm(i)  = um(i)
                vgm(i)  = vm(i)
                thlm(i) = 288.3
                rtm(i)  = (9.45)/1000.0
!               Scalars; change if not equal to th and rt
                if ( sclr_sounding ) then
                  sclrm(i, SCLR_THETA)   = thlm(i)
                  sclrm(i, SCLR_RT)      = rtm(i)
                  edsclrm(i, SCLR_THETA) = thlm(i)
                  edsclrm(i, SCLR_RT)    = rtm(i)
                end if
              ELSE
                um(i)   =  3.0 + (4.3*gr%zt(i))/1000.0
                vm(i)   = -9.0 + (5.6*gr%zt(i))/1000.0
                ugm(i)  = um(i)
                vgm(i)  = vm(i)
                thlm(i) = 295.0 + ( (gr%zt(i) - 795.0)**(1.0/3.0) )
                rtm(i)  = (  5.0 - 3.0 
     .          * ( 1.0 - EXP( (795.0 - gr%zt(i))/500.0 ) )  )/1000.0
                ! Same as above
                if ( sclr_sounding ) then
                  sclrm(i, SCLR_THETA)   = thlm(i)
                  sclrm(i, SCLR_RT)      = rtm(i)
                  edsclrm(i, SCLR_THETA) = thlm(i)
                  edsclrm(i, SCLR_RT)    = rtm(i)
                end if
              END IF

            END IF ! runtype

          end do ! do while ( z(k) < gr%zt(i) )

        end do   ! i=2, gr%nnzp

        return
        end subroutine read_sounding

!------------------------------------------------------------------------
        real function linint( zmid, ztop, zbot, vtop, vbot )

!       Description:
!       Function to linearly interpolate quantities to a height zmid
!       given the values vtop at a height ztop and vbot at a height
!       zbot.
!------------------------------------------------------------------------

        implicit none

        ! Input  Variables

        real, intent(in) :: zmid, ztop, zbot, vtop, vbot

        ! Compute linear interpolation

        linint = ( (zmid-zbot)/(ztop-zbot) ) * (vtop-vbot) + vbot

        return
        end function linint

!------------------------------------------------------------------------
        subroutine read_profile( fname, x )

!       Description:
!       Subroutine to initialize one generic model variable from file
!------------------------------------------------------------------------

        use grid_class

        implicit none

        ! Constant Parameter

        integer, parameter :: nmaxsnd = 200

        ! Input Variables

        character*(*) fname

        ! Output

        real, dimension(gr%nnzp) :: x

        ! Local variables

        ! Input variables from namelist

        integer :: nlevels

        real, dimension(nmaxsnd) :: z, var

        namelist /profile/ nlevels, z, var

        ! Loop indices
        integer :: i,k

        ! Read sounding namelist

        open(10, file = trim(fname), status = 'old' )
        read(10, nml = profile )
        close(10)

        if ( nlevels > nmaxsnd ) then
           write(*,*) 'Error in sounding: nlevels > nmaxsnd'
           write(*,*) 'nlevels = ',nlevels
           write(*,*) 'nmaxsnd = ',nmaxsnd
           stop 'STOP in sounding'
        end if

! Error check: if lowest themodynamic grid height is lower than the
! lowest value from the input sounding, then the linear interpolation
! scheme will fail

        if ( gr%zt(2) < z(1) ) then
           write(*,*) 'First level of input sounding must be',
     .     ' below first thermodynamic level'
           write(*,*) ' first sounding level z(1) = ',z(1)
           write(*,*) ' first thermodynamic level gr%zt(2) = ',gr%zt(2)
           stop 'STOP in sounding'
        end if

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
              end if
            x(i) = linint( gr%zt(i), z(k), z(k-1), var(k), var(k-1) )
          enddo ! while
        end do ! i=2, gr%nzzp

        return
        end subroutine read_profile

        end module sounding
!-----------------------------------------------------------------------
