! $Id$
module sounding

  implicit none

  public ::  & 
    read_sounding, & 
    read_profile ! Not currently used in HOC
       
  private ! Default Scope
       
  contains
!------------------------------------------------------------------------
  subroutine read_sounding( iunit, runfile, runtype,  & 
                            thlm, rtm, um, vm, ugm, vgm, &
                            sclrm, edsclrm )

!       Description:
!       Subroutine to initialize model variables from a namelist file
!       References:
!       None
!------------------------------------------------------------------------

        use grid_class, only:  & 
            gr ! Variable(s)

        use constants, only:  & 
            fstderr ! Constant

        use parameters_model, only: & 
            sclr_dim ! Variable(s)

        use std_atmosphere_mod, only:  & 
            std_atmosphere ! Procedure(s)

        use interpolation, only:  & 
            lin_int ! Procedure(s)       

        use array_index, only: & 
            iisclr_rt, &  ! Variable
            iisclr_thl
!           ,iisclr_CO2

        use error_code, only: &
          clubb_at_least_debug_level ! Function

        implicit none

        
        ! Constant parameter
        integer, parameter :: nmaxsnd = 600
        integer, parameter :: sclr_max = 10

        ! Input variables
        integer, intent(in) :: iunit ! File unit to use for namelist

        character(len=*), intent(in) ::  & 
        runfile,   & ! Contains namelists
        runtype      ! String for DYCOMS II RF02

        ! Output variables
        real, intent(out), dimension(gr%nnzp) ::  & 
        thlm,  & ! Liquid potential temperature    [K]
        rtm,   & ! Total water mixing ratio        [kg/kg]
        um,    & ! u wind                          [m/s]
        vm,    & ! v wind                          [m/s]
        ugm,   & ! u geostrophic wind              [m/s]
        vgm      ! v geostrophic wind              [m/s]

        ! Optional output variables
        real, intent(out), dimension(gr%nnzp, sclr_dim) ::  & 
        sclrm, edsclrm ! Passive scalar input   [units vary] 

        ! Local variables

        ! Input variables from namelist
        integer :: nlevels  ! Levels in the input sounding

        real, dimension(nmaxsnd) :: & 
        z,      & ! Altitude                               [m]
        theta,  & ! Liquid potential temperature sounding  [K]
        rt,     & ! Total water mixing ratio sounding      [kg/kg]
        u,      & ! u wind sounding                        [m/s] 
        v,      & ! v wind sounding                        [m/s]
        ug,     & ! u geostrophic wind sounding            [m/s]
        vg     ! v geostrophic wind sounding            [m/s]

        real, dimension(nmaxsnd, sclr_max) ::  & 
        sclr, edsclr ! Passive scalar input sounding    [units vary]

        ! Namelists
        namelist /sounding/ nlevels, z, theta, rt, u, v, ug, vg

        namelist /scalar_sounding/ sclr, edsclr

        integer :: i, j, k  ! Loop indices

        ! Is this model being extended by 1976 Standard Atmosphere?
        logical :: l_std_atmo
        
        ! Read sounding namelist

        l_std_atmo = .false.
        
        open(unit = iunit, file = runfile, status = 'old')
        read(unit = iunit, nml = sounding)

        ! Read in a passive scalar sounding, if enabled
        if ( sclr_dim > 0 ) then
          ! Initialize to zero
          sclr   = 0.0
          edsclr = 0.0
          ! Read values from namelist
          read(unit = iunit, nml = scalar_sounding)
        end if

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
           write(fstderr,*) 'First level of input sounding must be', & 
           ' below first thermodynamic level'
           write(fstderr,*) ' first sounding level z(1) = ',z(1)
           write(fstderr,*) ' first thermodynamic level gr%zt(2) = ', & 
             gr%zt(2)
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
          if ( sclr_dim > 0 ) then
            sclrm(1,1:sclr_dim)   = sclr(1,1:sclr_dim)
            edsclrm(1,1:sclr_dim) = edsclr(1,1:sclr_dim)
          end if
        end if

        if ( clubb_at_least_debug_level( 1 ) ) then
          print *, "Reading in sounding information"
!------------------Printing Model Inputs-------------------------------
          print *, "u = ", u(1:nlevels)
          print *, "v = ", v(1:nlevels)
          print *, "ug = ", ug(1:nlevels)
          print *, "vg = ", vg(1:nlevels)
          print *, "theta = ", theta(1:nlevels)
          print *, "rt = ", rt(1:nlevels)
          do i = 1, sclr_dim, 1
            write(6,'(a5,i2,a2)',advance='no') "sclr(", i,")="
            write(6,'(8g10.2)') sclr(1:nlevels,i)
            write(6,'(a7,i2,a2)',advance='no') "edsclr(", i, ")="
            write(6,'(8g10.3)') edsclr(1:nlevels,i)
          end do
        end if ! clubb_at_least_debug_level( 1 )
!----------------------------------------------------------------------

! Use linear interpolation from two nearest prescribed grid points
! (one above and one below) to initialize mean quantities in the model
! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )
        do i=2, gr%nnzp
          k=1
          do while ( z(k) < gr%zt(i) .and. .not. l_std_atmo )
            k=k+1
            if ( k > nlevels ) then
!              write(fstderr,*) 'STOP Not enough sounding data to ',
!     .          'initialize grid:'
!              write(fstderr,'(a,f7.1,/a,f7.1)') 
!     .          '  highest sounding level', z(nlevels),
!     .          '  should be higher than highest thermodynamic point',
!     .          gr%zt(gr%nnzp)
!              stop 'STOP in sounding'

               l_std_atmo = .true.
               exit
            end if  ! k > nlevels

            ! Regular situation w/ linear int.
            IF (      trim( runtype ) /= "dycoms2_rf02_do"  & 
                .AND. trim( runtype ) /= "dycoms2_rf02_ds"  & 
                .AND. trim( runtype ) /= "dycoms2_rf02_nd"  & 
                .AND. trim( runtype ) /= "dycoms2_rf02_so" ) THEN  

              um(i)   = lin_int( gr%zt(i), z(k), z(k-1), u(k), u(k-1) )
              vm(i)   = lin_int( gr%zt(i), z(k), z(k-1), v(k), v(k-1) )
              ugm(i)  = lin_int( gr%zt(i), z(k), z(k-1), ug(k), ug(k-1) )
              vgm(i)  = lin_int( gr%zt(i), z(k), z(k-1), vg(k), vg(k-1) )
              thlm(i) = lin_int( gr%zt(i), z(k), z(k-1),  & 
                            theta(k), theta(k-1) )
              rtm(i)  = lin_int( gr%zt(i), z(k), z(k-1), rt(k), rt(k-1) )
              
              if ( sclr_dim > 0 ) then
                do j = 1, sclr_dim 
                  sclrm(i,j) = lin_int( gr%zt(i), z(k), z(k-1),  & 
                                       sclr(k,j), sclr(k-1,j) )
                  edsclrm(i,j) = lin_int( gr%zt(i), z(k), z(k-1),  & 
                                         edsclr(k,j), edsclr(k-1,j) )
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
                ! Passive Scalars
                ! Change this if they are not equal to theta_l and rt in RF02
                if ( iisclr_thl > 0  ) then
                  sclrm(i, iisclr_thl)   = thlm(i)
                  edsclrm(i, iisclr_thl) = thlm(i)
                end if
                if ( iisclr_rt > 0 ) then
                  sclrm(i, iisclr_rt)    = rtm(i)
                  edsclrm(i, iisclr_rt)  = rtm(i)
                end if
              ELSE
                um(i)   =  3.0 + (4.3*gr%zt(i))/1000.0
                vm(i)   = -9.0 + (5.6*gr%zt(i))/1000.0
                ugm(i)  = um(i)
                vgm(i)  = vm(i)
                thlm(i) = 295.0 + ( (gr%zt(i) - 795.0)**(1.0/3.0) )
                rtm(i)  = (  5.0 - 3.0  & 
                * ( 1.0 - EXP( (795.0 - gr%zt(i))/500.0 ) )  )/1000.0
                ! Passive Scalars
                ! Same as above
                if ( iisclr_thl > 0  ) then
                  sclrm(i, iisclr_thl)   = thlm(i)
                  edsclrm(i, iisclr_thl) = thlm(i)
                end if
                if ( iisclr_rt > 0 ) then
                  sclrm(i, iisclr_rt)    = rtm(i)
                  edsclrm(i, iisclr_rt)  = rtm(i)
                end if
              END IF

            END IF ! runtype

          end do ! do while ( z(k) < gr%zt(i) )

          ! If the grid extends beyond the sounding data, use
          ! Standard Atmosphere
          ! Joshua Fasching April 2009
          if ( l_std_atmo ) then
            call std_atmosphere( gr%zt(i), thlm(i), rtm(i) )
                  
            um(i) = um(i-1)
            vm(i) = vm(i-1)
            ugm(i) = um(i)
            vgm(i) = vm(i)

          end if

        end do   ! i=2, gr%nnzp
        
        if ( l_std_atmo ) then
          write(fstderr,*) "Warning: 1976 Standard Atmosphere "// & 
            "was used to complete the grid."
        end if

  return
  end subroutine read_sounding

!------------------------------------------------------------------------
        subroutine read_profile( fname, x )

!       Description:
!       Subroutine to initialize one generic model variable from file
!------------------------------------------------------------------------

        use grid_class, only:  & 
            gr ! Variable(s)
        use interpolation, only:  & 
            lin_int ! Procedure
        
        implicit none

        ! Constant Parameter

        integer, parameter :: nmaxsnd = 200

        ! Input Variables

        character*(*), intent(in) :: fname

        ! Output

        real, dimension(gr%nnzp), intent(out) :: x

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
           write(*,*) 'First level of input sounding must be', & 
           ' below first thermodynamic level'
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
                write(*,*) 'STOP Not enough sounding data to ', & 
                           'initialize grid:'
                write(*,'(a,f7.1,/a,f7.1)') '  highest sounding level' & 
                ,z(nlevels) & 
                ,'  should be higher than highest thermodynamic point' & 
                ,gr%zt(gr%nnzp)

                write(*,*) ' filename: ',fname
                stop 'STOP in read_profile'
              end if
            x(i) = lin_int( gr%zt(i), z(k), z(k-1), var(k), var(k-1) )
          enddo ! while
        end do ! i=2, gr%nzzp

        return
        end subroutine read_profile
        
        end module sounding
!-----------------------------------------------------------------------
