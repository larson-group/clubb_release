! CVS:  $Id: driver_read.F,v 1.1 2005-10-27 20:06:50 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $


!-----------------------------------------------------------------------
      program driver_read

      use kinds, only:  int_kind, dbl_kind
      use bugsrad_physconst, only:  gravity, cp_dry_air, sol_const

      implicit none

!-----------------------------------------------------------------------
! driver_read is the main routine for running the CSU radiative transfer
! code offline (that is, apart from the CSU GCM).  It reads a profile
! from a file and also specifies variables that are not read in.  It
! then calls BUGSrad to do the radiative transfer.  This driver is not
! used when the code is compiled online with the CSU GCM.

! REFERENCES:
! Phil Partain /wombat (04-04-00).

! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     bugs_rad   :The radiative transfer code

! FUNCTIONS CALLED:
!     none.
 
! INCLUDED COMMON BLOCKS:
!     none.
 
! LOCAL VARIABLES:
      integer (kind=int_kind)::
     &  nlen      !Length of total domain.
     &, len       !Length of sub domain.
     &, nlm       !Number of layers.
     &, i,l

      real (kind=dbl_kind), dimension(:), allocatable::
     &  ts         !Surface temperature                             (K).
     &, tbmn       !PBL top temperature                             (K). ! Kurt
     &, pc         !Cloud base pressure of PBL clouds             (hPa). ! Kurt
     &, amu0       !Cosine of solar zenith angle                    (-).
     &, slr        !Fraction of daylight                            (-).
     &, alvdr      !Visible direct surface albedo                   (-).
     &, alndr      !Near-IR direct surface albedo                   (-).
     &, alvdf      !Visible diffuse surface albedo                  (-).
     &, alndf      !Near-IR diffuse surface albedo                  (-).

      real (kind=dbl_kind), dimension(:,:), allocatable::
     &  pl         !Layer pressure                                (hPa). ! Kurt changed to (mb)
     &, dpl        !Layer thickness                               (hPa).
     &, tl         !Temperature                                     (K).
     &, ql         !Specific humidity                           (kg/kg).
     &, qcwl       !Cloud water mixing ratio                    (kg/kg).
     &, qcil       !Cloud ice mixing ratio                      (kg/kg).
     &, qrwl       !Rain mixing ratio                           (kg/kg).
     &, qril       !Snow mixing ratio                           (kg/kg).
     &, o3l        !Ozone mixing ratio                          (kg/kg).
     &, acld       !Radiative cloud fraction                        (-).

      !Note that rain mixing ratio is unused by BUGSrad, but I've put it
      !into the std_profile.dat file for completeness

      real (kind=dbl_kind), dimension(:,:), allocatable::
     &  pl2        !Level pressure                                (hPa). ! Kurt changed to (mb)

      real (kind=dbl_kind), dimension(:,:), allocatable::
     &  atl        !All-sky LW radiative heating rate             (K/s).
     &, asl        !All-sky SW radiative heating rate             (K/s).
     &, fulw       !All-sky LW upwelling flux                   (W/m^2).
     &, fdlw       !All-sky LW downwelling flux                 (W/m^2).
     &, fusw       !All-sky SW upwelling flux                   (W/m^2).
     &, fdsw       !All-sky SW downwelling flux                 (W/m^2).
      !For timing
!      real, dimension(2) :: tarray
!      real :: dtime, elapsed
!-----------------------------------------------------------------------

!---- 1. READ PROFILE DATA FROM FILE:
!      print *,'Congratulations... Your profile.dat is being read in!!!'
      open(10,file='profile.dat',action='read')
      read(10,*) nlen,nlm
!      nlen = 1    !to do timing tests
      len = nlen

!---- ALLOCATE ARRAYS
       allocate(ts(nlen) , tbmn(nlen) , pc(nlen) , amu0(nlen) 
     &, slr(nlen))

       allocate(alvdr(nlen) , alndr(nlen) , alvdf(nlen), alndf(nlen))

       allocate(pl(nlen,nlm), dpl(nlen,nlm), tl(nlen,nlm), ql(nlen,nlm)
     &, qcwl(nlen,nlm) , qcil(nlen,nlm), qrwl(nlen,nlm), qril(nlen,nlm)
     &, o3l(nlen,nlm), acld(nlen,nlm), dpl(nlen,nlm))

       allocate(pl2(nlen,nlm+1))

       allocate(atl(nlen,nlm), asl(nlen,nlm), fulw(nlen,nlm+1)
     &, fdlw(nlen,nlm+1), fusw(nlen,nlm+1), fdsw(nlen,nlm+1))
!----
      do l=1,nlm

! Kurt 

!        read(10,*) i,pl(1,l),pl2(1,l),tl(1,l),ql(1,l),o3l(1,l),
!     &             qcwl(1,l), qcil(1,l), qrwl(1,l), qril(1,l), acld(1,l)
!
!        pl2(1,l) = pl2(1,l)/100.       !convert from Pascals to millibars
!        pl(1,l)  = pl(1,l)/100.        !convert from Pascals to millibars

        read(10,*) i,pl(1,l),pl2(1,l),tl(1,l),ql(1,l),o3l(1,l)
     &, qcwl(1,l),qcil(1,l),acld(1,l),dpl(1,l)

! Kurt says:
!        ql(1,l) = ql(1,l)/100.
        o3l(1,l) = o3l(1,l)/100000.

! Kurt made the profile read data in units of mb, so below lines are no
! Longer needed!  23 JUN 05

!        pl2(1,l) = pl2(1,l)/1.       !convert from Pascals to millibars
!        pl(1,l)  = pl(1,l)/1.        !convert from Pascals to millibars

      enddo

! Kurt Kotenberg does not want alvdr, etc. to be read in, so commented out:
      read(10,*) pl2(1,nlm+1), ts(1), amu0(1)
!      read(10,*) pl2(1,nlm+1),ts(1), amu0(1), alvdr(1), alvdf(1),
!     &           alndr(1), alndf(1)

! End comments

! Kurt made the profile read data in units of mb, so below lines are no
! Longer needed!  23 JUN 05

!      pl2(1,nlm+1) = pl2(1,nlm+1)/100. !convert from Pascals to millibars
!     pl2(1,nlm+1) = pl2(1,nlm+1)/1. !convert from Pascals to millibars
      
      close(10)
           

! Kurt says these values are inputted for the Nov11 case, therefore, below
! computations not necessary. 16 JUN 05

!      do l=1,nlm
!        dpl(1,l) = pl2(1,l+1)-pl2(1,l)
!      enddo

! clouds?  Hardcoded here, can read in if you want to.
!
!      qcwl(1,:) = 0.0
!      qcil(1,:) = 0.0
!      qril(1,:) = 0.0
!
!      acld(1,:) = 0.0
!      acld(1,8) = 0.2
!      acld(1,9) = 0.3
!      acld(1,12) = 0.25
!      acld(1,13) = 0.25
!      acld(1,17) = 0.4
!      acld(1,18) = 0.3
!
!      qcil(1,8) = 0.0003
!      qcil(1,9) = 0.0003
!      qcil(1,12) = 0.0003
!      qcil(1,13) = 0.0003
!      qcil(1,17) = 0.0003
!      qcil(1,18) = 0.0003
!

! End Kurt's comments


!---- 2. COPY PROFILE TO ALL COLUMNS:
! copy the same column to all columns (only useful if testing multiple
! identical columns for timing, otherwise, it doesn't hurt)
      do i=1,nlen       
        pl2(i,:) = pl2(1,:)
        pl(i,:) = pl(1,:)
        dpl(i,:) = dpl(1,:)
        tl(i,:) = tl(1,:)
        ql(i,:) = ql(1,:)
        o3l(i,:) = o3l(1,:)
        acld(i,:) = acld(1,:)
        qcwl(i,:) = qcwl(1,:)
        qcil(i,:) = qcil(1,:)
        qril(i,:) = qril(1,:)
        amu0(i) = amu0(1)
        alvdr(i) = alvdr(1)
        alvdf(i) = alvdf(1)
        alndr(i) = alndr(1)
        alndf(i) = alndf(1)
      enddo

! Kurt Kotenberg added values for 11 Nov 99 altocumulus case
!Nov 11th values...
!      ts(:) = 289.0
      tbmn(:) = 286.07585
      pc(:) = 985.634

! Kurt Kotenberg commented amu0 out, since it is being read
! in by the profile.dat

!        amu0 = 0.4329  !actual Nov 11th angle...
!      amu0(:) = 1.0
!      amu0 = 0.1

! End additions

!---- 3. SPECIFY OTHER VARIABLES:
! Kurt Kotenberg uncommented these values since they are not being read
! in by the profile.dat

!      alvdr(:) = 0.2
!      alvdf(:) = 0.2
!      alndr(:) = 0.2
!      alndf(:) = 0.2

!      slr(:) = 1.0

      ! Kurt Kotenberg then changed the values from above to the values below per Adam's calculations for the actual surface albedos.  11/03/04

      alvdr(:) = 0.1
      alvdf(:) = 0.1
      alndr(:) = 0.1
      alndf(:) = 0.1


      slr(:) = 1.0
! End Kurt Kotenberg's changes


!---- 4. CALL THE RADIATIVE TRANSFER CODE:
!      elapsed = dtime(tarray)
      call bugs_rad(nlen,len,nlm,pl2,pl,dpl,tl,ql,qcwl,qcil,qril
     +,              o3l,ts,amu0,slr,alvdf,alndf,alvdr,alndr,sol_const
     +,              gravity,cp_dry_air,asl,atl,fdsw,fusw,fdlw,fulw,
     +               acld)
!      elapsed = dtime(tarray)

!---- 5. OUTPUT RESULTS:
! print fluxes in W/m2, heating rates in K/day.
!      print *, "Dtime: ", elapsed
!      print *, 'Fluxes'
!      do l=1,nlm
!        print '(I4,2F13.4,(F13.7),(F13.4))',l,fdsw(1,l),fusw(1,l)
!     &,                                  fdlw(2,l),fulw(2,l)
!      enddo
!      print *, 'Heating Rates'
!      do l=1,nlm
!        print '(I4,F12.3,2(F15.5))',
!     &         l,pl(1,l),asl(1,l)*86400.,atl(1,l)*86400. !K/day
!      enddo

! Kurt Kotenberg commented out the above and added the following so that a
! *.dat file is written for the fluxes and heating rates:


! For Fluxes:

      open (unit=20, file="bugradFCLEX.dat")
        do l = 1, nlm
       write (20,'(I4,2F13.4,(F13.7),(F13.4))'),
     &  l,fdsw(1,l),fusw(1,l),fdlw(1,l),fulw(1,l)
!Cadam       write (9,*)l,pl2_loc(2,l),fdsw(2,l),fusw(2,l),fdlw(2,l),fulw(2,l)
      enddo
        close (unit=20)


! For heating rates:

      open (unit=22, file="bugradHCLEX.dat")
Cadam      open (unit=10, file="C\work\rad_model\bugsrad\bugoutrad2.dat")
      do l=1,nlm
       write (22, '(I4,F13.7,F12.2,2(F15.5))'),
     &         l,qcwl(1,l),pl2(1,l),asl(1,l)*86400.,atl(1,l)*86400. 
! Kurt Kotenberg says uncomment the following to add dpl to the outputted file
! (used in computations for single band heating rate)

!     &,         dpl(1,l) !K/day
! end changes

      enddo
!         write (10,*) cldt2d,SWinc,SWnettoa,SWnetsfc,SWdnsfc,LWtoa,
!     &         LWdnsfc,LWnetsfc,radvbc,radvdc,radnbc,radndc
      close (unit=22)

! End Kurt's changes


! Kurt Kotenberg added the following print statements to check to be 
! certain units and correct values are being read in properly.
! Once verified, should be commented out, since unnecessary.
!      do l=1,nlm
!      print '(I4,F15.6)', l,tl(1,l)

!       print *, o3l
!      enddo

      end program driver_read

!-----------------------------------------------------------------------
