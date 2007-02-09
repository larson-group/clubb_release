!------------------------------------------------------------------------
! $Id: numer_check.F,v 1.4 2007-02-09 22:13:43 dschanen Exp $
        module numerical_check

        implicit none

        public :: invalid_model_arrays
        private :: isnan2d

        contains

!------------------------------------------------------------------------
        logical function invalid_model_arrays( ) 

!       Description:
!       Checks for invalid floating point values in select model arrays.

!       References:
!       None
!------------------------------------------------------------------------

        use diagnostic_variables
        use prognostic_variables
        use model_flags, only: lcoamps_micro, kk_rain
        use constants, only: sclr_dim, hydromet_dim, fstderr

        implicit none

        ! Local Variables
        integer :: i

        invalid_model_arrays = .false.
           
        ! Check whether any variable array contains a NaN for
        ! um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp, 
        ! wp2, & wp3.       

        if ( isnan2d( um ) ) then
          write(fstderr,*) "NaN in um model array" 
!         write(fstderr,*) "um= ", um
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( vm ) ) then
          write(fstderr,*) "NaN in vm model array" 
!         write(fstderr,*) "vm= ", vm
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wp2 ) ) then
          write(fstderr,*) "NaN in wp2 model array" 
!         write(fstderr,*) "wp2= ", wp2
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wp3 ) ) then
          write(fstderr,*) "NaN in wp3 model array" 
!         write(fstderr,*) "wp3= ", wp3
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtm ) ) then
          write(fstderr,*) "NaN in rtm model array" 
!         write(fstderr,*) "rtm= ", rtm
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( thlm ) ) then
          write(fstderr,*) "NaN in thlm model array" 
!         write(fstderr,*) "thlm= ", thlm
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtp2 ) ) then
          write(fstderr,*) "NaN in rtp2 model array" 
!         write(fstderr,*) "rtp2= ", rtp2
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( thlp2 ) ) then
          write(fstderr,*) "NaN in thlp2 model array" 
!         write(fstderr,*) "thlp2= ", thlp2
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wprtp ) ) then
          write(fstderr,*) "NaN in wprtp model array" 
!         write(fstderr,*) "wprtp= ", wprtp
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( wpthlp ) ) then
          write(fstderr,*) "NaN in wpthlp model array" 
!         write(fstderr,*) "wpthlp= ", wpthlp
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtpthlp ) ) then
          write(fstderr,*) "NaN in rtpthlp model array" 
!         write(fstderr,*) "rtpthlp= ", rtpthlp
          invalid_model_arrays = .true.
!         return
        end if

        if ( kk_rain .or. lcoamps_micro ) then
          do i = 1, hydromet_dim, 1
            if ( isnan2d( hydromet(:,i) ) ) then
              write(fstderr,*) "NaN in a hydrometeor model array" 
!             write(fstderr,*) "hydromet= ", hydromet
              invalid_model_arrays = .true.
!             return
            end if
          end do
        end if

!       if ( isnan2d( wmt ) ) then
!         write(fstderr,*) "NaN in wmt model array" 
!         write(fstderr,*) "wmt= ", wmt
!         invalid_model_arrays = .true.
!         return
!       end if

        if ( isnan2d( wp2thvp ) ) then
          write(fstderr,*) "NaN in wp2thvp model array" 
!         write(fstderr,*) "wp2thvp = ", wp2thvp
          invalid_model_arrays = .true.
!         return
        end if

        if ( isnan2d( rtpthvp ) ) then
          write(fstderr,*) "NaN in rtpthvp model array" 
!         write(fstderr,*) "rtpthvp = ", rtpthvp
          invalid_model_arrays = .true.
        end if

        if ( isnan2d( thlpthvp ) ) then
          write(fstderr,*) "NaN in thlpthvp model array" 
!         write(fstderr,*) "thlpthvp = ", thlpthvp
          invalid_model_arrays = .true.
        end if

#ifdef SCALARS
        do i = 1, sclr_dim, 1
          if ( isnan2d( sclrm(:,i) ) ) then
            write(fstderr,*) "NaN in sclrm", i, "model array"
!           write(fstderr,'(a6,i2,a1)') "sclrm(", i, ")"
!           write(fstderr,*) sclrm(:,i)
            invalid_model_arrays = .true.
          end if
          if ( isnan2d( edsclrm(:,i) ) ) then
            write(fstderr,*) "NaN in edsclrm", i, "model array"
!           write(fstderr,'(a8,i2,a1)') "edsclrm(", i, ")"
!           write(fstderr,*) edsclrm(:,i)
            invalid_model_arrays = .true.
          end if
        end do
#endif /*SCALARS*/

        return
        end function invalid_model_arrays

!------------------------------------------------------------------------
        logical function isnan2d( x2d )

!       Description:
!       Checks if a given real is a NaN, +inf or -inf.

!       Notes: 
!       Got a tip from Andy Vaught to use transfer( ) for portability.
!       Ideally this would be done with intrinsic module ieee in 
!       Fortran 2003 or with the isnan( ) extension provided by some 
!       (not pgroup) compilers.

!       Certain compiler optimizations may cause variables with invalid
!       results to flush to zero.  Avoid these!
!------------------------------------------------------------------------
        
        implicit none

        ! External
        intrinsic :: transfer, any
    
        ! Input Variables
        real, dimension(:), intent(in) :: x2d

        ! Local Variables
        integer(kind=4) :: nanbits

        real(kind=4) :: PosInf, NegInf

        data nanbits /Z"7F800000"/

        PosInf = transfer( nanbits, PosInf )
        NegInf = -( transfer( nanbits, NegInf ) )

        ! This works on compilers with standardized floating point, 
        ! because the IEEE 754 spec defines that subnormals and NaNs 
        ! should not equal themselves.
        ! However, all compilers do not seem to follow this. 
        if ( any( x2d /= x2d ) ) then
          isnan2d = .true.

        ! This a second check, assuming the above does not work as 
        ! expected.
        else if ( any( x2d == PosInf ) ) then
          isnan2d = .true.

        ! This may never be needed, it's here just in case.
        else if ( any( x2d == NegInf ) ) then
          isnan2d = .true.

        else 
          isnan2d = .false. ! Yippee, real numbers.

        end if

        return
        end function isnan2d
!------------------------------------------------------------------------
        end module numerical_check
