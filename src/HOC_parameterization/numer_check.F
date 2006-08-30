!------------------------------------------------------------------------
! $Id: numer_check.F,v 1.3 2006-08-30 03:16:24 dschanen Exp $
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

        implicit none

        invalid_model_arrays = .false.
           
        ! Check whether any variable array contains a NaN for
        ! um, vm, thlm, rtm, rtp2, thlp2, wprtp, wpthlp, rtpthlp, 
        ! wp2, & wp3.       

        if ( isnan2d( um ) ) then
          print *, "NaN in um model array" 
!         print *, "um= ", um
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( vm ) ) then
          print *, "NaN in vm model array" 
!         print *, "vm= ", vm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wp2 ) ) then
          print *, "NaN in wp2 model array" 
!         print *, "wp2= ", wp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wp3 ) ) then
          print *, "NaN in wp3 model array" 
!         print *, "wp3= ", wp3
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtm ) ) then
          print *, "NaN in rtm model array" 
!         print *, "rtm= ", rtm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( thlm ) ) then
          print *, "NaN in thlm model array" 
!         print *, "thlm= ", thlm
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtp2 ) ) then
          print *, "NaN in rtp2 model array" 
!         print *, "rtp2= ", rtp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( thlp2 ) ) then
          print *, "NaN in thlp2 model array" 
!         print *, "thlp2= ", thlp2
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wprtp ) ) then
          print *, "NaN in wprtp model array" 
!         print *, "wprtp= ", wprtp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wpthlp ) ) then
          print *, "NaN in wpthlp model array" 
!         print *, "wpthlp= ", wpthlp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtpthlp ) ) then
          print *, "NaN in rtpthlp model array" 
!         print *, "rtpthlp= ", rtpthlp
          invalid_model_arrays = .true.
!         return
        endif

!       if ( isnan2d( Nrm ) ) then
!         print *, "NaN in Nrm model array" 
!         print *, "Nrm= ", Nrm
!         invalid_model_arrays = .true.
!         return
!       endif

!       if ( isnan2d( rrm ) ) then
!         print *, "NaN in rrm model array" 
!         print *, "rrm= ", rrm
!         invalid_model_arrays = .true.
!         return
!       endif

        if ( isnan2d( wmt ) ) then
          print *, "NaN in wmt model array" 
!         print *, "wmt= ", wmt
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( wp2thvp ) ) then
          print *, "NaN in wp2thvp model array" 
!         print *, "wp2thvp = ", wp2thvp
          invalid_model_arrays = .true.
!         return
        endif

        if ( isnan2d( rtpthvp ) ) then
          print *, "NaN in rtpthvp model array" 
!         print *, "rtpthvp = ", rtpthvp
          invalid_model_arrays = .true.
        endif

        if ( isnan2d( thlpthvp ) ) then
          print *, "NaN in thlpthvp model array" 
!         print *, "thlpthvp = ", thlpthvp
          invalid_model_arrays = .true.
        endif

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
