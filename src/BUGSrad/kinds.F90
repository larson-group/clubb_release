! CVS:  $Id$
! CVS:  $Name: not supported by cvs2svn $

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module kinds

!***********************************************************************
!
!     This module defines variable precision for all common data
!     types.
!
!     Modifications:
!     Modified by David Schanen for HOC to allow for compile time 
!     promotion (e.g. to double and quad precision)
!-----------------------------------------------------------------------

implicit none

!-----------------------------------------------------------------------

! Additions
integer, private :: ct_int
logical, private :: ct_log

real, private :: ct_real

double precision, private :: ct_dble

integer, parameter :: char_len  = 80,                    &
                      int_kind  = kind(ct_int),          &
                      log_kind  = kind(ct_log),          &
                      real_kind = kind(ct_real),         &
                      dbl_kind  = kind(ct_dble)

end module kinds

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
