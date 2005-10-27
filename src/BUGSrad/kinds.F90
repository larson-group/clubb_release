! CVS:  $Id: kinds.F90,v 1.1 2005-10-27 20:06:50 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module kinds

!***********************************************************************
!
!     This module defines variable precision for all common data
!     types.
!
!-----------------------------------------------------------------------

implicit none

!-----------------------------------------------------------------------

integer, parameter :: char_len  = 80,                       &
                         int_kind  = kind(1),               &
                         log_kind  = kind(.true.),          &
                         real_kind = selected_real_kind(6), &
                         dbl_kind  = selected_real_kind(13)    !13 for dble

end module kinds

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


