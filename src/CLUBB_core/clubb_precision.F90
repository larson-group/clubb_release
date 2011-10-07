!-------------------------------------------------------------------------------
! $Id$
module clubb_precision

  implicit none

  public :: stat_nknd, stat_rknd, time_precision, dp, sp

  private ! Default scope

  ! The precisions below are arbitrary, and could be adjusted as
  ! needed for long simulations or time averaging.  Note that on
  ! most machines 12 digits of precision will use a data type
  ! which is 8 bytes long.
  integer, parameter ::  & 
    stat_nknd = selected_int_kind( 8 ), & 
    stat_rknd = selected_real_kind( p=12 ), & 
    time_precision = selected_real_kind( p=12 ), &
    dp = selected_real_kind( p=12 ), & ! double precision
    sp = selected_real_kind( p=5 )     ! single precision

end module clubb_precision
!-------------------------------------------------------------------------------
