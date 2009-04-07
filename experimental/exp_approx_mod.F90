! $Id$
!-------------------------------------------------------------------------------
module exp_approx_mod

  implicit none

  public :: exp

  interface exp
    module procedure dp_exp, sp_exp
  end interface

  private :: dp_exp, sp_exp

  private

  integer, parameter :: &
    r4 = 4, r8 = 8, & ! 4 byte real, 8 byte real
    i4 = 4            ! 4 byte integer

  contains

!-------------------------------------------------------------------------------
  elemental function dp_exp( x )

! Description:
!   Compute e^x using some properties of the IEEE-754 floating point format
!   to obtain the result.

! References:
!   Original C language algorithm:
!   ``A Fast, Compact Approximation of the Exponential Function''
!   Neural Computation 11(3), Schraudolph.

!   Fortran Implementation:
!   Annals of Nuclear Energy 31 (2004) 1027-1037, Technical note, 
!   ``Computational efficiencies of approximated exponential functions 
!     for transport calclulations of the characteristics method'', 
!     Yamamoto, et al.
!   Modified for adhere to Fortran 95 standards by Dave Schanen, UWM.
!-------------------------------------------------------------------------------

    use endian, only: big_endian

    implicit none

    ! External
    intrinsic :: int, transfer, huge, tiny

    ! Constant parameters
    real(kind=r8), parameter :: &
      exp_a = 1048576.0_r8 / 0.69314718056_r8

    integer(kind=i4), parameter :: &
      exp_c = 1072693248_i4 - 60801_i4

    ! Input
    real(kind=r8), intent(in) :: x

    ! Output
    real(kind=r8) :: dp_exp

    ! Local variables
    integer(kind=i4), dimension(2) :: exi

    ! --- Begin Code ---

    ! Handle the overflow/underflow condition
    if ( x >= 700._r8 ) then
      dp_exp = huge( dp_exp )
      return
    else if ( x <= -700._r8 ) then
      dp_exp = tiny( dp_exp )
      return
    end if

    if ( big_endian ) then
      exi(1) = int( exp_a*x ) + exp_c      
      exi(2) = 0
    else ! little_endian
      exi(1) = 0
      exi(2) = int( exp_a*x ) + exp_c      
    end if

    dp_exp = transfer( exi, dp_exp )

    return
  end function dp_exp

!-------------------------------------------------------------------------------
  elemental function sp_exp( x )
! Description:
!   The same algorithm as above, but written for a single precision argument
!   and single precision result.
!
! References:
!   None
!-------------------------------------------------------------------------------

     use endian, only: big_endian

    implicit none

    ! External
    intrinsic :: int, transfer

    ! Constant parameters
    real(kind=r4), parameter :: &
      exp_a = 1048576.0_r4 / 0.69314718056_r4

    integer(kind=i4), parameter :: &
      exp_c = 1072693248_i4 - 60801_i4

    ! Input
    real(kind=r4), intent(in) :: x

    ! Output
    real(kind=r4) :: sp_exp

    ! Local variables
    integer(kind=i4), dimension(2) :: exi

    real(kind=r8) :: tmp_exp

    ! --- Begin Code ---

    ! Handle the overflow/underflow condition
    if ( x >= 700._r4 ) then
      sp_exp = huge( sp_exp )
      return
    else if ( x <= -700._r4 ) then
      sp_exp = tiny( sp_exp )
      return
    end if

    if ( big_endian ) then
      exi(1) = int( exp_a*x ) + exp_c      
      exi(2) = 0
    else ! little_endian
      exi(1) = 0
      exi(2) = int( exp_a*x ) + exp_c      
    end if

    tmp_exp = transfer( exi, tmp_exp )
    sp_exp = real( tmp_exp, kind=r4)

    return
  end function sp_exp

end module exp_approx_mod
