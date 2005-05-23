!----------------------------------------------------------------------
!  MODULE ENDIAN
!  contains:

!  big_endian( ) and little_endian( ) are
!  simple boolean functions to determine byte ordering on the system
!  One would be sufficent, the pair are available for convenience

!  native_4byte_real is a portable byte re-ordering subroutine
!----------------------------------------------------------------------
      module endian

      implicit none

      contains

!----------------------------------------------------------------------
!  Adapted from Chris Golaz's subroutine
!  Return .true. if the system uses most significant bit (MSB) byte ordering
!----------------------------------------------------------------------
      logical function big_endian( )

      implicit none

!---- Parameters ------------------------------------------------------
      integer, parameter :: ascii_0 = 48,
     .                      ascii_1 = 49,
     .                      ascii_2 = 50,
     .                      ascii_3 = 51

!---- Internal variables ----------------------------------------------
      integer(kind=4)    :: i
      logical            :: internal_endian

      i = ascii_0 + ascii_1*256 + ascii_2*(256**2) + ascii_3*(256**3)

      big_endian = internal_endian( i )

      return

      end function big_endian

!----------------------------------------------------------------------
!  Wrapper function
!  Return .true. if the system uses least significant bit (LSB) byte ordering
!----------------------------------------------------------------------
      logical function little_endian( )

      implicit none

      little_endian = .not. big_endian( )

      end function little_endian

!-----------------------------------------------------------------------
!     SUBPROGRAM: native_4byte_real
!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003 
!  LAST MODIFIED: 29 April 2003
!
!       SYNOPSIS: Converts a 32 bit, 4 byte, REAL from big Endian to
!                 little Endian, or conversely from little Endian to big
!                 Endian.
!
!    DESCRIPTION: This subprogram allows one to convert a 32 bit, 4 byte,
!                 REAL data element that was generated with, say, a big
!                 Endian processor (e.g. Sun/sparc, SGI/R10000, etc.) to its
!                 equivalent little Endian representation for use on little
!                 Endian processors (e.g. PC/Pentium running Linux). The
!                 converse, little Endian to big Endian, also holds.
!                 This conversion is accomplished by writing the 32 bits of
!                 the REAL data element into a generic 32 bit INTEGER space
!                 with the TRANSFER intrinsic, reordering the 4 bytes with
!                 the MVBITS intrinsic, and writing the reordered bytes into
!                 a new 32 bit REAL data element, again with the TRANSFER
!                 intrinsic. The following schematic illustrates the
!                 reordering process
!
!
!                  --------    --------    --------    --------
!                 |    D   |  |    C   |  |    B   |  |    A   |  4 Bytes
!                  --------    --------    --------    --------
!                                                             |
!                                                              -> 1 bit
!                                       ||
!                                     MVBITS
!                                       ||
!                                       \/
!
!                  --------    --------    --------    --------
!                 |    A   |  |    B   |  |    C   |  |    D   |  4 Bytes
!                  --------    --------    --------    --------
!                         |           |           |           |
!                         24          16          8           0   <- bit
!                                                                 position
!
!          INPUT: realIn,  a single 32 bit, 4 byte REAL data element.
!         OUTPUT: realOut, a single 32 bit, 4 byte REAL data element, with
!                 reverse byte order to that of realIn.
!    RESTRICTION: It is assumed that the default REAL data element is
!                 32 bits / 4 bytes.
!
!-----------------------------------------------------------------------
      SUBROUTINE native_4byte_real( realInOut )

      IMPLICIT NONE

      REAL, INTENT(INOUT):: realInOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
!      Modified 4/19/05 -dschanen
!
!      REAL, INTENT(IN):: realInOut
!      REAL, INTENT(OUT) :: realOut
!                                                   ! a single 32 bit, 4 byte
!                                                   ! REAL data element, with
!                                                   ! reverse byte order to
!                                                   ! that of realIn
!----------------------------------------------------------------------
! Local variables (generic 32 bit INTEGER spaces):

      INTEGER                                       :: i_element
      INTEGER                                       :: i_element_br
!----------------------------------------------------------------------
! Transfer 32 bits of realIn to generic 32 bit INTEGER space:
      i_element = TRANSFER( realInOut, 0 )
!----------------------------------------------------------------------
! Reverse order of 4 bytes in 32 bit INTEGER space:
      CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
      CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
      CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
      CALL MVBITS( i_element,  0, 8, i_element_br, 24 )
!----------------------------------------------------------------------
! Transfer reversed order bytes to 32 bit REAL space (realOut):
      realInOut = TRANSFER( i_element_br, 0.0 )

      RETURN
!----------------------------------------------------------------------
      END SUBROUTINE
!----------------------------------------------------------------------

      end module endian

!----------------------------------------------------------------------
!  This function is only called indirectly
!----------------------------------------------------------------------
      logical function internal_endian( i )
      implicit none

      character(len=4) :: i

      if ( trim(i) == '0123' ) then
        internal_endian = .false.  ! little_endian
      elseif ( trim(i) == '3210' ) then
        internal_endian = .true.   ! big_endian
      else
        stop "endian() failed"
      endif
 
      return

      end function internal_endian
!----------------------------------------------------------------------
