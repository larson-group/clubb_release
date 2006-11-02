!-----------------------------------------------------------------------
! $Id: recl.F,v 1.1 2006-11-02 01:36:10 dschanen Exp $
!       Description:
!       Preprocessing rules for determining how large an unformatted 
!       data record is when using Fortran write.

!       Notes:
!       New directives will need to be added to port HOC GrADS output
!       to new compilers that do not use byte size record lengths.

!       Early Alpha processors lacked the ability to work with anything
!       smaller than a 32 bit word, so DEC Fortran and its successors
!       (Compaq Visual Fortran, newer Intel Fortran, etc.) all use 4 
!       byte records.  Note that specifying byterecl on Alpha still 
!       results in a performance hit, even on newer chips.
!-----------------------------------------------------------------------
#if defined (__i386__) && __INTEL_COMPILER >= 900 /* Newer ifort on i386 */
#define F_RECL 1
#elif defined (__x86_64__) && __INTEL_COMPILER >= 900 /*Newer ifort on x86_64*/
#define F_RECL 1
#elif __alpha /*Assume 4 byte word on Alpha processors*/
#define F_RECL 1
#else
#define F_RECL 4
#endif
