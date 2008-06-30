      module random

      implicit none

      public :: rand_permute, ran2

      private ! Defualt Scope
      
      contains
!----------------------------------------------------------------------
! subroutine rand_permute( )

! Generates a vector of length n
!    containing the integers 0, ... , n-1 in random order.
! We do not use a new seed.
! Follow `Quasi-Monte Carlo sampling' by Art Owen, Section 1.3
! He follows, in turn, Luc Devroye 'Non-uniform random ...' (1986)

! Input: n = number of elements to permute

! Output: pvect = array of n numbers in random order
!----------------------------------------------------------------------
        subroutine rand_permute( n, pvect )

        implicit none

! Input
  
        integer, intent(in) :: n

! Output

        integer, intent(out) :: pvect(1:n)

! Local

        integer j, k, temp, seed

! Continue the old string of random numbers by choosing seed>0
        seed = 1;

! Start with an ordered vector, pvect
        do j=1,n
          pvect(j) = j
        enddo

! Now re-arrange the elements
        do j=n,2,-1
          temp = pvect(j)
          ! choose an element randomly between 1 and j
          k = int(dble(j)*ran2(seed)+1.0);
          ! swap elements j and k 
          pvect(j) = pvect(k)
          pvect(k) = temp
        enddo

! Convert range of array from 1:n to 0:n-1
        do j=1,n
          pvect(j) = pvect(j) - 1
        enddo

        return
        end subroutine rand_permute
!------------------------------------------------------------------------
        double precision FUNCTION ran2(idum)

        implicit none
        
        ! Input Variable(s)
        INTEGER, intent(inout) :: idum

        integer :: IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        double precision AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *  NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.-EPS)
        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2
!------------------------------------------------------------------------
      end module random
