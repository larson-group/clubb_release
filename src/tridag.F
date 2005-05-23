! 
! Solves for a vector u(1:n) of length n the tridiagonal linear set given
! by equation (2.4.1). a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors
! and are not modified.
! Parameter: NMAX is the maximum expected value of n.
!
! Reference: Numerical Recipes in Fortran, The Art of Scientific Computing
! 2nd Edtion, William H. Press et al., 1992
!

        subroutine tridag(a,b,c,rr,u,n, isValid)
        implicit none

        integer n
        real a(n),b(n),c(n),rr(n),u(n)

        integer j
        real bet
        real, allocatable :: gam(:)
        logical, intent(inout) :: isValid

!        Allocate memory

        allocate( gam(1:n) )

!        Solve system
        if (b(1).eq.0.) then
           do j=1,n
              write(*,*) j,b(j)
           end do
           isValid = .false.
           print *, "tridag: rewrite equations"
           deallocate( gam )
           return
!           stop 'tridag: rewrite equations'
        end if
        bet=b(1)
        u(1)=rr(1)/bet
        do j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
!          if (bet.eq.0.) stop 'tridag failed'
!         Return and handle the exception
          if (bet.eq.0.) then
             isValid = .false.
             print *, "tridag failed"
             deallocate( gam )
             return
          endif
          u(j)=(rr(j)-a(j)*u(j-1))/bet
        enddo
        do j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
        enddo

!        Deallocate memory

        deallocate( gam )
       
        return
        end subroutine tridag
