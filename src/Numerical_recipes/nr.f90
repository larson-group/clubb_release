! $Id: num_rec.f90,v 1.1 2008-07-24 17:31:45 dschanen Exp $
!   From _Numerical Recipes in Fortran 90_
!   (C) 1988-1996 Numerical Recipes Software

! This is an incomplete version containing only the interfaces needed to make
! use of the amoeba and amebsa routines in CLUBB. -dschanen 1 Apr 2010

MODULE nr

  INTERFACE
    SUBROUTINE amebsa(p,y,pb,yb,f_tol,func,iter,temptr)
    USE nrtype
    INTEGER(I4B), INTENT(INOUT) :: iter
    REAL(SP), INTENT(INOUT) :: yb
    REAL(SP), INTENT(IN) :: f_tol,temptr
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: y,pb
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
    INTERFACE
      FUNCTION func(x)
      USE nrtype
      REAL(SP), DIMENSION(:), INTENT(IN) :: x
      REAL(SP) :: func
      END FUNCTION func
    END INTERFACE
    END SUBROUTINE amebsa
  END INTERFACE

  INTERFACE
    SUBROUTINE amoeba(p,y,f_tol,func,iter)
    USE nrtype
    INTEGER(I4B), INTENT(OUT) :: iter
    REAL(SP), INTENT(IN) :: f_tol
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p
    INTERFACE
      FUNCTION func(x)
      USE nrtype
      REAL(SP), DIMENSION(:), INTENT(IN) :: x
      REAL(SP) :: func
      END FUNCTION func
    END INTERFACE
    END SUBROUTINE amoeba
  END INTERFACE

  INTERFACE ran1
    SUBROUTINE ran1_s(harvest)
    USE nrtype
    REAL(SP), INTENT(OUT) :: harvest
    END SUBROUTINE ran1_s

    SUBROUTINE ran1_v(harvest)
    USE nrtype
    REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
    END SUBROUTINE ran1_v
  END INTERFACE

END MODULE nr
