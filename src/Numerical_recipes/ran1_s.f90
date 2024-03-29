! $Id: num_rec.f90,v 1.1 2008-07-24 17:31:45 dschanen Exp $
!   From _Numerical Recipes in Fortran 90_
!   (C) 1988-1996 Numerical Recipes Software
SUBROUTINE ran1_s( harvest )
  
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
    iran0,jran0,kran0,nran0,mran0,rans

IMPLICIT NONE

REAL(SP), INTENT(OUT) :: harvest

if (lenran < 1) call ran_init(1)
rans = iran0-kran0
if (rans < 0) rans=rans+2147483579_k4b
iran0 = jran0
jran0 = kran0
kran0 = rans
nran0 = ieor(nran0,ishft(nran0,13))
nran0 = ieor(nran0,ishft(nran0,-17))
nran0 = ieor(nran0,ishft(nran0,5))
if (nran0 == 1) nran0=270369_k4b
mran0 = ieor(mran0,ishft(mran0,5))
mran0 = ieor(mran0,ishft(mran0,-13))
mran0 = ieor(mran0,ishft(mran0,6))
rans  = ieor(nran0,rans)+mran0

harvest = amm*merge(rans,not(rans), rans<0 )

return
END SUBROUTINE ran1_s
