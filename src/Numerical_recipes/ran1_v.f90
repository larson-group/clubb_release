
SUBROUTINE ran1_v(harvest)
USE nrtype
USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
  iran,jran,kran,nran,mran,ranv
IMPLICIT NONE

REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
INTEGER(K4B) :: n

n=size(harvest)
if (lenran < n+1) call ran_init(n+1)
ranv(1:n) = iran(1:n)-kran(1:n)
where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
iran(1:n) = jran(1:n)
jran(1:n) = kran(1:n)
kran(1:n) = ranv(1:n)
nran(1:n) = ieor(nran(1:n),ishft(nran(1:n),13))
nran(1:n) = ieor(nran(1:n),ishft(nran(1:n),-17))
nran(1:n) = ieor(nran(1:n),ishft(nran(1:n),5))
where (nran(1:n) == 1) nran(1:n)=270369_k4b
mran(1:n) = ieor(mran(1:n), ishft(mran(1:n),5))
mran(1:n) = ieor(mran(1:n), ishft(mran(1:n),-13))
mran(1:n) = ieor(mran(1:n), ishft(mran(1:n),6))
ranv(1:n) = ieor(nran(1:n), ranv(1:n))+mran(1:n)

harvest = amm * merge( ranv(1:n), not(ranv(1:n)), ranv(1:n)<0 )

return
END SUBROUTINE ran1_v
