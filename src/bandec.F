      SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d,isValid )
      INTEGER m1,m2,mp,mpl,n,np,indx(n)
      REAL d,a(np,mp),al(np,mpl),TINY
      PARAMETER (TINY=1.e-20)
      INTEGER i,j,k,l,mm
      REAL dum
      LOGICAL, INTENT(INOUT) :: isValid
      mm=m1+m2+1
      if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) stop 'bad args in bandec'
      l=m1
      do 13 i=1,m1
        do 11 j=m1+2-i,mm
          a(i,j-l)=a(i,j)
11      continue
        l=l-1
        do 12 j=mm-l,mm
          a(i,j)=0.
12      continue
13    continue
      d=1.
      l=m1
      do 18 k=1,n
        dum=a(k,1)
        i=k
        if(l.lt.n)l=l+1
        do 14 j=k+1,l
          if(abs(a(j,1)).gt.abs(dum))then
            dum=a(j,1)
            i=j
          endif
14      continue
        indx(k)=i
c-->Chris Golaz modification
c-->Modified by dschanen to return with isValid=.false. rather than stoping
        if(dum.eq.0.) a(k,1)=TINY
        if(dum.eq.0.) then
          write(*,*) 'Singular matrix encountered'
          isValid = .false.
          return
        endif
c<--end of modification
        if(i.ne.k)then
          d=-d
          do 15 j=1,mm
            dum=a(k,j)
            a(k,j)=a(i,j)
            a(i,j)=dum
15        continue
        endif
        do 17 i=k+1,l
          dum=a(i,1)/a(k,1)
          al(k,i-k)=dum
          do 16 j=2,mm
            a(i,j-1)=a(i,j)-dum*a(k,j)
16        continue
          a(i,mm)=0.
17      continue
18    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
