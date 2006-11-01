!------------------------------------------------------------------------------
! $Id: gamma_polpak.f90,v 1.1 2006-11-01 20:42:56 mjfalk Exp $

! Description:
! POLPAK is a library of FORTRAN90 routines, using double precision arithemetic, that evaluate a variety of mathematical functions.
!
! It includes routines to evaluate the recursively defined polynomial families of
!
!   * Bernoulli
!   * Bernstein
!   * Cardan
!   * Chebyshev
!   * Euler
!   * Gegenbauer
!   * Hermite
!   * Jacobi
!   * Laguerre
!   * Legendre
!   * Zernike

! This code includes only a modified version of the POLPAK approximation of the 
! Gamma function, and was put into this module by us.

! Distributed by John Burkardt.

! References:
! <http://www.csit.fsu.edu/~burkardt/> 

!------------------------------------------------------------------------------

module polpak_gamma
implicit none

public :: gamma, gamma_log

contains

function gamma( x )

!*******************************************************************************
!
!! GAMMA returns the value of the Gamma function at X.
!
!  Definition:
!
!    GAMMA(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) EXP(-T) dT
!
!  Recursion:
!
!    GAMMA(X+1) = X*GAMMA(X)
!
!  Restrictions:
!
!    0 < X ( a software restriction).
!    26 May 2005 worked around this restriction by incorporating and
!      updating the an older gamma approx. code into this subroutine 
!    -David Schanen
!
!  Special values:
!
!    GAMMA(0.5) = sqrt(PI)
!
!    For N a positive integer, GAMMA(N+1) = N!, the standard factorial.
!
!  Modified:
!    
!    14 April 1999 Burkardt
!    26 May 2005 Schanen
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which the Gamma function
!    is desired.
!
!    Output, real ( kind = 8 ) GAMMA, the Gamma function of X.
!
  implicit none

  real ( kind = 8 ) gamma
! real ( kind = 8 ) gamma_log
  real ( kind = 8 ) x
!--- negative x addition ---!
  real( kind = 8 ) :: g
  real( kind = 8 ) :: z, r
  integer          :: i, imax

  real( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real( kind = 8 ),  dimension(26), parameter :: data_g = (/ &
    1.0D0, 0.5772156649015329D0, &
    -0.6558780715202538D0, -0.420026350340952D-1, &
    0.1665386113822915D0,-.421977345555443D-1, &
    -.96219715278770D-2, .72189432466630D-2, &
    -.11651675918591D-2, -.2152416741149D-3, &
    0.1280502823882D-3, -.201348547807D-4, &
    -.12504934821D-5, .11330272320D-5, &
    -.2056338417D-6, .61160950D-8, &
    0.50020075D-8, -.11812746D-8, &
    0.1043427D-9, .77823D-11, &
    -.36968D-11, .51D-12, &
    -.206D-13, -.54D-14, .14D-14, .1D-15/)
!--- end negative x addition ---!

  if ( x >= 0. ) then
    gamma = exp ( gamma_log ( x ) ) 
  else
!--- negative x addition ---!
    if ( abs( x ) > 1.0D0 ) then
      z = abs( x )
      imax = int( z )
      r = 1.0D0
      do i=1, imax
        r = r * ( z - real(i) )
      enddo
      z = z - real(imax)
    else
      z = x
    endif
  
    g = data_g(26)
    do i = 25, 1, -1
      g = g * z + data_g(i)
    enddo
    gamma = 1.0D0 / (g * z)
    if ( abs( x ) > 1.0D0 ) then
      g = g * r
      if ( x < 0.0D0 ) g = -pi / (x * g * sin(pi * x) )
    endif
!--- end negative x addition ---!
  endif

  return
end function gamma

function gamma_log ( x )

!*******************************************************************************
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    log ( GAMMA(X) ) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from reference 3, while approximations
!    for X < 12.0 are similar to those in reference 1, but are unpublished.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    16 June 1999
!
!  Authors:
!
!    W. J. Cody and L. Stoltz
!    Argonne National Laboratory
!
!  Reference:
!
!    # 1)
!    W. J. Cody and K. E. Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    # 2)
!    K. E. Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    # 3)
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function. 
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma 
!    function of X.  If X <= 0.0, or if overflow would occur, the
!    program returns the value HUGE().
!
!  Machine-dependent constants:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!             in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    XINF   - largest machine representable floating-point number;
!             approximately BETA**MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!    Approximate values for some important machines are:
!
!                              BETA      MAXEXP         XBIG
!
!    CRAY-1        (S.P.)        2        8191       9.62D+2461
!    Cyber 180/855
!      under NOS   (S.P.)        2        1070       1.72D+319
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)        2         128       4.08D+36
!    IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!    IBM 3033      (D.P.)       16          63       4.29D+73
!    VAX D-Format  (D.P.)        2         127       2.05D+36
!    VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            FRTBIG
!
!    CRAY-1        (S.P.)   3.13D+615
!    Cyber 180/855
!      under NOS   (S.P.)   6.44D+79
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)   1.42D+9
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)   2.25D+76
!    IBM 3033      (D.P.)   2.56D+18
!    VAX D-Format  (D.P.)   1.20D+9
!    VAX G-Format  (D.P.)   1.89D+76
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 = - 5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  integer i
  real ( kind = 8 ) gamma_log
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  eps = epsilon ( eps )

  if ( x <= eps ) then

    res = - log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = - 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  gamma_log = res

  return
end function gamma_log

end module polpak_gamma
