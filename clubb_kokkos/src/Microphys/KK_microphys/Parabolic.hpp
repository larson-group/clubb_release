// Native Algorithm850 special functions used by the KK engine closure.
// Retains the dispatcher, U recurrence, McLaurin/elementary/Poincare algorithms
// and all original coefficients/stopping tests for a>=1.5 and mode=0. Every
// supported KK evaporation/autoconversion/accretion PDF call has this order;
// arbitrary negative-a/Airy/scaled direct helper calls are outside this closure.
// References: ACM TOMS32(1),2006,102–112; Native Parabolic.f90.
#pragma once
#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <limits>
#include <array>
#include <stdexcept>
namespace clubb_kokkos::parabolic {
KOKKOS_INLINE_FUNCTION double gamma(double x) {
  // Native Parabolic gamma: Cody's ACM TOMS Algorithm 715 approximation.
  // The public mixed-moment factorial uses this floating-point result rather
  // than an integer factorial, retaining Native coefficients and evaluation.
  constexpr double one = 1., half = .5, twelve = 12., two = 2., zero = 0.;
  constexpr double sqrtpi = .9189385332046727417803297;
  constexpr double pi     = 3.1415926535897932384626434;
  constexpr double xbig   = 171.624;
  constexpr double xminin = std::numeric_limits<double>::min();
  constexpr double eps    = std::numeric_limits<double>::epsilon();
  constexpr double xinf   = std::numeric_limits<double>::max();
  const double p[8]       = {-1.71618513886549492533811e+0, 2.47656508055759199108314e+1,
                             -3.79804256470945635097577e+2, 6.29331155312818442661052e+2,
                             8.66966202790413211295064e+2,  -3.14512729688483675254357e+4,
                             -3.61444134186911729807069e+4, 6.64561438202405440627855e+4};
  const double q[8]       = {-3.08402300119738975254353e+1, 3.15350626979604161529144e+2,
                             -1.01515636749021914166146e+3, -3.10777167157231109440444e+3,
                             2.25381184209801510330112e+4,  4.75584627752788110767815e+3,
                             -1.34659959864969306392456e+5, -1.15132259675553483497211e+5};
  const double c[7]       = {
      -1.910444077728e-03,     8.4171387781295e-04,          -5.952379913043012e-04,
      7.93650793500350248e-04, -2.777777777777681622553e-03, 8.333333333333333331554247e-02,
      5.7083835261e-03};
  bool parity = false;
  double fact = one;
  int n       = 0;
  double y    = x;
  int k       = 0;
  double y1, z, fn_val, xnum, xden, ysq, sum;
  if (y <= zero) {
    y      = -x;
    y1     = Kokkos::trunc(y);
    fn_val = y - y1;
    if (fn_val != zero) {
      if (y1 != Kokkos::trunc(y1 * half) * two)
        parity = true;
      fact = -pi / Kokkos::sin(pi * fn_val);
      y    = y + one;
    } else {
      fn_val = xinf;
      k      = 1;
    }
  }
  if (k == 0) {
    if (y < eps) {
      if (y >= xminin)
        fn_val = one / y;
      else {
        fn_val = xinf;
        k      = 1;
      }
    } else if (y < twelve) {
      y1 = y;
      if (y < one) {
        z = y;
        y = y + one;
      } else {
        n = static_cast<int>(y) - 1;
        y = y - n;
        z = y - one;
      }
      xnum = zero;
      xden = one;
      for (int i = 0; i < 8; ++i) {
        xnum = (xnum + p[i]) * z;
        xden = xden * z + q[i];
      }
      fn_val = xnum / xden + one;
      if (y1 < y)
        fn_val = fn_val / y1;
      else if (y1 > y)
        for (int i = 1; i < n + 1; ++i) {
          fn_val = fn_val * y;
          y      = y + one;
        }
    } else {
      if (y <= xbig) {
        ysq = y * y;
        sum = c[6];
        for (int i = 0; i < 6; ++i)
          sum = sum / ysq + c[i];
        sum    = sum / y - y + sqrtpi;
        sum    = sum + (y - half) * Kokkos::log(y);
        fn_val = Kokkos::exp(sum);
      } else {
        fn_val = xinf;
        k      = 1;
      }
    }
  }
  if (k == 0) {
    if (parity)
      fn_val = -fn_val;
    if (fact != one)
      fn_val = fact / fn_val;
  }
  return fn_val;
}

// Fixed local arrays use zero-based C++ slots and hold no case/model state.
// Algorithm850's physical recurrence orders and counters retain their values.
template <int N> struct ParabolicArray {
  std::array<double, N> values{};
  double& operator()(int i) { return values[i]; }
  const double& operator()(int i) const { return values[i]; }
};
inline constexpr double lnsqrttwopi   = 0.9189385332046727418;
inline constexpr double sqrttwopi     = 2.5066282746310005024;
inline constexpr double eulmasc       = 0.5772156649015328606;
inline constexpr double pihalf        = 1.5707963267948966192;
inline constexpr double pikwart       = 0.7853981633974483096;
inline constexpr double pidriekwart   = 2.3561944901923449288;
inline constexpr double sqrtpi        = 1.7724538509055160272;
inline constexpr double lnpi          = 1.1447298858494001741;
inline constexpr double pi            = 3.1415926535897932385;
inline constexpr double onepi         = .31830988618379067153;
inline constexpr double sqrt2         = 1.4142135623730950488;
inline constexpr double sqrt3         = 1.7320508075688772935;
inline constexpr double sqrt2opi      = 0.7978845608028653559;
inline constexpr double twopi         = 6.2831853071795864769;
inline constexpr double oneoversqrtpi = 0.5641895835477562869;
inline constexpr double twooversqrtpi = 1.1283791670955125739;
inline constexpr double onethird      = 0.33333333333333333;
inline constexpr double twothird      = 0.6666666666666666666667;
inline constexpr double onesix        = .166666666666666666666666666667;
inline constexpr double piquart       = 0.78539816339744830962;
inline constexpr double twoexp14      = 1.18920711500272106671749997;
inline constexpr double dwarf         = std::numeric_limits<double>::min() * 1000.;
inline constexpr double giant         = std::numeric_limits<double>::max() / 1000.;
inline constexpr double mac_tol       = std::numeric_limits<double>::epsilon();
inline void parab(double a, double x, int mode, ParabolicArray<2>& uaxx, ParabolicArray<2>& vaxx,
                  int& ierr, double epss);
inline void recuru(double a, double x, int mode, ParabolicArray<2>& uaxx, int& ierr, double epss);
inline void series(double a, double x, int mode, ParabolicArray<2>& uax, ParabolicArray<2>& vax,
                   int& ierr, double epss);
inline void expaelem(double a, double x, int mode, ParabolicArray<2>& uaxx, ParabolicArray<2>& vaxx,
                     int& ierr, double epss);
inline void expax(double a, double x, int mode, ParabolicArray<2>& uax, ParabolicArray<2>& vax,
                  int& ierr, double epss);
inline void parab(double a, double x, int mode, ParabolicArray<2>& uaxx, ParabolicArray<2>& vaxx,
                  int& ierr, double epss) {
  // Evaluate the Native parab algorithm in call-local storage.
  if (a < 1.5 || mode != 0)
    throw std::runtime_error(
        "Algorithm850 direct request is outside the supported KK call closure");
  ParabolicArray<2> uaxx1;
  ParabolicArray<2> vaxx1;
  double x2;
  double f1;
  double f2;
  double f3;
  double f4;
  double f5;
  double f6;
  double f7;
  double f8;
  double f9;
  double f10;
  ierr = 0;
  if ((x < 0)) {
    uaxx(0) = 0.0;
    uaxx(1) = 0.0;
    vaxx(0) = 0.0;
    vaxx(1) = 0.0;
    ierr    = 2;
  }
  if ((ierr == 0)) {
    if ((x < .1857815261497950467)) {
      x2 = (x * x);
      f1 = ((((-0.23) * x2) + (1.2 * x)) + 18.72);
      if ((a > f1)) {
        if (((a < (1.0 / x2)) && (x < 0.005))) {
          series(a, x, mode, uaxx1, vaxx, ierr, epss);
          expaelem(a, x, mode, uaxx, vaxx1, ierr, epss);
        } else {
          expaelem(a, x, mode, uaxx, vaxx, ierr, epss);
        }
      } else {
        series(a, x, mode, uaxx, vaxx, ierr, epss);
      }
    } else {
      if ((x < 3.0)) {
        f1 = (((((-0.23) * x) * x) + (1.2 * x)) + 18.72);
        f2 = ((3.75 / x) - 1.25);
        if ((a < f2)) {
          series(a, x, mode, uaxx, vaxx, ierr, epss);
        } else {
          if ((a < f1)) {
            series(a, x, mode, uaxx, vaxx, ierr, epss);
            recuru(a, x, mode, uaxx, ierr, epss);
          } else {
            expaelem(a, x, mode, uaxx, vaxx, ierr, epss);
          }
        }
      } else {
        if ((x < 12.0)) {
          f1 = (((((-0.23) * x) * x) + (1.2 * x)) + 18.72);
          if ((a < f1)) {
            series(a, x, mode, uaxx, vaxx, ierr, epss);
            recuru(a, x, mode, uaxx, ierr, epss);
          } else {
            expaelem(a, x, mode, uaxx, vaxx, ierr, epss);
          }
        } else {
          f7 = ((2.5 * x) - 30.0);
          if (((a < f7) && (a < 150))) {
            expax(a, x, mode, uaxx, vaxx, ierr, epss);
          } else {
            expaelem(a, x, mode, uaxx, vaxx, ierr, epss);
          }
        }
      }
    }
  }
}
inline void recuru(double a, double x, int mode, ParabolicArray<2>& uaxx, int& ierr, double epss) {
  // Evaluate the Native recuru algorithm in call-local storage.
  if (a < 1.5 || mode != 0)
    throw std::runtime_error(
        "Algorithm850 helper request is outside the supported KK call closure");
  double eps;
  ParabolicArray<2> uaxx1;
  ParabolicArray<2> vaxx1;
  ParabolicArray<2> uaxx2;
  ParabolicArray<2> vaxx2;
  double a1;
  double a2;
  double a3;
  double xargu;
  double xhalf;
  double gammaa;
  int ierr1;
  int ierr2;
  int mode2;
  eps   = epss;
  ierr  = 0;
  mode2 = 0;
  xhalf = (0.5 * x);
  a1    = ((-(static_cast<int>(a) - a)) + 21);
  a2    = (1 + a1);
  parab(a1, x, mode2, uaxx1, vaxx1, ierr1, epss);
  parab(a2, x, mode2, uaxx2, vaxx2, ierr2, epss);
  a3      = (a1 - 1);
  uaxx(0) = ((x * uaxx1(0)) + ((a3 + 1.5) * uaxx2(0)));
  uaxx(1) = (((-xhalf) * uaxx(0)) - ((a3 + 0.5) * uaxx1(0)));
  while ((Kokkos::abs((a3 - a)) > (eps * 100))) {
    uaxx2(0) = uaxx1(0);
    uaxx1(0) = uaxx(0);
    a2       = a1;
    a1       = a3;
    a3       = (a1 - 1);
    uaxx(0)  = ((x * uaxx1(0)) + ((a3 + 1.5) * uaxx2(0)));
    uaxx(1)  = (((-xhalf) * uaxx(0)) - ((a3 + 0.5) * uaxx1(0)));
  }
}
inline void series(double a, double x, int mode, ParabolicArray<2>& uax, ParabolicArray<2>& vax,
                   int& ierr, double epss) {
  // Evaluate the Native series algorithm in call-local storage.
  if (a < 1.5 || mode != 0)
    throw std::runtime_error(
        "Algorithm850 helper request is outside the supported KK call closure");
  double eps;
  double ah;
  double ahp14;
  double ahp34;
  double zerk;
  double f1;
  double f2;
  double err1;
  double err2;
  double x2;
  double am1;
  double a0;
  double bm1;
  double b0;
  double y1;
  double y2;
  double x2k;
  double x2km;
  double x2kd;
  double facto1;
  double facto2;
  double a2;
  double b2;
  double acof;
  double bcof;
  double twok;
  double argu;
  double y1p;
  double y2p;
  double facto1d;
  double facto2d;
  double acofd;
  double bcofd;
  double err1p;
  double err2p;
  double t;
  double t2;
  double ta;
  double xargu;
  double gammaa;
  double acofn;
  double bcofn;
  double acofdn;
  double bcofdn;
  double mubA;
  double mu2A;
  double mu4A;
  double mubB;
  double mu2B;
  double mu4B;
  double musA;
  double gmusA;
  double musB;
  double gmusB;
  double smusA;
  double smusB;
  double betaA;
  double betaB;
  double sa;
  ParabolicArray<9> coefs;
  double yk;
  double y;
  double ha;
  double dd;
  ParabolicArray<26> gk;
  ParabolicArray<13> wk;
  int k;
  int i;
  int li;
  gk(0)  = 1.0;
  gk(1)  = .41666666666666666667e-1;
  gk(2)  = 0;
  gk(3)  = (-.97463348765432098765e-2);
  gk(4)  = 0;
  gk(5)  = .12318482904510826965e-1;
  gk(6)  = 0;
  gk(7)  = (-.37822933917705539682e-1);
  gk(8)  = 0;
  gk(9)  = .21514326767209896474;
  gk(10) = 0;
  gk(11) = (-1.9630003732872175294);
  gk(12) = 0;
  gk(13) = 26.256830962378916652;
  gk(14) = 0;
  gk(15) = (-484.19061617504532506);
  gk(16) = 0;
  gk(17) = 11773.948564802034554;
  gk(18) = 0;
  gk(19) = (-365037.92569092371983);
  gk(20) = 0;
  gk(21) = 14054558.808383655048;
  gk(22) = 0;
  gk(23) = (-657894020.31009296820);
  gk(24) = 0;
  gk(25) = 36795321248.737494263;
  wk(0)  = 1.0;
  wk(1)  = (-.17361111111111111111e-2);
  wk(2)  = .81219457304526748971e-3;
  wk(3)  = (-.11215312855682914598e-2);
  wk(4)  = .33920312789254653178e-2;
  wk(5)  = (-.18817620620360951108e-1);
  wk(6)  = .16870892343666095644;
  wk(7)  = (-2.2330644165264046334);
  wk(8)  = 40.925670821571152909;
  wk(9)  = (-991.44221631790969480);
  wk(10) = 30664.092693023833555;
  wk(11) = (-1178670.6504836921080);
  wk(12) = 55108658.068376820265;
  eps    = epss;
  ierr   = 0;
  ah     = (a * 0.5);
  ahp14  = (0.25 + ah);
  ahp34  = (0.5 + ahp14);
  err1   = 1.0;
  err2   = 1.0;
  err1p  = 1.0;
  err2p  = 1.0;
  x2     = (x * x);
  am1    = 1.0;
  a0     = a;
  bm1    = 1.0;
  b0     = a;
  y1     = (am1 + ((a0 * x2) * 0.5));
  y1p    = (a0 * x);
  y2     = ((bm1 + ((b0 * x2) / 6.0)) * x);
  y2p    = y1;
  k      = 2;
  x2k    = x2;
  x2km   = (x2 * x);
  x2kd   = x;
  facto1 = 0.5;
  facto2 = (1.0 / 6.0);
  while (((((err1 > eps) || (err2 > eps)) || (err1p > eps)) || (err2p > eps))) {
    acofn  = 0;
    bcofn  = 0;
    acofdn = 0;
    bcofdn = 0;
    for (li = 1; li < 3; li += 1) {
      twok    = (2.0 * k);
      zerk    = (0.5 * (k - 1.0));
      a2      = ((a * a0) + ((zerk * (twok - 3.0)) * am1));
      b2      = ((a * b0) + ((zerk * (twok - 1.0)) * bm1));
      x2k     = (x2k * x2);
      x2km    = (x2km * x2);
      x2kd    = (x2kd * x2);
      facto1d = (facto1 / (twok - 1.0));
      facto2d = (facto2 / twok);
      facto1  = (facto1d / twok);
      facto2  = (facto2d / (twok + 1.0));
      acof    = (a2 * (x2k * facto1));
      bcof    = (b2 * (x2km * facto2));
      acofd   = (a2 * (x2kd * facto1d));
      bcofd   = (b2 * (x2k * facto2d));
      y1      = (y1 + acof);
      y2      = (y2 + bcof);
      y1p     = (y1p + acofd);
      y2p     = (y2p + bcofd);
      k       = ((k + 1));
      am1     = a0;
      a0      = a2;
      bm1     = b0;
      b0      = b2;
      acofn   = (acofn + acof);
      bcofn   = (bcofn + bcof);
      acofdn  = (acofdn + acofd);
      bcofdn  = (bcofdn + bcofd);
    }
    err1  = Kokkos::abs((acofn / ((y1) > (eps) ? (y1) : (eps))));
    err2  = Kokkos::abs((bcofn / ((y2) > (eps) ? (y2) : (eps))));
    err1p = Kokkos::abs((acofdn / ((y1p) > (eps) ? (y1p) : (eps))));
    err2p = Kokkos::abs((bcofdn / ((y2p) > (eps) ? (y2p) : (eps))));
  }
  f1     = (sqrtpi / (Kokkos::pow(2.0, ahp14) * gamma(ahp34)));
  f2     = ((-sqrtpi) / (Kokkos::pow(2.0, (ah - 0.25)) * gamma(ahp14)));
  uax(0) = ((f1 * y1) + (f2 * y2));
  uax(1) = ((f1 * y1p) + (f2 * y2p));
  f1 = ((((Kokkos::pow(2.0, ahp14) / pi) * gamma(ahp14)) * Kokkos::sin((pi * ((3.0 / 4.0) - ah)))) *
        Kokkos::sin((pi * ((3.0 / 4.0) - ah))));
  f2 = ((((Kokkos::pow(2.0, (ahp14 + 0.5)) / pi) * gamma((ah + (3.0 / 4.0)))) *
         Kokkos::sin((pi * (0.25 - ah)))) *
        Kokkos::sin((pi * (0.25 - ah))));
  vax(0) = ((f1 * y1) + (f2 * y2));
  vax(1) = ((f1 * y1p) + (f2 * y2p));
}
inline void expaelem(double a, double x, int mode, ParabolicArray<2>& uaxx, ParabolicArray<2>& vaxx,
                     int& ierr, double epss) {
  // Evaluate the Native expaelem algorithm in call-local storage.
  if (a < 1.5 || mode != 0)
    throw std::runtime_error(
        "Algorithm850 helper request is outside the supported KK call closure");
  double eps;
  double mu;
  double mu2;
  double t;
  double t2p1;
  double facu;
  double t2;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  double t10;
  double t11;
  double t12;
  double t13;
  double t14;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t25;
  double t26;
  double t27;
  double t28;
  double t29;
  double t30;
  double t31;
  double t32;
  double t33;
  double t34;
  double t35;
  double t36;
  double t37;
  double t38;
  double t39;
  double t40;
  double t41;
  double tau;
  double tau2;
  double tau3;
  double tau4;
  double tau5;
  double tau6;
  double tau7;
  double tau8;
  double tau9;
  double tau10;
  double tau11;
  double tau12;
  double tau13;
  double tau14;
  double tau15;
  double tau16;
  double tau17;
  double tau18;
  double tau19;
  double tau20;
  double tau21;
  double tau22;
  double tau23;
  double tau24;
  double tau25;
  double tau26;
  double tau27;
  double tau28;
  double tau29;
  double tau30;
  double tau31;
  double tau32;
  double tau33;
  double tau34;
  double tau35;
  double tau36;
  double ffa;
  double x2;
  double xargu;
  double sinpia;
  double smus;
  double musp;
  double musn;
  double sphip;
  double schip;
  double sphin;
  double schin;
  double gammaa2;
  double e1;
  double e2;
  double e3;
  double e4;
  double delta1;
  double delta2;
  double delta3;
  double delta4;
  double maxi;
  double mu4;
  double tts;
  double ttss;
  double eta;
  double mus;
  double gmus;
  double s1;
  double s2;
  double s3;
  double s4;
  double mus1;
  double mus2;
  double argu3;
  double dl;
  double dnew;
  double beta;
  ParabolicArray<16> phi;
  ParabolicArray<16> chi;
  ParabolicArray<16> u;
  ParabolicArray<16> v;
  ParabolicArray<26> gk;
  ParabolicArray<13> wk;
  int k;
  int itwok;
  gk(0)  = 1.0;
  gk(1)  = .41666666666666666667e-1;
  gk(2)  = 0;
  gk(3)  = (-.97463348765432098765e-2);
  gk(4)  = 0;
  gk(5)  = .12318482904510826965e-1;
  gk(6)  = 0;
  gk(7)  = (-.37822933917705539682e-1);
  gk(8)  = 0;
  gk(9)  = .21514326767209896474;
  gk(10) = 0;
  gk(11) = (-1.9630003732872175294);
  gk(12) = 0;
  gk(13) = 26.256830962378916652;
  gk(14) = 0;
  gk(15) = (-484.19061617504532506);
  gk(16) = 0;
  gk(17) = 11773.948564802034554;
  gk(18) = 0;
  gk(19) = (-365037.92569092371983);
  gk(20) = 0;
  gk(21) = 14054558.808383655048;
  gk(22) = 0;
  gk(23) = (-657894020.31009296820);
  gk(24) = 0;
  gk(25) = 36795321248.737494263;
  wk(0)  = 1.0;
  wk(1)  = (-.17361111111111111111e-2);
  wk(2)  = .81219457304526748971e-3;
  wk(3)  = (-.11215312855682914598e-2);
  wk(4)  = .33920312789254653178e-2;
  wk(5)  = (-.18817620620360951108e-1);
  wk(6)  = .16870892343666095644;
  wk(7)  = (-2.2330644165264046334);
  wk(8)  = 40.925670821571152909;
  wk(9)  = (-991.44221631790969480);
  wk(10) = 30664.092693023833555;
  wk(11) = (-1178670.6504836921080);
  wk(12) = 55108658.068376820265;
  eps    = epss;
  ierr   = 0;
  mu     = Kokkos::sqrt((2.0 * Kokkos::abs(a)));
  t      = (Kokkos::abs(x) / (mu * sqrt2));
  xargu  = Kokkos::abs((((x * x) * 0.25) + a));
  if (((((a * Kokkos::log(((x * 0.5) + Kokkos::sqrt(xargu)))) + ((x * 0.5) * Kokkos::sqrt(xargu))) -
        (a * 0.5)) > Kokkos::log(giant))) {
    ierr    = 1;
    uaxx(0) = dwarf;
    uaxx(1) = dwarf;
    vaxx(0) = giant;
    vaxx(1) = giant;
  }
  if ((ierr == 0)) {
    x2      = (x * x);
    tau     = (((-a) * 0.5) / (((x2 * 0.25) + a) + ((x * 0.5) * Kokkos::sqrt(((x2 * 0.25) + a)))));
    tau2    = (tau * tau);
    tau3    = (tau2 * tau);
    tau4    = (tau3 * tau);
    tau5    = (tau4 * tau);
    tau6    = (tau5 * tau);
    tau7    = (tau6 * tau);
    tau8    = (tau7 * tau);
    tau9    = (tau8 * tau);
    tau10   = (tau9 * tau);
    tau11   = (tau10 * tau);
    tau12   = (tau11 * tau);
    tau13   = (tau12 * tau);
    tau14   = (tau13 * tau);
    tau15   = (tau14 * tau);
    tau16   = (tau15 * tau);
    tau17   = (tau16 * tau);
    tau18   = (tau17 * tau);
    tau19   = (tau18 * tau);
    tau20   = (tau19 * tau);
    tau21   = (tau20 * tau);
    tau22   = (tau21 * tau);
    tau23   = (tau22 * tau);
    tau24   = (tau23 * tau);
    tau25   = (tau24 * tau);
    tau26   = (tau25 * tau);
    tau27   = (tau26 * tau);
    tau28   = (tau27 * tau);
    tau29   = (tau28 * tau);
    tau30   = (tau29 * tau);
    tau31   = (tau30 * tau);
    tau32   = (tau31 * tau);
    tau33   = (tau32 * tau);
    tau34   = (tau33 * tau);
    tau35   = (tau34 * tau);
    tau36   = (tau35 * tau);
    phi(0)  = 1.0;
    phi(1)  = ((((-1.666666666666666666666667) * tau3) - (2.5 * tau2)) - (.75 * tau));
    // BFB-031: Native contracts the low/high power groups separately before
    // adding them. Keep the same phi(2) polynomial and its original coefficients.
    const double phi_2_low_powers =
        Kokkos::fma(3.28125, tau2, Kokkos::fma(27.875, tau3, 67.375 * tau4));
    const double phi_2_high_powers =
        Kokkos::fma(64.16666666666666666666667, tau5, 21.38888888888888888888889 * tau6);
    phi(2) = phi_2_low_powers + phi_2_high_powers;
    phi(3)  = ((((((((-525.2160493827160493827160) * tau9) - (2363.472222222222222222222 * tau8)) -
                   (4254.25 * tau7)) -
                  (3861.229166666666666666667 * tau6)) -
                 (1814.5125 * tau5)) -
                (396.328125 * tau4)) -
               (27.0703125 * tau3));
    phi(4)  = (((((((((19126.61779835390946502058 * tau12) + (114759.7067901234567901235 * tau11)) +
                     (292637.2523148148148148148 * tau10)) +
                    (411244.1666666666666666667 * tau9)) +
                   (344694.4244791666666666667 * tau8)) +
                  (173584.125 * tau7)) +
                 (49877.23984375 * tau6)) +
                (7077.45703125 * tau5)) +
               (329.91943359375 * tau4));
    phi(5)  = ((((((((((((-5328.1988525390625) * tau5) - (1496065.721609933035714286 * tau7)) -
                       (153266.177490234375 * tau6)) -
                      (7447882.350390625 * tau8)) -
                     (22100856.07764756944444444 * tau9)) -
                    (42015571.25390625 * tau10)) -
                   (52708852.66637731481481481 * tau11)) -
                  (22880216.54128086419753086 * tau13)) -
                 (43565653.69020061728395062 * tau12)) -
                (6933398.951903292181069959 * tau14)) -
               (924453.1935871056241426612 * tau15));
    phi(6)  = (((((((((((((107230.0019073486328125 * tau6) + (3911743.44635009765625 * tau7)) +
                         (49592558.34827532087053571 * tau8)) +
                        (327327256.2444545200892857 * tau9)) +
                       (3526512665.739431423611111 * tau11)) +
                      (1320314667.667214471726190 * tau10)) +
                     (6506075517.774877025462963 * tau12)) +
                    (7726873009.8828125 * tau14)) +
                   (8449099135.763165509259259 * tau13)) +
                  (4881459532.087512860082305 * tau15)) +
                 (501515857.5210048010973937 * tau17)) +
                (2031139222.960069444444444 * tau16)) +
               (55723984.16900053345526597 * tau18));
    phi(7)  = ((((((((((((((((-2585008.974552154541015625) * tau7) -
                            (115109601.2886428833007812 * tau8)) -
                           (1814958496.40869140625 * tau9)) -
                          (78228459666.49902913411458 * tau11)) -
                         (15118535470.70343279157366 * tau10)) -
                        (273482436440.6177463107639 * tau12)) -
                       (1223138659577.850206163194 * tau14)) -
                      (677676841014.5960557725694 * tau13)) -
                     (1629074225394.340964988426 * tau15)) -
                    (1152270884857.328880529835 * tau17)) -
                   (1602934313802.891348379630 * tau16)) -
                  (588917533590.4779128086420 * tau18)) -
                 (42266641992.18690462581923 * tau20)) -
                (202879881562.4971422039323 * tau19)) -
               (4025394475.446371869125641 * tau21));
    phi(8)  = (((((((((((((((((338971759786.5465644792884 * tau24) +
                              (72622595.87882459163665771 * tau8)) +
                             (3839296326.682519912719727 * tau9)) +
                            (741010050728.2100511823382 * tau11)) +
                           (72972869796.83574485778809 * tau10)) +
                          (4728818470068.397703588577 * tau12)) +
                         (65097477007191.15089355469 * tau14)) +
                        (20667683341938.66336140951 * tau13)) +
                       (152479523376491.2169460720 * tau15)) +
                      (367482071900019.2868381076 * tau17)) +
                     (270645546827389.1750793457 * tau16)) +
                    (382220257297374.9787138846 * tau18)) +
                   (179020321106247.5114691394 * tau20)) +
                  (302446999586890.6026957948 * tau19)) +
                 (76811168492401.26177650130 * tau21)) +
                (4067661117438.558773751461 * tau23)) +
               (22575519201784.00119432061 * tau22));
    phi(9)  = ((((((((((((((((((((-10820487030026246.15802560) * tau24) -
                                (2329974951.112288981676102 * tau9)) -
                               (3204421618368.187572338364 * tau11)) -
                              (143186772419.2560920119286 * tau10)) -
                             (38658804680937.31829602378 * tau12)) -
                            (1564351256647306.276325171 * tau14)) -
                           (295696972623483.3523212160 * tau13)) -
                          (6033601440981970.440630913 * tau15)) -
                         (39340249994654496.14279334 * tau17)) -
                        (17548357451443068.50438477 * tau16)) -
                       (68913028780415733.17775943 * tau18)) -
                      (103213237253410851.2086602 * tau20)) -
                     (95001751263219690.44791428 * tau19)) -
                    (87982837747558695.44900832 * tau21)) -
                   (29268899384432143.34313245 * tau23)) -
                  (58193404452593856.87994868 * tau22)) -
                 (440154830082830.7139763560 * tau26)) -
                (2772975429521833.498051043 * tau25)) -
               (32604061487617.08992417452 * tau27));
    phi(10) = (((((((((((((((((((((20356704334716668369.03719 * tau24) +
                                  (3526672650910581.893464877 * tau30)) +
                                 (5906050907256.692442968488 * tau11)) +
                                (84053846361.37582501396537 * tau10)) +
                               (152794280429199.7814958678 * tau12)) +
                              (19259274438496755.70484603 * tau14)) +
                             (2146857055341595.890335684 * tau13)) +
                            (120401447638376916.9311805 * tau15)) +
                           (1937739604392837793.628372 * tau17)) +
                          (553463759919030815.3303338 * tau16)) +
                         (5291532561100032069.755284 * tau18)) +
                        (19849372125387206085.07073 * tau20)) +
                       (11452797885661395365.50547 * tau19)) +
                      (27701776760098776464.83803 * tau21)) +
                     (28202657288836392386.51597 * tau23)) +
                    (31172857703295082254.68334 * tau22)) +
                   (5055359283148227235.199949 * tau26)) +
                  (11563080908824599399.17167 * tau25)) +
                 (1641667477501437855.453314 * tau27)) +
                (52900089763658728.40197315 * tau29)) +
               (372945632833794035.2339107 * tau28));
    phi(11) = ((((((((((((((((((((((((-10421163556087996323925.44) * tau24) -
                                    (267313411091476194892.6407 * tau30)) -
                                   (3367884798525.126806809567 * tau11)) -
                                  (266994698245891.9833631720 * tau12)) -
                                 (126736372928041779.7376528 * tau14)) -
                                (7868339035378478.878049733 * tau13)) -
                               (1310856141939127651.253943 * tau15)) -
                              (51013867921762764148.51193 * tau17)) -
                             (9505095396241888043.965687 * tau16)) -
                            (210101665156711621956.3567 * tau18)) -
                           (1766396419169469155581.958 * tau20)) -
                          (680809067349466450984.1967 * tau19)) -
                         (3714461044542513155277.256 * tau21)) -
                        (8996413456126446920199.301 * tau23)) -
                       (6381319935521667021236.195 * tau22)) -
                      (7669111682905600049152.052 * tau26)) -
                     (9899229839297131395937.954 * tau25)) -
                    (4798668661333286324214.086 * tau27)) -
                   (923291756397106785011.7944 * tau29)) -
                  (2388134461189160429354.980 * tau28)) -
                 (6991628530430228603.794118 * tau32)) -
                (54534702537355783109.59412 * tau31)) -
               (423735062450316885.0784314 * tau33));
    phi(12) = (((((((((((((((((((((((((2265181499981701606442225.0 * tau24) +
                                      (1161521061286761581713892.0 * tau30)) +
                                     (148397423935013.3999250465 * tau12)) +
                                    (435451384030404205.1201922 * tau14)) +
                                   (13128818258663568.25439738 * tau13)) +
                                  (7937408485723114496.434515 * tau15)) +
                                 (773568243048000449287.6993 * tau17)) +
                                (93358506268269984746.98223 * tau16)) +
                               (4769376606728078799223.182 * tau18)) +
                              (85555881143101015387051.24 * tau20)) +
                             (22697751563763377591312.77 * tau19)) +
                            (260197218851336240835649.8 * tau21)) +
                           (1328296900615223882220191.0 * tau23)) +
                          (647062456346154782161169.5 * tau22)) +
                         (3822476128794570206783031.0 * tau26)) +
                        (3220716113222632840960304.0 * tau25)) +
                       (3781841350960097251345735.0 * tau27)) +
                      (2103782484922223264132877.0 * tau29)) +
                     (3106620169009037645045708.0 * tau28)) +
                    (178206675733161601307968.3 * tau32)) +
                   (514292924731074141830449.8 * tau31)) +
                  (46546067778486203906898.96 * tau33)) +
                 (1007853846038078711159.049 * tau35)) +
                (8617150383625572980409.870 * tau34)) +
               (55991880335448817286.61384 * tau36));
    chi(0)  = 1.0;
    chi(1)  = (((2.333333333333333333333333 * tau3) + (3.5 * tau2)) + (1.25 * tau));
    chi(2)  = ((((((-25.27777777777777777777778) * tau6) - (75.83333333333333333333333 * tau5)) -
                 (79.95833333333333333333333 * tau4)) -
                (33.625 * tau3)) -
               (4.21875 * tau2));
    chi(3)  = (((((((587.0061728395061728395062 * tau9) + (2641.527777777777777777778 * tau8)) +
                   (4759.027777777777777777778 * tau7)) +
                  (4330.520833333333333333333 * tau6)) +
                 (2047.1125 * tau5)) +
                (453.109375 * tau4)) +
               (31.9921875f * tau3));
    chi(4) =
        ((((((((((-7859.49609375) * tau5) - (373.90869140625 * tau4)) - (54833.22265625 * tau6)) -
              (189793.2208333333333333333 * tau7)) -
             (375736.0005208333333333333 * tau8)) -
            (447530.4166666666666666667 * tau9)) -
           (124738.8117283950617283951 * tau11)) -
          (318189.0131172839506172840 * tau10)) -
         (20789.80195473251028806584 * tau12));
    chi(5) = (((((((((((5889.0618896484375 * tau5) + (167081.378173828125 * tau6)) +
                      (1618334.674874441964285714 * tau7)) +
                     (8018883.473046875 * tau8)) +
                    (23725519.15985243055555556 * tau9)) +
                   (56413247.88802083333333333 * tau11)) +
                  (45020597.45963541666666667 * tau10)) +
                 (46595060.47183641975308642 * tau12)) +
                (7411564.396862139917695473 * tau14)) +
               (24461987.83320473251028807 * tau13)) +
              (988208.5862482853223593964 * tau15));
    chi(6) = ((((((((((((((-116554.3498992919921875) * tau6) - (4209790.52947998046875 * tau7)) -
                        (53046844.35861642020089286 * tau8)) -
                       (348738575.260791015625 * tau9)) -
                      (3740348588.554578993055556 * tau11)) -
                     (1402928860.327645399305556 * tau10)) -
                    (6891971112.248734085648148 * tau12)) -
                   (8173643819.459394290123457 * tau14)) -
                  (8942654096.485291280864198 * tau13)) -
                 (5161879657.283629115226337 * tau15)) -
                (530173906.5222050754458162 * tau17)) -
               (2147389212.053647976680384 * tau16)) -
              (58908211.83580056393842402 * tau18));
    chi(7) =
        (((((((((((((((2776491.120815277099609375 * tau7) + (122708978.7942123413085938 * tau8)) +
                     (1925312041.345973423549107 * tau9)) +
                    (82508725486.56345738002232 * tau11)) +
                   (15983788916.83610055106027 * tau10)) +
                  (287952809783.2672276475694 * tau12)) +
                 (1285107966580.498679832176 * tau14)) +
                (712650296975.9644274450231 * tau13)) +
               (1710504064563.104781539352 * tau15)) +
              (1208933528167.414560828189 * tau17)) +
             (1682287248047.253970550412 * tau16)) +
            (617751113141.8105388374486 * tau18)) +
           (44328429406.43992436366408 * tau20)) +
          (212787605947.7454370522786 * tau19)) +
         (4221755181.565707082253721 * tau21));
    chi(8)  = ((((((((((((((((((-353396089990.2293970103219) * tau24) -
                              (77307924.64520037174224854 * tau8)) -
                             (4062968417.365064620971680 * tau9)) -
                            (778781085218.4028647286551 * tau11)) -
                           (76913846048.55869483947754 * tau10)) -
                          (4959356633780.345308757964 * tau12)) -
                         (68082808844916.25546820747 * tau14)) -
                        (21641132068083.78743394640 * tau13)) -
                       (159329371436569.0371229384 * tau15)) -
                      (383535732852395.0897533275 * tau17)) -
                     (282612499203202.3440902145 * tau16)) -
                    (398767118978192.2197332605 * tau18)) -
                   (186682859410003.5980701839 * tau20)) -
                  (315453595729638.2593476723 * tau19)) -
                 (80088585766097.86195177945 * tau21)) -
                (4240753079882.752764123863 * tau23)) -
               (23536984672244.36711526127 * tau22));
    chi(9)  = (((((((((((((((((((11229654792749015.22701717 * tau24) +
                                (2463116376.890134066343307 * tau9)) +
                               (3360186490494.481565139510 * tau11)) +
                              (150652583603.0128768086433 * tau10)) +
                             (40437382791704.19147014618 * tau12)) +
                            (1630902159338730.110003549 * tau14)) +
                           (308720787738688.0545137042 * tau13)) +
                          (6283277020019534.490934813 * tau15)) +
                         (40904941176357109.48990767 * tau17)) +
                        (18258673913072200.31158583 * tau16)) +
                       (71616091309247803.17800093 * tau18)) +
                      (107183846394401245.5057093 * tau20)) +
                     (98687634471617789.29578043 * tau19)) +
                    (91345622223600985.09434355 * tau21)) +
                   (30378295319311232.27795892 * tau23)) +
                  (60406850807305980.38153996 * tau22)) +
                 (456764446312371.4956358411 * tau26)) +
                (2877683806119897.731818695 * tau25)) +
               (33834403430546.03671376601 * tau27));
    chi(10) = ((((((((((((((((((((((-21053597837396406673.0365) * tau24) -
                                  (3646220876365177.889853517 * tau30)) -
                                 (6184760673789.150449261069 * tau11)) -
                                (88364300020.93355963006616 * tau10)) -
                               (159549747541818.1452511725 * tau12)) -
                              (20033854775307941.92047071 * tau14)) -
                             (2236938226816488.952569230 * tau13)) -
                            (125079628228935167.7962316 * tau15)) -
                           (2009277564547413368.390572 * tau17)) -
                          (574371677952986224.0321184 * tau16)) -
                         (5483258169122136680.737239 * tau18)) -
                        (20548963029842481258.16041 * tau20)) -
                       (11861476998383736648.72274 * tau19)) -
                      (28668438826053413980.12674 * tau21)) -
                     (29172870122949767579.89601 * tau23)) -
                    (32252013255464105175.54592 * tau22)) -
                   (5227324779598662081.545248 * tau26)) -
                  (11957474746245008204.2076 * tau25)) -
                 (1697408519086965241.663463 * tau27)) -
                (54693313145477668.34780275 * tau29)) -
               (385594378487915085.2699942 * tau28));
    chi(11) = (((((((((((((((((((((((10749119360152201659095.30 * tau24) +
                                    (275549372270078767749.3455 * tau30)) +
                                   (3524530603107.690844335593 * tau11)) +
                                  (278502821903198.7234563101 * tau12)) +
                                 (131627099146616557.0788647 * tau14)) +
                                (8187462234898173.302466749 * tau13)) +
                               (1359407981935523997.779927 * tau15)) +
                              (52789252821617587071.26402 * tau17)) +
                             (9845374898471241989.909180 * tau16)) +
                            (217243684855444804880.3725 * tau18)) +
                           (1824328572800792450321.910 * tau20)) +
                          (703503512901659205452.3662 * tau19)) +
                         (3834651403368051540134.680 * tau21)) +
                        (9281659895577361070060.539 * tau23)) +
                       (6585513896734531838538.180 * tau22)) +
                      (7907844921517711462070.132 * tau26)) +
                     (10208879014255827826055.24 * tau25)) +
                    (4947481686715095751572.891 * tau27)) +
                   (951780479796384621502.0100 * tau29)) +
                  (2461979171794923193819.241 * tau28)) +
                 (7206755562135774099.295476 * tau32)) +
                (56213398719189220090.88340 * tau31)) +
               (436773064371865096.9269985 * tau33));
    chi(12) = ((((((((((((((((((((((((((-2331551239435886577025091.0) * tau24) -
                                      (1194449360583363232233128.0 * tau30)) -
                                     (154712207932248.0126878145 * tau12)) -
                                    (451767349373865480.3252141 * tau14)) -
                                   (13649617202765419.19403282 * tau13)) -
                                  (8221062051469232071.896864 * tau15)) -
                                 (799266188640123245085.8777 * tau17)) -
                                (96565089437574073764.48836 * tau16)) -
                               (4923418321891009149545.962 * tau18)) -
                              (88200546293836821900475.57 * tau20)) -
                             (23413619423670496869405.09 * tau19)) -
                            (268104858375645682223269.2 * tau21)) -
                           (1367612737884965288180114.0 * tau23)) -
                          (666447368381836165309217.2 * tau22)) -
                         (3932747815369540922524828.0 * tau26)) -
                        (3314285837643629657858286.0 * tau25)) -
                       (3890307474197961385216392.0 * tau27)) -
                      (2163604453746493613327353.0 * tau29)) -
                     (3195297616564396799146471.0 * tau28)) -
                    (183238172201351976450037.1 * tau32)) -
                   (528838703410831864182048.1 * tau31)) -
                  (47858665362437589747700.5 * tau33)) -
                 (1036244095222249942459.304 * tau35)) -
                (8859971761162727071404.065 * tau34)) -
               (57569116401236107914.40578 * tau36));
    mu2     = (2.0 * a);
    maxi    = 1.0;
    k       = 0;
    while (((k < 13) && (maxi > eps))) {
      if ((k == 0)) {
        musp  = 1.0;
        musn  = 1.0;
        sphip = 0.0;
        schip = 0.0;
        sphin = 0.0;
        schin = 0.0;
      } else {
        musn = (musn * mu2);
        musp = ((-musp) * mu2);
      }
      delta1 = (phi(k) / musp);
      delta2 = (chi(k) / musp);
      delta3 = (phi(k) / musn);
      delta4 = (chi(k) / musn);
      sphip  = (sphip + delta1);
      schip  = (schip + delta2);
      sphin  = (sphin + delta3);
      schin  = (schin + delta4);
      e1     = Kokkos::abs((delta1 / sphip));
      e2     = Kokkos::abs((delta2 / schip));
      e3     = Kokkos::abs((delta3 / sphin));
      e4     = Kokkos::abs((delta4 / schin));
      k      = ((k + 1));
      maxi = ((((((e1) > (e2) ? (e1) : (e2))) > (e3) ? (((e1) > (e2) ? (e1) : (e2))) : (e3))) > (e4)
                  ? (((((e1) > (e2) ? (e1) : (e2))) > (e3) ? (((e1) > (e2) ? (e1) : (e2))) : (e3)))
                  : (e4));
    }
    xargu   = (((x * x) * 0.25) + a);
    facu    = Kokkos::pow(xargu, 0.25);
    uaxx(0) = ((1.0 / (facu * sqrt2)) * sphip);
    uaxx(1) = (((-facu) / sqrt2) * schip);
    sinpia  = Kokkos::sin((pi * a));
    mu2     = (mu * mu);
    mu4     = (mu2 * mu2);
    for (k = 0; k < 26; k += 1) {
      if ((k == 0)) {
        mus  = 1.0;
        gmus = 0.0;
      } else {
        mus = (mus * mu2);
      }
      gmus = (gmus + (gk(k) / mus));
    }
    for (k = 0; k < 13; k += 1) {
      if ((k == 0)) {
        mus  = 1.0;
        smus = 0.0;
      } else {
        mus = (mus * mu4);
      }
      smus = (smus + (wk(k) / mus));
    }
    dl = (((-2.0) * a) *
          ((t * Kokkos::sqrt(((t * t) + 1.0))) + Kokkos::log((t + Kokkos::sqrt(((t * t) + 1.0))))));
    if ((dl < Kokkos::log(dwarf))) {
      dnew = 0.0;
    } else {
      dnew = (((Kokkos::exp(dl) * sqrttwopi) * smus) / (gmus * gmus));
    }
    vaxx(0) =
        ((1.0 / facu) * ((((sinpia * dnew) * sphip) / (sqrt2 * pi)) + (sphin / Kokkos::sqrt(pi))));
    vaxx(1) = ((-facu) * ((((sinpia * dnew) * schip) / (sqrt2 * pi)) - (schin / Kokkos::sqrt(pi))));
    ffa     = Kokkos::exp((
        ((a * Kokkos::log(((x * 0.5) + Kokkos::sqrt(xargu)))) + ((x * 0.5) * Kokkos::sqrt(xargu))) -
        (a * 0.5)));
    uaxx(0) = (uaxx(0) / ffa);
    uaxx(1) = (uaxx(1) / ffa);
    vaxx(0) = (vaxx(0) * ffa);
    vaxx(1) = (vaxx(1) * ffa);
  }
}
inline void expax(double a, double x, int mode, ParabolicArray<2>& uax, ParabolicArray<2>& vax,
                  int& ierr, double epss) {
  // Evaluate the Native expax algorithm in call-local storage.
  if (a < 1.5 || mode != 0)
    throw std::runtime_error(
        "Algorithm850 helper request is outside the supported KK call closure");
  double eps;
  double aph;
  double amh;
  double aphl;
  double amhl;
  double err1;
  double err2;
  double err1p;
  double err2p;
  double x2;
  double a2;
  double x2k;
  double b2;
  double y1;
  double y2;
  double acof;
  double bcof;
  double facto1;
  double facto2;
  double y1p;
  double y2p;
  double a2p;
  double b2p;
  double acofd;
  double bcofd;
  double lambda;
  double facto1p;
  double facto2p;
  double sqrtx;
  double facu;
  double afacu;
  double phiax;
  double over;
  double xargu;
  double ffa;
  int k;
  int l;
  int m;
  eps   = epss;
  over  = (giant * 1.e-5);
  m     = 0;
  ierr  = 0;
  aph   = (a + 0.5);
  amh   = (a - 0.5);
  x2    = (x * x);
  sqrtx = Kokkos::sqrt(x);
  xargu = (((x * x) * 0.25) + a);
  ffa =
      (((a * Kokkos::log(((x * 0.5) + Kokkos::sqrt(xargu)))) + ((x * 0.5) * Kokkos::sqrt(xargu))) -
       (a * 0.5));
  if (((ffa > Kokkos::log(giant)) || (ffa < Kokkos::log(dwarf)))) {
    uax(0) = dwarf;
    uax(1) = dwarf;
    vax(0) = giant;
    vax(1) = giant;
    ierr   = 1;
  }
  if ((ierr == 0)) {
    phiax   = Kokkos::exp((((0.25 * x) * x) + (a * Kokkos::log(x))));
    facto1  = (1.0 / (sqrtx * phiax));
    facto1p = (((-sqrtx) * 0.5) / phiax);
    facto2  = ((sqrt2 / (sqrtpi * sqrtx)) * phiax);
    facto2p = ((sqrtx / (sqrtpi * sqrt2)) * phiax);
    err1    = 1.0;
    err2    = 1.0;
    err1p   = 1.0;
    err2p   = 1.0;
    y1      = 1.0;
    y2      = 1.0;
    y1p     = 1.0;
    y2p     = 1.0;
    k       = 1;
    a2      = 1.0;
    b2      = 1.0;
    x2k     = 1.0;
    while ((((((err1 > eps) || (err2 > eps)) || (err1p > eps)) || (err2p > eps)) && (m < 500))) {
      l    = ((2 * k));
      aphl = (aph + l);
      amhl = (amh - l);
      a2p  = (((-amhl) * ((amh + l) - 1.0)) * a2);
      a2   = (((-(aphl - 2.0)) * (aphl - 1.0)) * a2);
      b2   = (((amhl + 2.0) * (amhl + 1.0)) * b2);
      if ((Kokkos::abs((aph - l)) < dwarf)) {
        b2p = 0.0;
      } else {
        b2p = ((aphl / (aph - l)) * b2);
      }
      x2k   = ((l * x2) * x2k);
      acof  = (a2 / x2k);
      bcof  = (b2 / x2k);
      acofd = (a2p / x2k);
      bcofd = (b2p / x2k);
      y1    = (y1 + acof);
      y2    = (y2 + bcof);
      y1p   = (y1p + acofd);
      y2p   = (y2p + bcofd);
      k     = ((k + 1));
      err1  = Kokkos::abs((acof / y1));
      err2  = Kokkos::abs((bcof / y2));
      err1p = Kokkos::abs((acofd / y1p));
      err2p = Kokkos::abs((bcofd / y2p));
      if ((Kokkos::abs(a2) > over)) {
        m = 500;
      }
      if ((Kokkos::abs(b2) > over)) {
        m = 501;
      }
    }
    if ((m == 500)) {
      ierr   = 1;
      uax(0) = giant;
      uax(1) = giant;
    } else {
      uax(0) = (facto1 * y1);
      uax(1) = (facto1p * y1p);
    }
    if ((m == 501)) {
      ierr   = 1;
      vax(0) = giant;
      vax(1) = giant;
    } else {
      vax(0) = (facto2 * y2);
      vax(1) = (facto2p * y2p);
    }
  }
}
} // namespace clubb_kokkos::parabolic
