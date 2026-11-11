// Solve real-coefficient cubic and quadratic equations, retaining complex roots.
//
// cubic_solve(nz,...) and quadratic_solve(nz,...) return owned (level,root)
// LayoutLeft Views in the coefficient Views' memory space. Their leading
// coefficients must be nonzero and inputs conformable, as in Native.
// All coefficient Views must be accessible to the first View's execution space.
// cube_root returns the real signed cube root through the Native real power.
// Scalar cubic/quadratic overloads evaluate one polynomial inside a device
// kernel; they retain the same equations and root order as the profile APIs.
// No storage persists between calls.
// Reference: src/CLUBB_core/calc_roots.F90 and Cardano's cubic formula.
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"
#include <Kokkos_Complex.hpp>

namespace clubb_kokkos {

using complex_t = Kokkos::complex<real_t>;

// Device scalar overloads express Native array equations one level at a time;
// forward declarations preserve the Native profile procedure order below.
KOKKOS_INLINE_FUNCTION Kokkos::Array<complex_t, 3> cubic_solve(real_t a_coef, real_t b_coef,
                                                               real_t c_coef, real_t d_coef);
KOKKOS_INLINE_FUNCTION Kokkos::Array<complex_t, 2> quadratic_solve(real_t a_coef, real_t b_coef,
                                                                   real_t c_coef);

// =============================================================================
template <class A, class B, class C, class D>
Kokkos::View<complex_t **, Kokkos::LayoutLeft, typename A::memory_space>
cubic_solve(int nz, A a_coef, B b_coef, C c_coef, D d_coef) {
  // Solve a*x^3+b*x^2+c*x+d=0 at every vertical level with Cardano's formula,
  // returning all three complex roots in Native order. a must be nonzero.
  // Reference: http://mathworld.wolfram.com/CubicFormula.html
  static_assert(A::rank == 1 && B::rank == 1 && C::rank == 1 && D::rank == 1);
  using Execution = typename A::execution_space;
  Kokkos::View<complex_t **, Kokkos::LayoutLeft, typename A::memory_space> roots(
      "cubic_solve::roots", nz, 3);
  Kokkos::parallel_for(
      "cubic_solve", Kokkos::RangePolicy<Execution>(0, nz), KOKKOS_LAMBDA(const int k) {
        const auto level_roots = cubic_solve(a_coef(k), b_coef(k), c_coef(k), d_coef(k));
        roots(k, 0) = level_roots[0];
        roots(k, 1) = level_roots[1];
        roots(k, 2) = level_roots[2];
      });
  return roots;
}

// =============================================================================
template <class A, class B, class C>
Kokkos::View<complex_t **, Kokkos::LayoutLeft, typename A::memory_space>
quadratic_solve(int nz, A a_coef, B b_coef, C c_coef) {
  // Solve a*x^2+b*x+c=0 at every vertical level, returning the plus-square-root
  // solution first and minus-square-root solution second. a must be nonzero.
  static_assert(A::rank == 1 && B::rank == 1 && C::rank == 1);
  using Execution = typename A::execution_space;
  Kokkos::View<complex_t **, Kokkos::LayoutLeft, typename A::memory_space> roots(
      "quadratic_solve::roots", nz, 2);
  Kokkos::parallel_for(
      "quadratic_solve", Kokkos::RangePolicy<Execution>(0, nz), KOKKOS_LAMBDA(const int k) {
        const auto level_roots = quadratic_solve(a_coef(k), b_coef(k), c_coef(k));
        roots(k, 0) = level_roots[0];
        roots(k, 1) = level_roots[1];
      });
  return roots;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION real_t cube_root(real_t x) {
  // Calculate the real cube root. A negative base is evaluated as -abs(x)^(1/3)
  // to avoid NaNs caused by the finite representation of the exponent.
  real_t result;
  if (x >= zero) {
    result = Kokkos::pow(x, one_third);
  } else {
    result = -Kokkos::pow(Kokkos::abs(x), one_third);
  }
  return result;
}

// Scalar execution adapters for the two Native profile routines. These are
// overloads of the same APIs, not alternate polynomial/root algorithms.
// =============================================================================
KOKKOS_INLINE_FUNCTION Kokkos::Array<complex_t, 3> cubic_solve(real_t a_coef, real_t b_coef,
                                                               real_t c_coef, real_t d_coef) {
  // Evaluate the Native Cardano equations for one level. For D<0 all roots are
  // real; D=0 gives repeated real roots; D>0 gives one real root and a conjugate
  // pair. Retain every complex conversion and the Native S/T sign correction.
  // Reference: http://mathworld.wolfram.com/CubicFormula.html
  real_t cap_Q_coef, cap_R_coef, discriminant, sgn_S, sgn_T;
  complex_t sqrt_det, cap_S_coef, cap_T_coef;
  // Value materialization avoids taking host constexpr addresses in CUDA complex constructors.
  const complex_t i_cmplx{real_t(zero), real_t(one)};
  complex_t sqrt_3, one_half_cmplx, one_third_cmplx;
  Kokkos::Array<complex_t, 3> roots;

  // Declare constants as complex values to retain the Native conversions.
  sqrt_3 = complex_t(Kokkos::sqrt(three), real_t(zero));
  one_half_cmplx = complex_t(real_t(one_half), real_t(zero));
  one_third_cmplx = complex_t(real_t(one_third), real_t(zero));

  // Q = (3*(c/a)-(b/a)^2)/9.
  cap_Q_coef = (three * (c_coef / a_coef) - (b_coef / a_coef) * (b_coef / a_coef)) / real_t(9);

  // R = (9*(b/a)*(c/a)-27*(d/a)-2*(b/a)^3)/54.
  cap_R_coef = (real_t(9) * (b_coef / a_coef) * (c_coef / a_coef) - real_t(27) * (d_coef / a_coef) -
                two * ((b_coef / a_coef) * (b_coef / a_coef) * (b_coef / a_coef))) /
               real_t(54);

  // D = Q^3+R^2. Its square root can be complex.
  discriminant = (cap_Q_coef * cap_Q_coef * cap_Q_coef) + (cap_R_coef * cap_R_coef);
  sqrt_det = Kokkos::sqrt(complex_t(discriminant, real_t(zero)));

  // S=(R+sqrt(D))^(1/3), T=(R-sqrt(D))^(1/3).
  // Take the cube root of -w and negate it when Re(w)<0 so a negative real
  // base has the real cube root, preserving the S*T=-Q Cardano pairing.
  cap_S_coef = complex_t(cap_R_coef, real_t(zero)) + sqrt_det;
  cap_T_coef = complex_t(cap_R_coef, real_t(zero)) - sqrt_det;
  sgn_S = Kokkos::copysign(one, cap_S_coef.real());
  sgn_T = Kokkos::copysign(one, cap_T_coef.real());
  cap_S_coef = sgn_S * Kokkos::pow(sgn_S * cap_S_coef, one_third_cmplx);
  cap_T_coef = sgn_T * Kokkos::pow(sgn_T * cap_T_coef, one_third_cmplx);

  // The first root is always real-valued: -(1/3)*(b/a)+(S+T).
  roots[0] =
      -one_third_cmplx * complex_t(b_coef / a_coef, real_t(zero)) + (cap_S_coef + cap_T_coef);
  // Root2 is real for D<=0, otherwise the conjugate of root3.
  roots[1] = -one_third_cmplx * complex_t(b_coef / a_coef, real_t(zero)) -
             one_half_cmplx * (cap_S_coef + cap_T_coef) +
             one_half_cmplx * i_cmplx * sqrt_3 * (cap_S_coef - cap_T_coef);
  // Root3 retains the opposite sign on the imaginary contribution.
  roots[2] = -one_third_cmplx * complex_t(b_coef / a_coef, real_t(zero)) -
             one_half_cmplx * (cap_S_coef + cap_T_coef) -
             one_half_cmplx * i_cmplx * sqrt_3 * (cap_S_coef - cap_T_coef);
  return roots;
}

// =============================================================================
KOKKOS_INLINE_FUNCTION Kokkos::Array<complex_t, 2> quadratic_solve(real_t a_coef, real_t b_coef,
                                                                   real_t c_coef) {
  // Evaluate the Native quadratic formula for one level. D>0 gives two real
  // roots, D=0 a double root, and D<0 a conjugate pair. Keep the source order.
  real_t discriminant;
  complex_t sqrt_det;
  Kokkos::Array<complex_t, 2> roots;
  discriminant = b_coef * b_coef - four * a_coef * c_coef;
  sqrt_det = Kokkos::sqrt(complex_t(discriminant, real_t(zero)));
  roots[0] = (-complex_t(b_coef, real_t(zero)) + sqrt_det) / complex_t(two * a_coef, real_t(zero));
  roots[1] = (-complex_t(b_coef, real_t(zero)) - sqrt_det) / complex_t(two * a_coef, real_t(zero));
  return roots;
}

} // namespace clubb_kokkos
