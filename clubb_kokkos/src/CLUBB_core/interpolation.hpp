// Interpolate scalar values and vertical profiles, preserving native searches.
//
// lin_interpolate_two_points and linear_interp_factor evaluate linear formulas;
// mono_cubic_interp provides Steffen's monotone cubic and the native quintic
// option. These scalar formulas are device callable. binary_search returns a
// zero-based upper bracket, or -1 outside the source range.
// plinterp_fnc/zlinterp_fnc interpolate pressure/height profiles with zero fill;
// lin_interp_between_grids clamps outside endpoints and honors the native
// 1e-6 altitude tolerance. Profile results use their input values' memory space.
// lin_interpolate_on_grid_api is the host-facing scalar/diagnostic interface.
// Coordinates are heights [m] or pressures [Pa]; interpolated units follow input.
// References: WRF-HOC LIN_INT; Steffen (1990), A&A 239, 443–450.
#pragma once

#include "error_code.hpp"
#include "array_types.hpp"
#include "constants_clubb.hpp"
#include "model_flags.hpp"

#include <cstdint>
#include <type_traits>
#include <iostream>
#include <stdexcept>

namespace clubb_kokkos {

KOKKOS_INLINE_FUNCTION real_t lin_interpolate_two_points(real_t height_int, real_t height_high,
                                                         real_t height_low, real_t var_high,
                                                         real_t var_low) {
  // Compute a linear interpolation between two known heights and values.
  //
  // ################################ Height high, known variable value
  // -------------------------------- Height to be interpolated to
  // ################################ Height low, known variable value
  //
  // variable(@ Height interpolation) =
  // [(variable(@ Height high)-variable(@ Height low))/(Height high-Height low)]
  // * (Height interpolation-Height low) + variable(@ Height low).
  // Comments from WRF-HOC, Brian Griffin.

  // Check for valid input. Native CLUBB_GPU omits this host STOP check.
  KOKKOS_IF_ON_HOST(
      (if (Kokkos::abs(height_low - height_high) < real_t(1.0e-12)) throw std::runtime_error(
           "lin_interpolate_two_points: height_high and height_low cannot be equal.");))

  // Compute linear interpolation.
  real_t lin_interpolate_two_points =
      ((height_int - height_low) / (height_high - height_low)) * (var_high - var_low) + var_low;
  return lin_interpolate_two_points;
}

KOKKOS_INLINE_FUNCTION real_t linear_interp_factor(real_t factor, real_t var_high, real_t var_low) {
  // Determine the value of a linear interpolation using its dimensionless factor.
  real_t linear_interp_factor = factor * (var_high - var_low) + var_low;
  return linear_interp_factor;
}

KOKKOS_INLINE_FUNCTION real_t mono_cubic_interp(real_t z_in, int km1, int k00, int kp1, int kp2,
                                                real_t zm1, real_t z00, real_t zp1, real_t zp2,
                                                real_t fm1, real_t f00, real_t fp1, real_t fp2) {
  // Steffen's monotone cubic interpolation method; return the value between
  // z00 and zp1. Levels choose endpoint/interior/extrapolation branches; their
  // relative order is unchanged by zero-based indexing. Field units are retained.
  // Original author: Takanobu Yamaguchi, tak.yamaguchi@noaa.gov.
  // Modified for CLUBB, adding the 3/2 from eq.21 (dschanen,26 Oct2011) and
  // a quintic option. Reference: Steffen (1990), A&A 239,443–450.
  real_t f_out;

  // Local Variables
  real_t hm1, h00, hp1;            // Neighboring height intervals [m]
  real_t sm1, s00, sp1;            // Neighboring slopes [field units/m]
  real_t p00, pp1, dfdx00, dfdxp1; // Estimated/limited derivatives [field units/m]
  real_t c1, c2, c3, c4;           // Cubic coefficients [field units/m^3, /m^2, /m, field units]
  real_t w00, wp1;                 // Extrapolation weights [-]
  real_t coef1, coef2;             // Derivative-limiter coefficients [-]
  real_t zprime, beta, alpha, zn;  // Shifted height [m] and quintic coefficients

  // ---- Begin Code ----
  coef1 = real_t(1.0);
  coef2 = real_t(1.0);

  // Explicit fma follows the current NVHPC -O0 contraction: in a sum of two
  // products, round the first product and fuse the second with the sum. CUDA
  // can choose the opposite contraction for the same source expression.

  if (km1 <= k00) {
    hm1 = z00 - zm1;
    h00 = zp1 - z00;
    hp1 = zp2 - zp1;

    if (km1 == k00) {
      s00    = (fp1 - f00) / (zp1 - z00);
      sp1    = (fp2 - fp1) / (zp2 - zp1);
      dfdx00 = s00;
      pp1    = Kokkos::fma(sp1, h00, s00 * hp1) / (h00 + hp1);
      dfdxp1 = coef1 * (Kokkos::copysign(real_t(1.0), s00) + Kokkos::copysign(real_t(1.0), sp1)) *
               Kokkos::min(Kokkos::min(Kokkos::abs(s00), Kokkos::abs(sp1)),
                           coef2 * real_t(0.5) * Kokkos::abs(pp1));
    } else if (kp1 == kp2) {
      sm1    = (f00 - fm1) / (z00 - zm1);
      s00    = (fp1 - f00) / (zp1 - z00);
      p00    = Kokkos::fma(s00, hm1, sm1 * h00) / (hm1 + h00);
      dfdx00 = coef1 * (Kokkos::copysign(real_t(1.0), sm1) + Kokkos::copysign(real_t(1.0), s00)) *
               Kokkos::min(Kokkos::min(Kokkos::abs(sm1), Kokkos::abs(s00)),
                           coef2 * real_t(0.5) * Kokkos::abs(p00));
      dfdxp1 = s00;
    } else {
      sm1    = (f00 - fm1) / (z00 - zm1);
      s00    = (fp1 - f00) / (zp1 - z00);
      sp1    = (fp2 - fp1) / (zp2 - zp1);
      p00    = Kokkos::fma(s00, hm1, sm1 * h00) / (hm1 + h00);
      pp1    = Kokkos::fma(sp1, h00, s00 * hp1) / (h00 + hp1);
      dfdx00 = coef1 * (Kokkos::copysign(real_t(1.0), sm1) + Kokkos::copysign(real_t(1.0), s00)) *
               Kokkos::min(Kokkos::min(Kokkos::abs(sm1), Kokkos::abs(s00)),
                           coef2 * real_t(0.5) * Kokkos::abs(p00));
      dfdxp1 = coef1 * (Kokkos::copysign(real_t(1.0), s00) + Kokkos::copysign(real_t(1.0), sp1)) *
               Kokkos::min(Kokkos::min(Kokkos::abs(s00), Kokkos::abs(sp1)),
                           coef2 * real_t(0.5) * Kokkos::abs(pp1));
    }

    c1 = Kokkos::fma(real_t(2.), -s00, dfdx00 + dfdxp1) / (h00 * h00);
    c2 = (Kokkos::fma(real_t(3.), s00, -(real_t(2.) * dfdx00)) - dfdxp1) / h00;
    c3 = dfdx00;
    c4 = f00;

    if (!l_quintic_poly_interp) {
      // Old formula:
      // f_out = c1*((z_in-z00)**3) + c2*((z_in-z00)**2) + c3*(z_in-z00) + c4
      // Faster nested multiplication.
      zprime = z_in - z00;
      f_out  = Kokkos::fma(zprime, Kokkos::fma(zprime, Kokkos::fma(zprime, c1, c2), c3), c4);
    } else {
      // Use a quintic polynomial instead of the Steffen formula.
      // Unlike the formula above, this does not guarantee monotonicity.
      beta = real_t(120.) * Kokkos::fma(real_t(0.5) * h00, -(dfdx00 + dfdxp1), fp1 - f00);

      // Prevent an underflow by using a linear interpolation.
      if (Kokkos::abs(beta) < eps) {
        f_out = lin_interpolate_two_points(z00, zp1, zm1, // In
                                           fp1, fm1);     // In
      } else {
        alpha = Kokkos::fma((real_t(6.) / beta) * h00, dfdxp1 - dfdx00, real_t(0.5));
        zn    = (z_in - z00) / h00;
        f_out = Kokkos::fma(
            Kokkos::fma(dfdx00, h00,
                        Kokkos::fma(Kokkos::fma(beta / real_t(20.), zn,
                                                -(beta * (real_t(1.) + alpha) / real_t(12.))),
                                    zn, beta * alpha / real_t(6.)) *
                            (zn * zn)),
            zn, f00);
      }
    }
  } else {
    // Linear extrapolation.
    wp1   = (z_in - z00) / (zp1 - z00);
    w00   = real_t(1.) - wp1;
    f_out = Kokkos::fma(w00, f00, wp1 * fp1);
  }
  return f_out;
}

template <class View> KOKKOS_INLINE_FUNCTION int binary_search(int n, View array, real_t var) {
  // Return the upper bracket of a sorted profile, or -1 outside its bounds.
  // Native low=2 prevents using array(0) when the lower endpoint is selected;
  // C++ starts at1 and returns zero-based indices (including an exact first point).
  int i, high, low;
  low  = 1;
  high = n - 1;

  // Ensure variable is within range and array has more than1 element.
  if (n < 2 || var < array(0) || var > array(n - 1))
    return -1;

  // Special case to cover var == array(1) in native indexing.
  if (var >= array(0) && var <= array(1))
    return 1;

  while (low <= high) {
    // Convert the native midpoint, including its one-based integer division.
    i = (low + high + 2) / 2 - 1;
    if (var > array(i - 1) && var <= array(i))
      return i;
    else if (var < array(i))
      high = i - 1;
    else if (var > array(i))
      low = i + 1;
  }
  KOKKOS_IF_ON_HOST((std::cerr << "Logic error in function binary_search.\n";))
  KOKKOS_IF_ON_DEVICE((Kokkos::printf("Logic error in function binary_search.\n");))
  return -1;
}

// Forward declaration preserves the native plinterp/zlinterp definition order.
template <class Out, class Source, class Values>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>
zlinterp_fnc(int dim_out, int dim_src, Out grid_out, Source grid_src, Values var_src);

template <class Out, class Source, class Values>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>
plinterp_fnc(int dim_out, int dim_src, Out grid_out, Source grid_src, Values var_src) {
  // Interpolate in pressure by negating both grids, then use native height
  // interpolation. Values outside the source range are zero (WRF-HOC LIN_INT).
  using Result    = Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>;
  using Execution = typename Result::execution_space;
  Result grid_out_negative(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "plinterp_fnc::negative_grid_out"), dim_out);
  Result grid_src_negative(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "plinterp_fnc::negative_grid_src"), dim_src);
  Kokkos::parallel_for(
      "plinterp_fnc::negate_out", Kokkos::RangePolicy<Execution>(0, dim_out),
      KOKKOS_LAMBDA(const int k) { grid_out_negative(k) = -grid_out(k); });
  Kokkos::parallel_for(
      "plinterp_fnc::negate_src", Kokkos::RangePolicy<Execution>(0, dim_src),
      KOKKOS_LAMBDA(const int k) { grid_src_negative(k) = -grid_src(k); });
  auto var_out = zlinterp_fnc(dim_out, dim_src, grid_out_negative, // In
                              grid_src_negative, var_src);         // In
  return var_out;
}

template <class Out, class Source, class Values, class Result>
void zlinterp_fnc(int dim_out, int dim_src, Out grid_out, Source grid_src, Values var_src,
                  Result var_out) {
  // Interpolate a height profile with zero fill below/above its range.
  // Preserve the native zero-tail early exit after the first missing upper
  // bracket, even if a caller supplies nonmonotone output coordinates.
  // Reference: WRF-HOC function LIN_INT.
  static_assert(Out::rank == 1 && Source::rank == 1 && Values::rank == 1 && Result::rank == 1);
  using Execution = typename Values::execution_space;
  static_assert(Kokkos::SpaceAccessibility<Execution, typename Result::memory_space>::accessible);

  // ---- Begin Code ----
  Kokkos::parallel_for(
      "zlinterp_fnc", Kokkos::RangePolicy<Execution>(0, 1), KOKKOS_LAMBDA(const int) {
        int k = 0;
        for (int kint = 0; kint < dim_out; ++kint) {
          // Set to0 if below the input data's lowest point.
          if (grid_out(kint) < grid_src(0)) {
            var_out(kint) = real_t(0.0);
            continue;
          }
          // Native sequential search was replaced by binary search (Fasching,2008).
          k = binary_search(dim_src, grid_src, grid_out(kint));

          // Above the data, set this point and all those above it to zero.
          if (k == -1) {
            for (int remaining = kint; remaining < dim_out; ++remaining)
              var_out(remaining) = real_t(0.0);
            break;
          }
          const int km1 = Kokkos::max(0, k - 1);
          // Interpolate. Operand placement changed the historic WRF-HOC results;
          // retain the standard function's ordering (Joshua Fasching,April2008).
          var_out(kint) = lin_interpolate_two_points(grid_out(kint), grid_src(k),              // In
                                                     grid_src(km1), var_src(k), var_src(km1)); // In
        }
      });
}

template <class Out, class Source, class Values>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>
zlinterp_fnc(int dim_out, int dim_src, Out grid_out, Source grid_src, Values var_src) {
  // Retain an independently owned result for callers that keep the interpolated profile.
  using Result = Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>;
  Result var_out(Kokkos::view_alloc(Kokkos::WithoutInitializing, "zlinterp_fnc::var_out"), dim_out);
  zlinterp_fnc(dim_out, dim_src, grid_out, grid_src, var_src, var_out);
  return var_out;
}

template <class X, class T>
void lin_interpolate_on_grid_api(int nparam,       // In
                                 X xlist, T tlist, // InOut
                                 real_t xvalue,    // In
                                 real_t& tvalue) { // InOut
  // Linear interpolation for the25 June1996 altocumulus case (Michael Falk,
  // COAMPS). xlist gives sorted independent coordinates, tlist their values;
  // xvalue selects tvalue. Native checks stop on unsorted/out-of-range input.
  const auto xlist_h = to_host_space(xlist);
  const auto tlist_h = to_host_space(tlist);
  int bottombound, topbound;

  // Assure elements are ordered so interpolation uses the two closest points.
  if (clubb_at_least_debug_level_api(0))
    for (int i = 1; i < nparam; ++i)
      if (xlist_h(i) <= xlist_h(i - 1))
        throw std::runtime_error("xlist must be sorted for lin_interpolate_on_grid_api.");

  // Outside the largest/smallest x-value, crash.
  if (clubb_at_least_debug_level_api(0))
    if (xvalue < xlist_h(0) || xvalue > xlist_h(nparam - 1))
      throw std::runtime_error("lin_interpolate_on_grid_api: Value out of range");

  // Find top/bottom bounds, interpolate and return the value.
  topbound    = -1;
  bottombound = -1;
  for (int i = 1; i < nparam; ++i)
    if (xvalue >= xlist_h(i - 1) && xvalue <= xlist_h(i)) {
      bottombound = i - 1;
      topbound    = i;
    }

  if (clubb_at_least_debug_level_api(1) && (topbound == -1 || bottombound == -1))
    std::cerr
        << "Sanity check failed! xlist is not properly sorted\nin lin_interpolate_on_grid_api.\n";

  tvalue = lin_interpolate_two_points(xvalue, xlist_h(topbound), xlist_h(bottombound), // In
                                      tlist_h(topbound), tlist_h(bottombound));        // In
}

template <class Out, class Source, class Values, class Output>
void lin_interp_between_grids(ScratchWorkspace& workspace,                         // InOut
                              int size_interpolate, int size_current,              // In
                              Out interpolate_altitudes, Source current_altitudes, // In
                              Values current_values,                               // In
                              Output lin_interp_between_grids) {                   // Out
  // Interpolate current_values onto interpolate_altitudes; retain the first
  // source match within tol and use the closest endpoint outside the grid.
  using Result    = Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>;
  using Execution = typename Result::execution_space;
  auto workspace_scope = workspace.scope();
  const auto overlaps  = [&](const auto& input) {
    if constexpr (!std::is_same_v<typename std::remove_cvref_t<decltype(input)>::memory_space,
                                  typename Output::memory_space>) {
      return false;
    } else {
      if (!input.size() || !lin_interp_between_grids.size())
        return false;
      const auto begin  = reinterpret_cast<std::uintptr_t>(input.data());
      const auto output = reinterpret_cast<std::uintptr_t>(lin_interp_between_grids.data());
      return begin <= output ? output - begin < input.span() * sizeof(real_t)
                             : begin - output < lin_interp_between_grids.span() * sizeof(real_t);
    }
  };
  if (overlaps(interpolate_altitudes) || overlaps(current_altitudes) || overlaps(current_values)) {
    auto result = workspace.view<Result>(Kokkos::WithoutInitializing,
                                         "lin_interp_between_grids::overlap", size_interpolate);
    clubb_kokkos::lin_interp_between_grids(workspace, size_interpolate, size_current,
                                           interpolate_altitudes, current_altitudes, current_values,
                                           result);
    Kokkos::deep_copy(lin_interp_between_grids, result);
    return;
  }
  const real_t tol = real_t(1.0e-6);

  // Assure elements are ordered so interpolation uses the two closest points.
  if (clubb_at_least_debug_level_api(0)) {
    const auto host = to_host_space(workspace, current_altitudes);
    for (int i = 1; i < size_current; ++i)
      if (host(i) <= host(i - 1))
        throw std::runtime_error("current_altitudes must be sorted for lin_interp_between_grids.");

    // Outside the largest/smallest x-value, give a warning.
    const auto current = host;
    const auto target  = to_host_space(workspace, interpolate_altitudes);
    bool outside       = false;
    for (int i = 0; i < size_interpolate; ++i)
      outside = outside || target(i) < current(0) || target(i) > current(size_current - 1);
    if (outside)
      std::cerr << "lin_interp_between_grids: Value is out of range. Extrapolating...\n";
  }

  Kokkos::parallel_for(
      "lin_interp_between_grids", Kokkos::RangePolicy<Execution>(0, size_interpolate),
      KOKKOS_LAMBDA(const int i) {
        bool calc_done = false;
        // Find nearest source level. If no level is above/below, use the endpoint.
        int k = 0;
        while (!calc_done && k < size_current) {
          if (Kokkos::abs(interpolate_altitudes(i) - current_altitudes(k)) < tol) {
            lin_interp_between_grids(i) = current_values(k);
            calc_done                   = true;
          } else if (interpolate_altitudes(i) < current_altitudes(k)) {
            if (k > 0)
              lin_interp_between_grids(i) =
                  lin_interpolate_two_points(interpolate_altitudes(i), // In
                                             current_altitudes(k), current_altitudes(k - 1),
                                             current_values(k), current_values(k - 1)); // In
            else
              lin_interp_between_grids(i) = current_values(0);
            calc_done = true;
          }
          ++k;
        }
        if (!calc_done && k >= size_current) {
          lin_interp_between_grids(i) = current_values(size_current - 1);
          calc_done                   = true;
        }
      });
}

template <class Out, class Source, class Values>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>
lin_interp_between_grids(ScratchWorkspace& workspace,                         // InOut
                         int size_interpolate, int size_current,              // In
                         Out interpolate_altitudes, Source current_altitudes, // In
                         Values current_values) {
  // Return a caller-frame profile without releasing its reservation on return.
  using Result = Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>;
  auto result  = workspace.view<Result>(Kokkos::WithoutInitializing, "lin_interp_between_grids",
                                        size_interpolate);
  lin_interp_between_grids(workspace, size_interpolate, size_current, interpolate_altitudes,
                           current_altitudes, current_values, result);
  return result;
}

template <class Out, class Source, class Values>
Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>
lin_interp_between_grids(int size_interpolate, int size_current,              // In
                         Out interpolate_altitudes, Source current_altitudes, // In
                         Values current_values) {
  // Preserve an independent owned profile for callers using the original interface.
  using Result = Kokkos::View<real_t*, Kokkos::LayoutLeft, typename Values::memory_space>;
  Result result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "lin_interp_between_grids"),
                size_interpolate);
  ScratchWorkspace workspace;
  auto workspace_scope = workspace.scope();
  struct Completion {
    // The legacy owner is local: complete queued scratch use on every exit.
    ~Completion() { Kokkos::fence("lin_interp_between_grids::owned_workspace"); }
  } completion;
  lin_interp_between_grids(workspace, size_interpolate, size_current, interpolate_altitudes,
                           current_altitudes, current_values, result);
  return result;
}

} // namespace clubb_kokkos
