// Convert liquid water potential temperature to absolute temperature and back.
//
// thlm2T_in_K_api selects the native scalar, vertical-profile or column/profile
// formula. T_in_K2thlm_api is elemental: its scalar form is device callable,
// and its View overload applies the same expression to each element. Temperatures
// are in kelvin, exner is dimensionless and rcm is liquid mixing ratio [kg/kg].
// Ice is excluded. Arrays use column-first LayoutLeft in the execution space.
// Reference: Cotton and Anthes (1989), Storm and Cloud Dynamics, Eq. (2.51).
#pragma once

#include "array_types.hpp"
#include "constants_clubb.hpp"

namespace clubb_kokkos {

KOKKOS_INLINE_FUNCTION real_t thlm2T_in_K_k(real_t thlm, real_t exner, real_t rcm) {
  // Calculate absolute temperature from liquid water potential temperature.
  // Does not include ice. Cotton and Anthes (1989), Eq. (2.51).
  real_t T_in_K; // Result temperature [K]

  // ---------------------- Begin Code ----------------------
  T_in_K = thlm * exner + Lv * rcm / Cp;

  return T_in_K;
}

inline Real1D thlm2T_in_K_1D(int nz, Real1D thlm, Real1D exner, Real1D rcm) {
  // Calculate absolute temperature [K] on all vertical levels, excluding ice.
  // Cotton and Anthes (1989), Eq. (2.51).
  Real1D T_in_K(Kokkos::view_alloc(Kokkos::WithoutInitializing, "thlm2T_in_K_1D::T_in_K"), nz);

  // ---------------------- Begin Code ----------------------
  Kokkos::parallel_for(
      "thlm2T_in_K_1D", Kokkos::RangePolicy<>(0, nz),
      KOKKOS_LAMBDA(const int k) { T_in_K(k) = thlm(k) * exner(k) + Lv * rcm(k) / Cp; });

  return T_in_K;
}

inline void thlm2T_in_K_2D(int nz, int ngrdcol, Real2D thlm, Real2D exner, Real2D rcm,
                           Real2D T_in_K) {
  // Calculate absolute temperature [K] on all columns and levels, excluding ice.
  // Cotton and Anthes (1989), Eq. (2.51).

  // ---- Begin Code ----
  Kokkos::parallel_for(
      "thlm2T_in_K_2D", Kokkos::RangePolicy<>(0, ngrdcol * nz), KOKKOS_LAMBDA(const int cell) {
        const int i  = cell % ngrdcol;
        const int k  = cell / ngrdcol;
        T_in_K(i, k) = thlm(i, k) * exner(i, k) + Lv * rcm(i, k) / Cp;
      });

  return;
}

inline Real2D thlm2T_in_K_2D(int nz, int ngrdcol, Real2D thlm, Real2D exner, Real2D rcm) {
  // Return an independently owned temperature field through the same kernel.
  Real2D T_in_K(Kokkos::view_alloc(Kokkos::WithoutInitializing, "thlm2T_in_K_2D::T_in_K"), ngrdcol,
                nz);
  thlm2T_in_K_2D(nz, ngrdcol, thlm, exner, rcm, T_in_K);
  return T_in_K;
}

KOKKOS_INLINE_FUNCTION real_t T_in_K2thlm_api(real_t T_in_K, real_t exner, real_t rcm) {
  // Calculate liquid water potential temperature [K] from absolute temperature.
  real_t thlm; // Liquid potential temperature [K]

  // ---- Begin Code ----
  thlm = (T_in_K - Lv / Cp * rcm) / exner;

  return thlm;
}

// C++ overloads implement the native generic/elemental interfaces. The named
// routines above retain their source order; these overloads only dispatch ranks.
KOKKOS_INLINE_FUNCTION real_t thlm2T_in_K_api(real_t thlm, real_t exner, real_t rcm) {
  // Select the scalar native formula.
  return thlm2T_in_K_k(thlm, exner, rcm);
}

inline Real1D thlm2T_in_K_api(int nz, Real1D thlm, Real1D exner, Real1D rcm) {
  // Select the native vertical-profile formula.
  return thlm2T_in_K_1D(nz, thlm, exner, rcm);
}

inline Real2D thlm2T_in_K_api(int nz, int ngrdcol, Real2D thlm, Real2D exner, Real2D rcm) {
  // Select the native column/profile formula.
  return thlm2T_in_K_2D(nz, ngrdcol, thlm, exner, rcm);
}

template <class View> View T_in_K2thlm_api(View T_in_K, View exner, View rcm) {
  // Apply the native elemental formula to a rank-one or rank-two device View.
  static_assert(View::rank == 1 || View::rank == 2);
  static_assert(Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace,
                                           typename View::memory_space>::accessible,
                "Temperature inputs must be execution-space Views");
  View thlm;
  if constexpr (View::rank == 1)
    thlm = View(Kokkos::view_alloc(Kokkos::WithoutInitializing, "T_in_K2thlm_api::thlm"),
                T_in_K.extent(0));
  else
    thlm = View(Kokkos::view_alloc(Kokkos::WithoutInitializing, "T_in_K2thlm_api::thlm"),
                T_in_K.extent(0), T_in_K.extent(1));

  // Select the rank before launching: CUDA extended lambdas cannot first-capture
  // variables solely inside an if-constexpr branch.
  if constexpr (View::rank == 1) {
    Kokkos::parallel_for(
        "T_in_K2thlm_api::profile", Kokkos::RangePolicy<>(0, T_in_K.size()),
        KOKKOS_LAMBDA(const int cell) {
          thlm(cell) = T_in_K2thlm_api(T_in_K(cell), exner(cell), rcm(cell));
        });
  } else {
    const int ngrdcol = T_in_K.extent(0);
    Kokkos::parallel_for(
        "T_in_K2thlm_api::columns", Kokkos::RangePolicy<>(0, T_in_K.size()),
        KOKKOS_LAMBDA(const int cell) {
          const int i = cell % ngrdcol;
          const int k = cell / ngrdcol;
          thlm(i, k)  = T_in_K2thlm_api(T_in_K(i, k), exner(i, k), rcm(i, k));
        });
  }

  return thlm;
}

inline void thlm2T_in_K_api(int nz, int ngrdcol, Real2D thlm, Real2D exner, Real2D rcm,
                            Real2D T_in_K) {
  // Fill the caller's temperature field without a returned-array copy.
  thlm2T_in_K_2D(nz, ngrdcol, thlm, exner, rcm, T_in_K);
}
} // namespace clubb_kokkos
