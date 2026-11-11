# Integration contracts

Read `DRAFT_STATUS.md` before using this draft. The complete native module's nine bodies are present; the compact two-routine stage remains available in `incremental_adg1/` for incremental hookup.

For a full-file hookup, copy the draft source/header to their native-named `clubb_kokkos/src/CLUBB_core/` owners, add the C++ source to the existing core library, include its public header at callers, then replace the reached ADG1/ADG2/Luhar native calls with C++ calls. No shared source, CMake or generator edits have been made by this draft. Do not replace working production code before the root's current integrated section is green.

The full draft no longer uses the native ADG1_w_closure bridge: that routine is now C++. To hook only the compact stage first, keep its explicit public16-argument ADG1_w_closure bridge, as documented in `/tmp/clubb_adg1_next/INTEGRATION.md`.

## Missing public complex root adapters

`f2k_dependencies.hpp` declares pending transient C++ adapters; it is not an implementation or fallback. Root should move the boundary into the normal `src/f2k/` owner when wiring it. Current shared bridge generation has no complex primitive support; this draft does not broaden that parser.

The actual lower source is public `src/CLUBB_core/calc_roots.F90`:

```fortran
function cubic_solve(nz, a_coef, b_coef, c_coef, d_coef) result(roots)
  ! nz: integer IN
  ! a_coef,b_coef,c_coef,d_coef: real(core_rknd), shape(nz), IN
  ! roots: complex(core_rknd), shape(nz,3)
end function

function quadratic_solve(nz, a_coef, b_coef, c_coef) result(roots)
  ! nz: integer IN
  ! a_coef,b_coef,c_coef: real(core_rknd), shape(nz), IN
  ! roots: complex(core_rknd), shape(nz,2)
end function
```

All arguments are required, with no optionals or stats/UDT state. `max_cubic_root` invokes these with nz=1 and owned singleton host coefficient Views. The declared result is a column-first HostSpace2D View of `std::complex<real_t>`. A native C shim can return real and imaginary buffers separately and assemble this owned result; do not assume C++/Fortran complex object ABI identity. Nothing needs global native storage.

Existing `f2k::clubb_at_least_debug_level_api` is the remaining shared metadata call. Existing constants, error codes and UDT conversions suffice. The profile HostReal1D alias and three staging overloads are local draft adaptations; they should be reconsidered when complex roots become GPU-native.

## Source review evidence

The final audit verifies all nine definition signatures and order, all six public signatures and three profile adapter argument lists, and every argument of all29 native subroutine calls (6 ADG1,5 ADG2,18 Luhar). It also checks token balance/100-column formatting and records manual formula/branch/comment/ownership review. The single root-function invocation is deliberately evaluated once into a local result before the native MAX operation; likewise the fractional-power candidate is evaluated once before its floor. These local temporaries express native intrinsic-operand evaluation without introducing reusable scratch/global storage.

No runtime evidence is claimed. Root owns compilation, dependency wiring, whole-model comparison, later source cleanup after integration, and permanent BFB-catalog updates when justified by numerical evidence.
