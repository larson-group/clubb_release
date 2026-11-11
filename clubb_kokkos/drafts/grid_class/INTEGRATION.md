# Remaining grid construction

This draft adds `setup_grid`, `cleanup_grid_api`, `setup_grid_heights`,
`calc_zt2zm_weights` and `calc_zm2zt_weights` to the existing grid owner.
The authoritative generated `grid` schema and the 17 existing operator bodies
remain unchanged. `read_grid_heights` stays an explicit public Fortran file
service, including its schema and monotonicity checks, as agreed for this phase.

Merge the fragments into the current live owner; do not replace the whole file:

- Add `procedure_includes.hpp.inc` and `procedure_declarations.hpp.inc` around
  the existing generated type and before existing operators.
- Add the missing includes in `cpp_includes.inc` and the private weight forward
  declarations in `private_declarations.cpp.inc` to the current source preamble.
- Put `prefix_definitions.cpp.inc` before `redirect_interpolated_azm_k`.
- Put `calc_zt2zm_weights.cpp.inc` after `cubic_interpolated_azm_2D`, before
  `linear_interpolated_azt_2D`.
- Put `calc_zm2zt_weights.cpp.inc` after `cubic_interpolated_azt_2D`, before `gradzm_2D`.
- Merge the four overloads in `clubb_api_module.hpp.inc` at their Native wrapper
  positions, with `api_includes.hpp.inc`. Public names follow the Native generic
  interfaces; private `_single_col`/`_multi_col` specific names are rank overloads.
- Replace reached construction/height/cleanup callsites. Remove those three
  obsolete whole-procedure bridge manifest entries and retain `read_grid_heights`.
  The two weight helpers are private C++ bodies, never transient public bridges.
- Update the existing grid opening description/TODO to describe construction,
  signed spacings and weights as C++ behavior and grid-file input as explicit
  Native I/O. Retain its existing interpolation/BFB notes.

Native order and behavior are retained: first-column-only extent selection,
ascending/descending searches, physical metadata assignment after weight calls,
all-column fatal range errors, source debug gates, and the first failing surface
column only. The surface check intentionally uses storage level 1 for either
storage direction. A stretched thermodynamic grid needs at least two levels;
no JAX denominator floor or one-level fallback has been added. Independent
interior and boundary weight writes share a launch with direct divisions and
unchanged expression association.

Kokkos Views release reference-counted storage rather than returning a Fortran
`DEALLOCATE STAT`. Cleanup therefore cannot reproduce the invalid-lifecycle
warning `Grid deallocation failed.`; no invented status is added. Allocation
checks follow Native component order and convert ordinary C++ allocation
exceptions to the Native all-column fatal path. Recovery from a Kokkos runtime
abort is not claimed. Valid caller-owned lifecycle behavior awaits integration.

`SOURCE_OUTLINE.json` records source and actual JAX scope. `source_audit.py`
produces separate `SOURCE_AUDIT.json` and `CLEANUP_AUDIT.json` records, including
exact fragment hashes. These are owner source checks, not independent review or
runtime evidence. No build, model comparison, CPU/GPU or BFB claim is made.
