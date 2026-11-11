# Grid dependencies

All scientific lower dependencies already exist: `array_types.hpp`, the
canonical `grid` mirror, `err_info_type`, `one`, `clubb_fatal_error`, and the
existing public transient debug query `f2k::clubb_at_least_debug_level_api`.
No new Fortran service, constant, mutable state, type schema or matrix owner is
required. `read_grid_heights` remains the existing explicit Native I/O boundary.

The private weights use the existing zero-based `t_above`, `t_below`, `m_above`
and `m_below` constants. Native physical grid metadata remains one-based; only
actual View indices are zero-based. Generic input staging uses the existing
shared array owner. Single-column API packing is call-local Kokkos storage.

Promotion must preserve generated grid regions and all existing operator bodies,
add the fragments at their stated Native positions, update reached driver/API
calls and retire obsolete bridge entries together. All runtime branches and
resource-failure behavior are untested in this isolated draft.
