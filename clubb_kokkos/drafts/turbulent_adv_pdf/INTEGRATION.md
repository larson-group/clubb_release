# PDF turbulent-advection draft

This draft fills the implicit band operators and explicit tendencies used by
variance, covariance and flux equations. It implements all four Native routines
in source order; JAX supports the same four. It is isolated from the build.

Promote `turbulent_adv_pdf.cpp` and `.hpp` into their Native-named core directory,
add the source to the core target, and replace the four matching lower calls in
the reviewed moment-advance drafts. The current source has no lower procedure
calls, stats, error state or host diagnostics. It uses existing `grid`, `Real2D`,
`Real3D`, `zero` from `constants_clubb.hpp` and the shared
`array_types` adapters. No new public bridge, selected constant or global owner
is required. Source analysis is recorded in `SOURCE_AUDIT.json`.

Physical fields are `(column,level)` LayoutLeft. LHS output deliberately retains
Native `(band,column,momentum level)`, with super/main/sub bands; pack only at
the existing band-solver boundary. The weight category axis is Native `1:2`
(upper/lower), translated to C++ `0:1`. It is not a signed physical offset.
The two zero endpoints are storage indices `1,nzm`, independently of grid
orientation. Centered and upwind branches use the original storage interior;
`grid_dir` retains the Native directional factors. Godunov RHS sign input is
thermodynamic-level, unlike the ordinary upwind sign input.

Concrete overloads run on the execution space that owns the core grid. Generic
View overloads stage through the existing shared owner and restore the caller's
output allocation, including strided sections; they do not transpose bands.
Only component Views/scalars cross kernel capture boundaries. Native loop
indices are retained locally and converted at View subscripts. Native independent
loops become flattened kernels, with columns adjacent for interior work and
bands adjacent for boundary writes. No recurrence or floating reduction changes.

The complete scientific descriptions and code comments remain in Native order.
Target comments correct the Godunov LHS output description to one rank-3 array
including zero endpoints, remove obsolete optional-range descriptions (no such
arguments exist), correct copied spelling/naming, complete the truncated
upper-zero-boundary comment, and include zero in the upwind else description.
No speculative FMA, explicit contraction, reassociation or signed-zero workaround
is introduced. Ordinary Kokkos min/max follows the literal Native argument order.

Compilation and all numerical evidence are pending. Neither CPU nor GPU, either
grid orientation, centered/upwind/Godunov branches, strided callers nor full model
outputs have been run for this draft. Parent integration owns maintained CPU
Fortran-vs-Kokkos zero-tolerance comparisons and occasional GPU milestones.
