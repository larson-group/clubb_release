# Integration contracts

[`diffusion.hpp`](diffusion.hpp) and [`diffusion.cpp`](diffusion.cpp) implement
both Native public routines in order. `diffusion_zt_lhs` and
`diffusion_zm_lhs` each retain all ten required arguments. Neither introduces
an optional argument, error/status interface or new stats behavior.

Physical fields use column-first LayoutLeft Views. Both Out operators keep
Native `(band,column,level)` shape, with three bands ordered superdiagonal,
main diagonal and subdiagonal. `diffusion_zt_lhs` writes
`(3,ngrdcol,nzt)`; `diffusion_zm_lhs` writes `(3,ngrdcol,nzm)`. Local upwind
operator scratch retains that same bands-first layout. There is no transpose
here: the reviewed moment parents already expect this Native boundary and
transpose only when they reach the established matrix solver API.

All scratch is call local and uses `WithoutInitializing`. Each option writes
the complete output operator. Public generic overloads preserve actual caller
sections using the existing array owner once per entry: compatible Views alias
storage; Out adaptation omits copying old output and restores the caller's
allocation after the call. No axis permutation or alternate physics path is
hidden in those entry overloads.

The actual live grid APIs `ddzm(nzm,nzt,ngrdcol,gr,field)` and
`ddzt(nzm,nzt,ngrdcol,gr,field)` already exist. Their owned execution-space
results are copied into the original Native-named scratch arrays. No additional
Fortran lower bridge, private binding, scalar crossing or fallback is needed.
Existing `zero` is supplied by the production constants owner.

One selected constant is missing from the current generated model-flags owner:
`l_upwind_Kh_dp_term = false`, from Native `model_flags.F90:151`. Add it at that
owner before compiling this draft; do not define a shadow constant in this
physics file. Retain both branches for future Native flag changes.

The source order differs from JAX in ways relevant to fidelity:

- Thermodynamic derivative/upwind scratch is computed only when the option is
  on. Momentum derivative/upwind scratch is computed unconditionally in Native,
  even when the final operator uses the option-off formula. The draft retains
  both behaviors; JAX skips the unused momentum work.
- Main diagonal expressions are evaluated exactly as Native writes them,
  including density factors and their original grouped sum. They are not
  reconstructed by negating computed off-diagonals as in parts of JAX.
- Explicit `zero + lhs_upwind` terms remain present. No zero/sign simplification
  or shared common-factor rewrite is applied.
- Boundary rows use Native numeric indices `1` and `nzt`/`nzm`; grid direction
  and the existing signed spacing operators supply orientation. They are not
  silently exchanged for physical grid bounds.

There are no source recurrences, statistics calls or explicit error paths in
this module. Each cell's band coefficients can execute independently; Native
lower/interior/upper block order and all derivation/conservation comments remain.
Source assumptions, including valid staggered grid dimensions and spacings,
are retained without new guards or numerical repairs.

Root owns production registration, replacing reached `f2k::diffusion_*` calls,
removing the obsolete lower bridges, and CPU bit-for-bit verification before
subsequent integration advances. See [DRAFT_STATUS.md](DRAFT_STATUS.md).
