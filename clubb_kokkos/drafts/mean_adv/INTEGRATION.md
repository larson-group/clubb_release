# Integration contracts

[`mean_adv.hpp`](mean_adv.hpp) and [`mean_adv.cpp`](mean_adv.cpp) retain both
Native public routines in order. `term_ma_zt_lhs` keeps ten required arguments;
`term_ma_zm_lhs` keeps seven. No optional arguments, error interface, diagnostic
state or extra lower science calls are introduced.

Physical velocities, spacings and weights use column-first LayoutLeft Views.
The weight category is the third dimension, with Native `t_above:t_below` and
`m_above:m_below` bounds `1:2` mapped directly to View categories `0:1`.
The Out operators remain Native `(band,column,level)`: thermodynamic
`(3,ngrdcol,nzt)` and momentum `(3,ngrdcol,nzm)`, with super/main/sub bands.
There is no operator-boundary transpose. The reviewed moment parents already
expect Native bands-first arrays and own the later solver API conversion.

Outputs retain caller allocations. Every branch writes the complete operator;
public generic entry overloads reuse the existing array owner to preserve caller
sections once per call, omitting an old-value copy for the Native Out field.
Compatible execution-space Views alias existing allocations. No hidden axis
permutation or alternate stencil is added.

There are no missing public lower routines or constants. The production
constants owner already supplies `one` and `zero`. Native module-private
`ndiags3=3` remains private to this implementation, and per-routine diagonal and
weight selectors keep their exact Native one-based values.

The complete branch and boundary behavior is preserved:

- Thermodynamic centered differencing runs interior, numeric last boundary,
  then numeric first boundary blocks in Native order. The weight differences
  retain their original signs and expression grouping.
- Thermodynamic upwind differencing retains `grid_dir*wm_zt >= zero` in all
  three blocks. Zero and signed zero select the same source branch; a NaN
  comparison follows the source else branch. Boundary coefficients remain
  exactly those selected by Native, with no extra physical-bound reindexing.
- Momentum differencing remains centered only. All three bands are zero at
  numeric level `1`, then the interior is computed, then all bands are zero
  at numeric level `nzm`. Physical upper/lower meaning reverses with the grid.
- Native equations are evaluated directly, without copying JAX's named common
  factors or adding an upwind momentum option that the source does not expose.

All loops are independent coefficient calculations. There is no recurrence or
shared scientific state, and no call-local physics scratch array is needed.
Source grid/dimension assumptions remain unchanged. Numerical diagrams,
argument descriptions and units remain attached to their Native owner.

Root owns shared source registration, replacing reached `f2k::term_ma_*` calls,
removing their obsolete bridges, and the CPU bit-for-bit comparison before the
next integration. See [DRAFT_STATUS.md](DRAFT_STATUS.md).
