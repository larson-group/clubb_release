# Saturation arithmetic and source distinctions

No forced FMA, reassociation, alternative transcendental approximation or other
speculative BFB workaround was introduced. The draft has no compiled or numerical
evidence. Root must review any later arithmetic change against the live native
reference and log the exact expression/evidence in `BFB_ADAPTATIONS.md`.

Ordinary translation choices retained for review:

- Native `T_in_C**2` becomes `T_in_C*T_in_C`. The Flatau factored product and ice
  Horner expression preserve the native operation tree and coefficient scaling.
  Ice coefficients remain `real_t(100)*real_t(source_literal)`.
- Non-integer/nested powers use `Kokkos::pow`; exp/log/log10 use the matching
  Kokkos functions. No alternative base-10 exponent or cube-root algebra is used.
  CPU/GPU libraries and contraction may still differ; no BFB claim follows from
  matching source expression trees or using an unoptimized build.
- Native ordered two-operand MIN/MAX map to ordinary comparisons and ternaries;
  the second operand is selected on a tie. This is not a measured cross-compiler
  claim for NaNs. Adjustment saturation deficits are evaluated once into a local
  `excess` before the comparison, avoiding repeated function evaluation.
- Lookup literals are unchanged and halfway rounding remains away from zero.
  NaN/out-of-integer-range temperature conversion has no new defined behavior.

Inherited source behavior deliberately retained:

- Liquid mixing-ratio GFDL uses a 173.15 K lower bound at native source lines
  149/365 despite a nearby -70 deg_C comment. The separate public GFDL vapor-
  pressure helper uses 203.15 K at line 791. Neither value was harmonized.
- Scalar liquid specific humidity is `ep*(esat/(p-(1-ep)*esat))` at lines
  264–265. Liquid 2D instead has `ep*(esat/(p-(1-ep*esat)))` at lines 529–530.
  Ice uses the scalar liquid form. `I_sat_sphum` is currently false in native;
  the two liquid branches are preserved independently and remain untested.
- Liquid partial-range results are undefined outside the selected levels.
  `WithoutInitializing` output does not invent values there; no BFB expectation
  applies to those elements. This is source-definedness, not a numerical fix.
- Native ice has no lookup selector case. Lookup and invalid ice selectors
  retain the original -99999.999 vapor-pressure sentinel rather than throwing.
- Adjustment's native counted loop leaves `iteration=itermax+1` on natural
  exhaustion, while an EXIT on the final iteration leaves `iteration=itermax`.
  Native tests only equality with itermax at line 1172. The draft preserves this
  exact distinction; it does not repair that potentially surprising condition.

The device table getter, scalar ice extraction, loop fusion, memory staging and
device abort are classified in `INTEGRATION.md`. They are language/execution
adaptations, not demonstrated BFB workarounds or measured optimizations.
