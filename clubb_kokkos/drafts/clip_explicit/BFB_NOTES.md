# Arithmetic and verification notes

This isolated draft has no compilation, runtime or bit-for-bit evidence. Root
must establish CPU Fortran versus CPU Kokkos equality with zero tolerance before
integration advances. GPU comparisons are occasional milestone evidence, with
any field differences documented independently.

No compiler-specific FMA, forced-rounding correction or reassociation was added.
Integer squares and cubes use ordinary explicit multiplication, retaining source
operand grouping. Native `SIGN` maps to `Kokkos::copysign`, preserving the sign
bit of the source operand. `EPSILON(rtm)` maps to
`std::numeric_limits<real_t>::epsilon()` because Native EPSILON depends on kind,
not the value or local floating spacing. The Native maximum in cloud-water
clipping uses its original operand order and second-operand tie choice.

The covariance bound keeps `max_mag_corr * sqrt(xp2*yp2)` and strict comparison
branches; no JAX floor or unconditional clip is introduced. Both skewness
limiter expressions, the absolute cap, variance comparison bounds, warning
values and normal STOP boundary remain source faithful. Host debug copies only
compare, assign and print variance values; they introduce no alternate floating
physics formula. The invalid-threshold reduction sums integer hits only.

[INTEGRATION.md](INTEGRATION.md) records Native/JAX distinctions and source
limitations. Source review covers all non-default branches, optional maximum,
blank budget names, perturbed winds and both grid directions. It cannot establish
CPU/GPU arithmetic agreement or performance. The maintained whole-model
comparison remains the integration and debugging gate; no separate numerical
oracle suite was run for this draft.
