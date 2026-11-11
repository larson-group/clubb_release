# Arithmetic and verification notes

This isolated draft has no compilation, runtime or bit-for-bit evidence.
Root must pass CPU Fortran versus CPU Kokkos with zero tolerance before
integration advances. GPU comparisons remain occasional milestone evidence.

No explicit FMA, forced rounding, arithmetic simplification or compiler-specific
correction was added. Every velocity/spacing/weight product and grouped weight
difference follows Native. No diagonal is reconstructed from other diagonals,
and the sign of the Native zero-flow comparison is unchanged.

The stencil is independent by cell and writes a complete bands-first operator.
Native numeric boundary block order and weight category indexing are retained.
There are no reductions, vertical recurrences, new lower calls or statistics
state. The public array entry overloads change storage ownership/layout only
through the existing array owner; they never transpose Native bands.

Source review covers centered and upwind thermodynamic branches, both velocity
signs, both numeric boundaries, and the centered momentum operator. This does
not establish runtime branch coverage or CPU/GPU numerical agreement. The
maintained whole-model comparison is the integration and debugging gate.
