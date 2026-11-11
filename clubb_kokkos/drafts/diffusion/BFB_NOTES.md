# Arithmetic and verification notes

No compilation or numerical execution was performed for this isolated draft.
Source checks do not establish bit-for-bit agreement. Root must pass the regular
CPU Fortran versus CPU Kokkos zero-tolerance gate; GPU comparisons remain
occasional milestone evidence with any differences documented separately.

No FMA, forced rounding, compiler workaround or speculative reassociation was
introduced. Density products, main-diagonal grouped sums, signs, numeric boundary
indices and explicit additions of zero retain Native expression structure.
The ordinary MIN/MAX conditional lowering uses the original operand order and
second operand on ties; no additional signed-zero adaptation is asserted.

The momentum derivative/upwind scratch still executes in the Native option-off
path. This source-definedness/work distinction and JAX's different factoring
are described in [INTEGRATION.md](INTEGRATION.md). No floating-point reduction,
vertical recurrence or inferred off-diagonal conservation replacement is added.

The two existing grid derivative calls remain explicit public C++ API calls.
Their results are owned execution-space Views copied into call-local scratch.
All output bands remain Native ordered at this boundary. Later solver-layout
conversion and its numerical effects belong to the reviewed parent/solver
owners, not these diffusion routines.
