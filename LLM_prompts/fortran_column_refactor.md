# Refactor Fortran column loops

Use when moving column loops into routines or changing multi-column interfaces.
Also follow [routine formatting](fortran_routine_formatting.md). Read the full
routine and callers, and agree which changes must remain BFB before editing.

Move the actual computation and column loop to the intended routine, preserving
its existing meaningful name where practical. A thin wrapper around an old
single-column routine does not complete a requested column-loop refactor.
Keep call sites compact and pass the full array with the intended dimensions;
audit argument order, optional groups and ownership through every caller.

Trace each array's column and level dimensions, including stats, diagnostics,
surface quantities, packing, temporaries, PDF state and microphysics paths.
Where semantics allow, put the contiguous column index `i` in the inner loop;
do not change arithmetic/reduction order merely to reach that layout. Shared
level-only state is not automatically per-column state. Use whole-column rank-2
profile stats and column-vector surface stats where appropriate; remove duplicate
caller updates when ownership moves. Do not retain copy-only legacy objects or
generic indexed buffers merely to keep a scalar implementation underneath.
Keep a tight column loop at an imported scalar-scheme boundary when that core
is outside the requested scope.

Keep behavior-preserving interface/layout work separate from intentional
physics or indexing fixes. In particular, changing a historical `(1,k)` to
`(i,k)` can change results: a task may explicitly preserve it with a TODO for
BFB, then fix it in a separately approved change. Record the boundary rather
than treating every hard-coded first-column index as an authorized repair.

At CPU/GPU boundaries, identify the last writer and next consumers for each
field, and synchronize only where required; stale host/device copies can overwrite
new values. Explain retained copies and record unavailable-device validation.

Validate each meaningful stage with the relevant build/BFB/case checks, including
matched `-O0` outputs for a BFB loop push. Use exercising cases, then
check real multi-column behavior and stats shapes. Unit tests or a new wrapper
alone do not establish equivalence. Inspect existing regression controls before
claiming a new mismatch; report old failures separately. Keep the diff focused,
retain source comments/formatting, and document necessary temporaries briefly.
