---
name: clubb-porting-validation
description: Establish a port's numerical oracle and acceptance contract, use integrated CLUBB comparisons and investigate drift without redundant leaf-test suites.
---

# Numerical and capability acceptance

Use with [the porting strategy](../SKILL.md). The source, numerical oracle,
target runtime/device and required capabilities are explicit planning choices.
Do not assume the source or reference must be Fortran. Existing implementation
criteria remain in force when continuing a port unless the user accepts a change.

## Establish the contract early

Record both revisions, build freshness, compiler/runtime/library versions,
optimization and math flags, dtype/precision, device/backend, inputs, cases,
steps/timesteps, overrides, randomness and statistics. Distinguish structural
source authority from the executable numerical reference. Verify matching input
and configuration before diagnosing output.

Try BFB when practical: it makes integration especially decisive. A debug or
low-optimization build can help, but neither a language pair nor `-O0` guarantees
exactness. Likewise Python/JAX are not language-wide proofs of impossibility.
Probe representative operations and the actual compiled execution path before
choosing a feasible contract. Preserve an already required BFB target.

Where exact agreement is infeasible, justify absolute/relative/ULP criteria,
units, near-zero behavior, comparison windows and any case-local exceptions.
Read the maintained comparator's actual semantics (for example percentage
points versus fractional relative error). Existing JAX settings are evidence
for that implementation, not universal defaults for a new port or device.
Check schemas, missing/extra fields, shapes, masks, health/return status and
NaNs separately; numerical tolerances must not hide structural output errors.

## Use real CLUBB execution as the primary oracle

Reuse existing comparison runners and shared `run_scm`, namelist, output and
bindiff owners. For a compatible new executable, prefer `run_scm.py -exe`.
If a thin new adapter is needed, put it at the existing owner rather than
duplicating case preparation or the runner stack under the port.

Use short, active multicolumn cases for iteration. Periodically broaden to
standard/all relevant statistics and supported nondefault flags, scalar/batch
shapes, forcing/cloud onset, restart/reset and longer evolution as required.
Final checks use the agreed duration/settings; a curated uncapped case can still
contain its own cap. Do not universalize a prior milestone's local step limit.

Confirm the changed algorithm actually runs: inspect active options, registered
fields and nonzero physical activity. Filled/masked output can hide real zeros;
a default case may bypass a new branch entirely. Label source-reviewed but
runtime-unexercised branches honestly. Invalid input rejected by both models is
negative-path evidence, not a science pass.

## Avoid unnecessary unit harnesses

An integrated comparison of an active ordinary numerical path is usually the
right test. For example, a custom LU solver with no uncovered options does not
need a new isolated numerical suite if plugging it into the exercised model
would reveal its relevant errors. Read its callers rather than assuming reach.

Use small focused checks for rare boundaries, unexercised supported options,
type/ownership/shape risks or defects hidden by ordinary outputs. Prefer a live
trusted-source oracle over copied formulas or large frozen fixtures. Draft
workers need reasonable confidence and clear risks, not a new proof framework.
Source review remains required; do not claim that every conceivable numerical
bug must manifest in the default model run. Follow the [pytest owner](../../clubb-pytests/SKILL.md) for placement and
human-review rules when adding a warranted test.

## Investigate every new failure

First separate setup, stale binaries/objects, ABI, input/configuration and
schema failures from arithmetic. Compare the starting revision under equivalent
settings before attributing a regression. Narrow a failing integration group
and locate the first diverging step/field/routine, then replay the same operands
to distinguish a formula, state, mask, index or call-order bug from rounding.

Small final differences alone do not prove harmless drift. Examine typed
operands, precision, promotion, reductions, FMA/contraction, math libraries and
RNG mapping. Controlled optimization/precision experiments or local instruction
inspection can establish cause. Record deliberate arithmetic adaptations with
their evidence at the owning boundary; do not apply global workarounds blindly.

In a toleranced port, later correct pieces can push accumulated drift over an
earlier limit. Reinvestigate each new failure rather than assuming it has the
same cause. Explain the first perturbation, propagation and ruled-out defects
before proposing a comparison adjustment. Obtain acceptance for a changed
contract, and preserve unaffected checks. Do not silently increase tolerances,
shorten the final run, disable statistics or alter physics to obtain a pass.

Use [numerical diagnosis](../../clubb-numerical-diagnosis/SKILL.md) for deeper
investigation. Its Fortran examples apply when Fortran is the selected oracle;
they do not override this port's explicitly chosen reference.

## Test the objective as the graph grows

For a device goal, compile and execute real target-owned kernels on that device
from the first meaningful pieces; CPU success does not prove device support.
Agree whether cross-device numerical agreement is an acceptance gate or a
separate checkpoint. Do not inherit another port's GPU policy implicitly.

For JIT/gradients or other transformations, run representative compiled and
feature probes incrementally on the normal path. Check meaningful derivatives,
not just that a call returns; use suitable controlled derivative checks where
needed. A bridge to untransformed legacy code cannot prove whole-model feature
support. JIT-friendly work can reduce test time even when speed is not the goal.
Deeper tuning still follows the requested scope and measured evidence.

## Close out with scoped evidence

Finish the [source audit](../fidelity/SKILL.md) separately. Identify tested
sources/binaries and effective options; distinguish full-duration, short, direct
edge, gradient and device results. Source review, token equality, compilation
and a numerical pass establish different facts. Mixed frozen-lane evidence
needs a real compatibility review plus relevant combined-model checks.

Use [testing/reporting](../../clubb-testing-and-reporting/SKILL.md) for final
evidence. Optimized workload/device measurements support performance claims;
debug passes, allocation/request counts or one contended trace alone do not
establish general speedups. Keep implementation-specific thresholds and command
examples with the maintained implementation owner.

If validation was not run for an audit-only request, state why and do not
claim runtime equivalence.
