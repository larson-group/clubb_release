---
name: clubb-porting-fortran
description: Understand Fortran kinds, arrays, derived types, state and C interoperability when Fortran is a port source, target or trusted backend.
---

# Fortran semantics and interoperability

Use for Fortran as source, target or backend with the
[porting strategy](../../SKILL.md) and [fidelity](../../fidelity/SKILL.md).
These are language considerations, not a list of any implementation's options.

## Recover the actual source semantics

- Inspect kind definitions, implicit conversions and literal types. A default
  real literal or intrinsic argument can differ from the field's real kind.
  Preserve integer-versus-real exponent behavior, rounding boundaries and
  source operand association when the numerical contract requires them.
- Array lower bounds are not necessarily one; include negative/zero bounds,
  sections, halos, inclusive endpoints and assumed-shape dummy bounds. Separate
  storage positions from physical levels, counts, option codes and serialized
  IDs. Use target-native storage indices and convert only at true boundaries.
- The first array index is contiguous in ordinary Fortran storage. Read actual
  sections/strides and compiler-created temporaries before assuming a zero-copy
  interface. A target layout/permutation should follow access needs rather than
  mechanically retaining this ordering.
- `intent`, optional presence, allocatable/pointer attributes and derived-type
  assignment determine ownership and mutation. Trace allocation/deallocation,
  pointer association, deep field copies and callers' expectations.
- Module state and declaration-initialized locals can persist between calls.
  A local initializer can imply SAVE. Preserve intended case/process lifetime
  and reset behavior without accidentally carrying one case into the next.
- Read active preprocessor/build gates, short-circuit assumptions and error
  flow. Logical expressions do not generally promise short-circuit evaluation;
  do not rely on an unselected operand being safe without checking the code.

## Make a deliberate C boundary

Use `iso_c_binding`, interoperable scalar kinds and `bind(C)` where appropriate.
General Fortran derived types with allocatable components are not simply C
structs. Choose an interoperable descriptor/handle or fieldwise accessor/copy
contract and verify precision, extents, layout, optional presence and mutation.
Do not expose private compiler descriptors as a portable ABI.

Specify value versus reference arguments, string encoding/length/termination,
array order/strides and who owns each allocation. C pointers are not automatically
valid Fortran arrays; use a shape-aware association with a documented lifetime.
Avoid a hidden singleton store when call-owned conversion suffices. Keep any
necessary opaque backend state behind the designated bridge owner.

During the early type milestone, exercise nested/allocated/unallocated fields,
returned/inout data, errors and repeated cleanup/reset as needed by the scope.
Keep statistics state opaque until the final statistics milestone. Remove
temporary accessors/converters when their only migration consumer disappears.

For Python-facing F2PY behavior, also read
[Python and F2PY](../python/SKILL.md).
