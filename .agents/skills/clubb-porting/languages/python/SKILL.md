---
name: clubb-porting-python
description: Map Python/NumPy objects, array ownership, dtype and F2PY boundaries when porting between runtimes.
---

# Python, NumPy and binding boundaries

Use for Python as source/target or a bridge with the
[porting strategy](../../SKILL.md). For JAX, also read its
[transformation guidance](../jax/SKILL.md); NumPy's eager execution and mutability
are not substitutes for JAX semantics.

## Objects, arrays and state

Python assignment usually aliases an object; it does not copy its fields or
array data. Distinguish a view, shallow object copy and independent snapshot.
Choose explicit owned model state and case reset/lifetime instead of relying on
module globals. Preserve the source's inout/return contract and public argument
order while making target mutation or returned replacements clear.

Specify dtypes and inspect scalar promotion: Python integers/reals, NumPy scalar
types and compiled-language literals can choose different arithmetic. Normalize
at a real API boundary, not through repeated per-routine `_f64`/`_np64` helpers.
Use native storage indices, and map identifiers/serialized ordinals separately.

NumPy arrays have explicit strides; C/F order, transposes and slices affect
contiguity and may trigger binding copies. Decide logical axes and access order
before selecting a layout. Don't repeatedly transpose/copy at every call merely
to preserve the source's storage representation. Return shape and optional-value
contracts should be stable, not a set of forms accepted by local fallback code.

Use array operations for independent large work where suitable; preserve true
recurrences and reduction order when required. Object-per-cell representations,
Python scalar loops and repeated temporary allocations can dominate execution.
Avoid converting each field element into a Python object during hot calls.

## F2PY and compiled backends

Inspect the generated wrapper and actual callable signature. Intent-dependent
returns, hidden dimensions, optional arguments, dtype/contiguity checks and
copy-in/copy-out behavior must be understood rather than guessed from the
original routine declaration. A successful build alone does not prove ABI or
runtime compatibility with the selected interpreter and numerical libraries.

Ordinary F2PY wrappers may not directly expose a general derived type. A bridge
can use fieldwise accessors or an explicit backend-owned handle/store: copy
fields in, call, and copy mutated/returned fields out. If a wrapper requires a
shared store, isolate it, define serialization/reentrancy and lifetime, and
avoid leaking it into every physics routine. Do not choose global storage by
habit when the binding can support per-call ownership.

Keep generated bindings behind one public wrapper owner. It owns dtype/layout,
type conversion, errors and one result contract. Repair missing arguments there;
do not scatter raw-binding imports or adapters for several possible returns
through the port. Check actual mutation/copy behavior and repeated requests at
the early interoperability milestone. Statistics can use a reference-owned
opaque handle until its coherent final replacement.

Use the project's existing dependency launcher/bootstrap and shared run tools.
A future compiled/runtime backend should not require a second case-generation
stack or a hard-coded path to another environment.
