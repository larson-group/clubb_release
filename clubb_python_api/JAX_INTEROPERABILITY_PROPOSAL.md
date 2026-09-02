# Tentative JAX Derived-Type Interoperability Proposal

## Status

This document describes a possible debugging-oriented extension to
`clubb_python_api`. It is a tentative design, not an implementation plan that
has been approved. The details should be reviewed and discussed before code is
changed.

## Problem and Goal

The JAX CLUBB implementation owns JAX pytree versions of Fortran derived types,
such as the grid, PDF parameters, configuration flags, and error information.
The Python API owns Python/NumPy mirrors of those same types and already converts
between its objects and derived types held in Fortran storage.

For occasional debugging, it is useful to replace one JAX routine with the
corresponding Fortran routine exposed through F2PY. The desired call site should
remain simple:

```python
output_1, pdf_params = function_f2py_version(input_1, pdf_params)
```

If `pdf_params` is a normal Python API object, the call should behave exactly as
it does today. If it is a canonical JAX object, the API should push its fields to
Fortran and reconstruct a canonical JAX object from the returned Fortran state.

This is intended as a rare, manually edited debugging aid. The immediate goal is
not a runtime configuration system that can switch arbitrary routines between
JAX and Fortran.

## Division of Responsibilities

`clubb_python_api` should remain the sole broker between Python-family callers
and Fortran. JAX core modules should not import F2PY, Python API derived types, or
conversion utilities. The API may accept structurally compatible caller-owned
objects, but all conversion to and from Fortran belongs at the Python API
boundary.

Beyond that separation, there are no hard architectural requirements. The first
version should favor a small, understandable implementation that can be edited
manually while investigating a routine. General registration, dynamic dispatch,
or a configurable mixed JAX/Fortran execution framework should be added only if
real debugging use demonstrates a need.

## Tentative Design

The existing `set_fortran_*` converters already read object fields and push
basic values into Fortran storage. They can continue accepting the current
Python API objects. Where their field requirements match, they can also accept
canonical JAX objects without first constructing an API-specific mirror.

The corresponding `get_fortran_*` converters could accept an optional `like`
template:

```python
def get_fortran_pdf_params(*, like=None):
    values = _read_pdf_params_from_fortran()
    if like is None:
        return pdf_parameter(**values)
    return _restore_like(like, values)
```

Public wrappers would explicitly carry the input template through the call:

```python
def function_f2py_version(input_1, pdf_params):
    set_fortran_pdf_params(pdf_params)
    output_1 = clubb_f2py.function(input_1)
    return output_1, get_fortran_pdf_params(like=pdf_params)
```

When `like` is omitted, behavior remains unchanged and the converter returns the
normal Python API type. When supplied, `_restore_like` could:

1. Rebuild the same object type with `_replace` or its constructor, without
   importing `clubb_jax` into the Python API.
2. Convert returned array fields through the corresponding template field's
   `__array_namespace__()`. NumPy inputs therefore remain NumPy, while JAX inputs
   return JAX arrays.
3. Preserve template-only fields that have no F2PY representation. The current
   example is `ErrInfo.reason_code` on the JAX side.
4. Reject incompatible required fields explicitly instead of silently dropping
   data.

Templates should be passed explicitly rather than remembered in converter
module state. Explicit templates avoid ambiguity when a routine carries more
than one instance of a type and avoid adding hidden state to an API that already
coordinates with global Fortran storage.

## Expected Limits

- F2PY cannot execute during JAX tracing. A substituted routine would be used in
  an eager debugging run, normally with `jax.disable_jit()`.
- Derived-type restoration does not by itself reconcile routine signatures.
  JAX-only stats arguments, additional diagnostics, and differing return order
  may still require small routine-specific wrapper edits.
- Numeric arrays returned outside a derived object may also need restoration to
  the caller's array namespace if a particular debugging path requires strict
  JAX output types.
- Objects created by a Python API routine without a caller-provided template
  should continue to use the normal Python API type unless a concrete use case
  motivates an explicit backend argument.
- This mechanism is for debugging, not differentiable execution through
  Fortran.

## Review Before Implementation

Before implementing this proposal, confirm the first routine or small set of
routines that must be substitutable. Then inspect their exact JAX and Python API
signatures and return contracts. The initial implementation should cover only
the derived types and output arrays those routines actually use.

The review should also decide whether `__array_namespace__` is sufficient for
all required fields, how `ErrInfo.reason_code` should be preserved or reset, and
whether any routine needs a deliberately documented return-order adapter.

Validation should start with round-trip tests for each affected derived type,
followed by a direct JAX-versus-F2PY routine comparison and one short integrated
eager debugging run. No automatic mixed-backend driver should be introduced as
part of that first implementation.
