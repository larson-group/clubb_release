This draft implements all five Native calendar routines and both immutable month
tables in source order. The complete 254-line Native file and actual JAX mirror
were read before the source outline. Gregorian/Julian integer divisions truncate
toward zero as in Native; elapsed days use FLOOR before integer conversion,
including negative elapsed time. Calendar seconds use `double`, matching Native
`time_precision` independently of core single/double precision.

The conversion assignments, leap predicate, one-based month prefix, leap-day
correction and upper-day fatal check retain the Native expressions and order.
The integer SUM is an ordered loop over the preceding months, followed by the
original day addition. The exact Native fatal literal raises a C++ runtime error;
JAX's omitted fatal branch is not imported. No additional date validation,
integer widening, clipping or fallback is introduced. Month indices and integer
ranges outside the defined Native contract remain outside the supported contract.

Both public tables are immutable constants in their original owner; all mutable
values are call-local scalars. The C++ interface retains Native procedure names,
arguments and output references. It needs only the existing `sec_per_day` value
and standard host math, with no new scientific bridge, global cache or class.
The actual Native callers are radiation solar geometry, input NetCDF date checks,
stats date formatting and existing API wrappers. Those packages continue to use
Native calendar until their call boundaries are ported; this source draft does
not claim those transitive paths are already C++.

At promotion, copy the same-named `.cpp/.hpp` into the core owner and register the
translation unit there. `SOURCE_AUDIT.json` records all five contracts and ordered
assignment expression trees, constants and error behavior. `CLEANUP_AUDIT.json`
records a separate argument/comment/whitespace pass with exact executable-token
equality. No deliberate BFB arithmetic adaptation was added. Compilation,
independent review and all date/rollover/error runtime evidence remain pending.
