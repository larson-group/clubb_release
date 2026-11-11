# New-hybrid PDF draft

This directory contains the complete Native-derived draft of
[new_hybrid_pdf_main.F90](../../../src/CLUBB_core/new_hybrid_pdf_main.F90) and
[new_hybrid_pdf.F90](../../../src/CLUBB_core/new_hybrid_pdf.F90): ten routines
across 2453 Native source lines. The draft is source-audited and formatted,
but **has not been compiled, integrated, executed or checked for BFB**.

[new_hybrid_pdf_main.cpp](new_hybrid_pdf_main.cpp) implements the 56-argument
public driver and its two private elemental helpers, in Native order.
[new_hybrid_pdf.hpp](new_hybrid_pdf.hpp) implements all seven public lower
elemental routines as device-callable scalar interfaces. Reference outputs
preserve Native `intent(out)` and `intent(inout)` behavior. No new optional
arguments, stats interface, global storage or scalar Fortran boundary is used.

## Promotion

1. Promote the three Native-named source/header files into
   `clubb_kokkos/src/CLUBB_core/` and add `new_hybrid_pdf_main.cpp` to the
   existing C++ library. The lower header needs no compilation unit.
2. Through the existing `src/f2k/generate_derived_types.py` model-constant owner,
   export `l_explicit_turbulent_adv_wp3` and
   `l_explicit_turbulent_adv_wpxp`. Both are Native compile-time constants with
   value `.false.` at `model_flags.F90:120–122`. The third flag,
   `l_explicit_turbulent_adv_xpyp`, already exists in the current C++ mirror.
   All three branches are implemented; no draft fallback definitions are used.
3. Include the driver header in the existing PDF closure caller and replace
   its reached `f2k::new_hybrid_pdf_driver` call with the same C++ name and
   argument list. Remove the whole-driver generated bridge when its final C++
   caller is replaced. Keep other PDF choices at their current boundaries.
4. Build and use the maintained short multicolumn model comparison for
   acceptance and diagnosis. Root owns promotion, dependencies and runtime
   validation. Nothing in this directory changes production build inputs.

The only lower routine outside this draft is the already-ported
`compute_gamma_Skw`. Existing `constants_clubb.hpp`, `model_flags.hpp`,
`array_types.hpp` and `pdf_parameter_module.hpp` supply the required types and
constants. The caller must supply allocated, conformable implicit-coefficient
Views, as Native initialization does. The draft updates all 27 coefficient/
explicit-term fields used by Native, preserving their current allocations and
metadata. Passive-scalar fields are touched only when `sclr_dim > 0`.

## Preserved behavior and ownership

The driver retains the Native column loop, lower-call order, strict comparisons,
per-scalar copy/work/copy order and all three explicit/implicit advection
branches. Profile calculations use GPU kernels; adjacent independent array
assignments share a kernel without changing operations within a cell. Source-
named scratch Views belong to each call. Model fields remain column-first
`LayoutLeft`; kernels capture component Views rather than the complete UDT.
There is no performance claim for this unmeasured draft.

The invalid-mixture check remains `ANY(mixt_frac < 0)` immediately after the
setter calculation and component-variance assignment for each column. An
integer GPU reduction implements Native `ANY`; the exact two Native diagnostic
strings are printed before a host exception implements Native `ERROR STOP`.
No additional realizability test or error-info argument was added.

The private responder retains both sign-dependent skewness-bound formulas,
ordered clipping, and the subsequent responder call. The inactive exponential
`F_w` blend and its explanations remain visible; the active gamma formulation
is unchanged. The lower routines preserve invalid-mixture sentinels, Native
zero floors, zero-flux/single-Gaussian cases, and every semi-implicit fallback.
The last mixed-moment fallback still tests only `F_w > 0`, exactly as Native.

[source_outline.json](source_outline.json) records the full initial inventory.
[source_audit.json](source_audit.json) records the separate comparison of every
signature, assignment-expression tree, branch predicate and subroutine-call
argument/order. Integer powers are audited with their complete bases and
surrounding groups. Scientific descriptions, full method equations/citations,
scratch ownership and the two fatal diagnostic strings were also reviewed.

JAX supplied scope context only. Its reduced driver signature, omitted fatal
path, formula aliases, gradient guards and clipped square roots were not
imported. Preparation/source-audit scripts remain draft tools and should not
be promoted into the production physics directory.
