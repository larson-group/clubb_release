This draft adds the eight missing Native procedures while preserving the four
already-live definitions. It is uncompiled and untested. Copy the reviewed source
into the permanent Native-named owner only when the integration owner is ready
for the CPU Fortran versus CPU Kokkos gate.

The public interfaces and all required arguments follow Native. Private
`check_nan` scalar/View dispatch expresses the Native generic; the two existing
View implementations are unchanged. Do not reduce the detailed PDF/settings
checks to the smaller JAX predicates.

Pending public lower calls are declared, without implementations, in
[`f2k_dependencies.hpp`](f2k_dependencies.hpp), with exact contracts in
[`lower_contracts.json`](lower_contracts.json):

- `error_code::initialize_error_headers()` has no arguments. Native lines 69–94
  use a bare `subroutine initialize_error_headers` declaration. The current
  `fortran_procedures` parser requires parentheses and rejects it. Extend that
  existing owner or provide a shape-free BIND(C) adapter through the existing
  transient owner; do not change Native or invent a new global accessor. This
  call owns Native debug-header metadata, including its OpenMP/process/CAM gates.
- `parameters_tunable::check_parameters_api(ngrdcol,clubb_params,err_info)` has
  integer In `ngrdcol`, real In `clubb_params(ngrdcol,nparams)` and UDT InOut
  `err_info`. It has no optional arguments. The current parser reads this exact
  contract. Keep it before the mixed configuration checks, without an invented
  early return after it.

The only currently missing selected model constant is
`l_explicit_turbulent_adv_wpxp = false`, from Native `model_flags.F90:122`.
Add it through the existing derived-type/constants generator owner. The PDF,
saturation, placement, advancement-order, `Lscale_tau` and LAPACK selectors used
here are already public C++ constants. Remove the staging dependency include
when the exact public transient interfaces exist.

Preserve the `CLUBB_GPU` compile gate on Native LAPACK warnings and the
`l_input_fields` fatal stop. The C++ and Native definitions of this macro must
agree; selecting a CUDA Kokkos package alone is not evidence that the macro is
set. C++ expresses the Native ERROR STOP as `std::runtime_error` at the same host
branch, retaining its exact message. Debug-level queries and stats field gates
reuse the existing public Native metadata calls.

Memory and caller contracts:

- `length_check` and `pdf_closure_check` accept const execution-space
  `LayoutStride` profiles, retaining the actual per-column subviews without
  copying full physical fields. PDF scalar profiles are `(level,scalar)`.
  The complete `pdf_parameter` UDT remains owned by the caller.
- Native PDF-parameter checks explicitly use `(1,:)`; the C++ subviews therefore
  use column zero even when a caller iterates another column. Do not substitute
  the caller's column. Literal field labels keep Native `%` and the inherited
  scalar label `wpsclrtlp`.
- `sfc_varnce_check` keeps all six scalar and three vector inputs required. The
  actual surface caller already stages scalar diagnostics. Its HostSpace
  vectors use the justified public template, which stages through the existing
  array owner once before selecting the concrete execution-space overload.
- `rad_check` has column-first `(column,level)` fields and one per-call vapor
  scratch allocation. `invalid_model_arrays` has strided execution-space
  profiles and `(level,species/scalar)` collections. Its hydrometeor names are
  `(character,species)` with the exact Native width of ten characters; only
  failing name diagnostics stage that list.
- New finite/parameter checks use integer reductions, then host decisions.
  Owning UDTs, Native calls and host printing do not enter device lambdas.
  Complete parameter values stage only inside their failing/warning branches.
- The six module sponge flags come from the existing
  `f2k::driver_environment()` snapshot, after the two original initialization
  and parameter-validation calls. This is existing shared Native metadata,
  not model array storage or a new accessor.

Native diagnostics and state distinctions remain visible. Errors set every
column fatal. The both-false relaxed-clipping branch prints its Native fatal
wording but leaves the commented-out fatal assignment/return inactive. The ten
mixing-length parameter warnings use `> one || < one`, preserving NaN behavior.
The eight host-implementation errors accumulate without early return. The four
unchanged live bodies retain NaN-before-negative return, beginning-step reset
and cold-field warning order. `invalid_model_arrays` reports every bad field
and returns a bool; it does not mutate an error structure.

All diagnostic labels and payload values are retained. C++ full-precision
streams replace Fortran list-directed formatting, so whitespace, boolean and
number rendering are not byte-format claims. Hydrometeor name trimming removes
only trailing spaces. No hidden global scratch, fallback or runtime option is
introduced.

Run only the integration owner's maintained comparison after promotion. This
draft adds no build targets or separate oracle suite; CPU BFB and occasional GPU
milestone evidence remain pending.
