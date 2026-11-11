# Parameter dependencies

Existing permanent owners provide `nu_vertical_res_dep`, `grid`,
`err_info_type`, `Real1D`, the shared array staging functions, the complete
zero-based `parameter_indices` catalog, `one`, `one_half`, `zero`, `eps` and
`clubb_fatal_error`. No new physical constant, type schema, state argument,
mutable global or class is required.

The one remaining lower boundary is genuine namelist I/O. The staged insertion
extends the existing `generate_driver_io.py` reader producer and uses its wire
contracts. Inputs are `runfile`, `ngrdcol`, required `iunit`, and all 105
initialized rank-1 real arrays in Native namelist order; those arrays are InOut
and `read_status` is Out. The service uses a call-local unit and no Native
parameter science/packing call. Its generated source previews must be merged by
the existing producer, not compiled as competing owners.

Promotion also needs generated-type region preservation, the permanent
initializer translation unit and both existing `clubb_api_module` wrappers.
Switch driver init/two derived calls and numerical_check's parameter validation
call, then remove their obsolete whole-procedure bridge manifest entries.

Pending direct-API limits match Native: missing/partial nu allocation is not
repaired, invalid grid type without prescribed spacing leaves an undefined
average, and valid parameter/profile dimensions are required. Resource errors,
file failures, text spacing and all numerical branches are untested in this
isolated source draft. Source checks are not runtime/BFB evidence.
