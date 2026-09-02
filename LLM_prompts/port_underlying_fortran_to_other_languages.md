Port underlying Fortran code to another language while preserving the Fortran
file's structure, ordering, names, comments, routine boundaries, and behavioral
surface.

Use this prompt when a user asks to port, re-port, audit, or "similarize" a
Python, JAX, or other target-language file against the underlying Fortran. The
goal is not an idiomatic rewrite in the target language. The goal is a
maintainable mirror: a reviewer should be able to compare the Fortran and target
files block-by-block and see the same routines, calls, comments, and logic in
the same order.

This prompt is specific to CLUBB-style Fortran-to-Python,
Fortran-to-JAX, or Fortran-to-other-language ports. In the rest of this prompt,
"source" means the underlying Fortran file or routine, and "target" means the
ported file or routine in the other language.

Core principle:

- Treat the underlying Fortran file as the authority.
- Keep the target file as close to the Fortran as the target language reasonably
  allows.
- A port is not complete merely because it runs or passes comparison tests.
  Structural mirror quality is part of the acceptance criteria.
- Do not introduce abstractions, helpers, optionals, aliases, or reordered logic
  just because they are convenient in the target language.
- Do not satisfy a port request with a wrapper around the existing Fortran/F2PY
  top-level routine unless the user explicitly asks for a wrapper. A wrapper can
  be useful as a temporary baseline, but it is not a port.
- It is acceptable to call already-existing lower-level kernels while porting a
  large routine, but the boundary must be clear: port the requested routine's own
  control flow, helper routines, clipping, stats, diagnostics, and call ordering.
- Do not bypass the repo's public target-language API to call raw generated
  Fortran/F2PY bindings from a ported target file. If the public API is missing
  needed behavior, fix or extend the API wrapper first.
- When a compliance audit is intentionally conservative, fix clear low-risk
  drift and document unavoidable or out-of-scope conflicts in the target file.
  Do not silently leave known structural drift, reduced diagnostics,
  unsupported branches, target-only helpers, or transitional API boundaries.

Before editing:

1. Identify the canonical underlying Fortran file and the target port file.
   - Read the full Fortran routine(s), not only the failing section.
   - Read the target file enough to understand existing local conventions.
   - If another target-language port exists, use it only as secondary context.
     The canonical Fortran still wins.

2. Build a Fortran outline before changing code.
   - Routine names in order.
   - Public and private helper routines in order.
   - Main control-flow blocks in order.
   - Major comments and divider blocks in order.
   - Routine calls in order.
   - Argument lists, including required versus optional arguments.
   - Output, input/output, and return-value ownership.

   Do this as a real intermediate artifact for substantial ports. Do not skip
   directly to coding based only on the first visible block or a similarly named
   target wrapper.

3. Decide what counts as a necessary language adaptation.
   - Index-base conversion is necessary, for example Fortran 1-based parameter
     indices mapped to Python 0-based array indices.
   - Array-layout conversion may be necessary.
   - Different syntax for inout/out returns may be necessary.
   - A target-language API wrapper may require a small adapter.
   - Document any retained adaptation that prevents a literal mirror with a
     nearby in-code comment or TODO.
   - For substantial ports, keep a short adaptation note in the target file or
     final response that names each intentional deviation. Retained deviations
     should be documented in the target file, not only in the final response.
   - If a deviation is intentionally left in place because fixing it would be
     too risky, too broad, or blocked by unported upstream/downstream code, add a
     concise in-code TODO at the exact boundary where the deviation appears. The
     TODO should name what differs from the source, why it exists now, and what
     future change would remove it.
   - If a deviation cannot be justified concretely, remove it before declaring
     the port complete. If it is fixed immediately, no lingering TODO is needed.
   - Prefer local, actionable TODOs over vague module-wide warnings. For
     example, put a TODO above a target-only helper, unsupported branch, API
     conversion import, or bulk input-conversion block rather than only saying
     the file is "transitional."

4. Audit the target API boundary before porting.
   - List the lower-level routines the source calls or the target must delegate
     to.
   - Confirm each delegated routine is available through the target API, for
     example `clubb_api` for Python/JAX driver ports.
   - If a needed routine is only available through raw F2PY/generated bindings,
     update the API wrapper first. Do not work around the missing API inside the
     ported physics file.
   - Confirm the API wrapper owns low-level dtype, memory-layout, and return
     contract details. Those details should not be reimplemented locally in each
     port.

Routine inventory:

1. The target file should have a 1-to-1 routine inventory with the Fortran file
   for the portion being ported.
   - Do not omit Fortran routines because they are private helpers.
   - Do not create new helpers that do not exist in the Fortran.
   - If the Fortran has a helper routine, port that helper using the same name and
     relative position where possible.
   - If a helper truly cannot be ported because it belongs in a lower-level API,
     state that explicitly and keep the target call site as close as possible.
   - Preserve routine order. Python and many other target languages can call
     routines defined later, so moving the public driver routine to the bottom is
     not an acceptable convenience.
   - New target-only helpers are a defect unless they are explicitly justified as
     unavoidable language adapters. They must not hide field-specific call lists,
     branch-specific behavior, diagnostics, stats, or clipping differences.

2. Preserve routine names.
   - Use the Fortran routine name unless the target language has a hard naming
     restriction.
   - Avoid "compat", "adapter", "impl", or underscored aliases unless they are
     already established in the Fortran or required by the target binding layer.
   - Do not use temporary renamed variables such as `_name` only to unpack and
     then assign back to `name`. Unpack or assign directly to the Fortran name.

3. Apply a strict helper acceptance standard.
   - A helper is acceptable when it mirrors a Fortran helper, isolates a
     deliberately reduced diagnostic block, handles an unavoidable language
     boundary, or removes a large repeated block that is structurally identical
     and explicitly documented.
   - A helper is not acceptable merely because it makes the target language look
     cleaner or more idiomatic.
   - Remove helpers that only assign or reorder fields, call one public API
     function, check an error flag, wrap stats calls, wrap clipping calls, hide
     simple indexing, or perform dtype/layout conversion without a real boundary
     reason.
   - Prefer visible repeated code when it makes the target file easier to
     compare against the Fortran source. Do not abstract away a Fortran mapping
     that is clearer inline.

Argument lists:

1. Match the Fortran argument list as closely as the target language allows.
   - Same argument names.
   - Same argument order.
   - Same grouping on continuation lines when practical.
   - Same required versus optional status.
   - Same public inputs and same returned outputs or inout values.

2. Do not add superfluous optional arguments in the target port.
   - Do not expose local debug flags, order flags, compatibility flags, or dummy
     stats/config arguments unless the Fortran routine actually takes them.
   - If the Fortran uses global state, use the target equivalent of global state
     instead of adding a new routine argument.
   - Optional arguments should appear after required arguments when the Fortran
     convention requires that.
   - Avoid keyword-only separators or other target-language signature features
     that make the port's call surface differ from the Fortran/API it is meant to
     replace.

3. Keep target wrappers honest.
   - If an argument is required in the Fortran, it should not become silently
     optional in the target wrapper.
   - If an argument is not present in the Fortran, it should not appear as an
     output or public input in the target wrapper.
   - If a Fortran derived type or global state is replaced by a transitional
     scalar, object, or side-effect helper, document the drift at the signature
     and at the call site that supplies it.
   - If an audit script exists, use it to catch extra, missing, reordered, or
     falsely optional arguments. Treat this as a quick contract check, not as the
     main porting workflow.

Call sites:

1. Preserve call order exactly unless the Fortran behavior changed or the user
   explicitly asks for a behavior change.
   - Do not move diagnostics, stats calls, clipping, initialization, or fatal
     checks across Fortran block boundaries.
   - If a routine computes a value internally in the Fortran, do not recompute it
     in the port unless the Fortran also does.
   - If Fortran code moved a stats or diagnostic update inside a helper, remove
     duplicate target-side updates.
   - Re-check field-specific arguments even when several fields share one solve
     path. For example, pressure, clipping, stats, and turbulent-advection terms
     may look interchangeable but still differ by field in the Fortran call list.

2. Format target-language calls to mirror the Fortran call formatting.
   - Use positional arguments when the target routine can safely accept them.
   - Group call arguments on lines that correspond to the Fortran continuation
     lines.
   - Keep related Fortran arguments together even if another grouping would be
     more idiomatic in the target language.
   - If a target call truly needs keywords, leave only that call keyworded and
     add a short TODO explaining why positional arguments cannot yet be used.

3. Keep called routine signatures aligned too.
   - If a call cannot be made positional because the called wrapper has drifted,
     fix the called wrapper when that is in scope.
   - If fixing the wrapper is not in scope, identify it as a target API mismatch
     rather than changing the port's structure to hide the mismatch.
   - Do not factor several similar Fortran calls into a new generic helper until
     every field-specific argument has been compared and the helper is documented
     as an unavoidable adaptation. Repetition is preferred when it makes the
     port easier to compare with the Fortran.

4. Do not create local pseudo-APIs.
   - Do not add one-line local helpers that only call `clubb_api`, convert an
     array, check a blank name, or forward arguments.
   - Stats, clipping, solvers, diagnostics, and constants should be called
     through the established target API or shared target abstraction unless the
     source Fortran has an equivalent helper routine. In pure JAX kernels, avoid
     direct Fortran-backed side effects; represent them through the JAX-side
     diagnostic/stat mechanism or document the transitional boundary.
   - If a blank-name or no-op condition is needed for a specific stats call,
     keep that condition explicit at the call site or move the behavior into the
     shared API after confirming it matches the source semantics.
   - Do not add defensive adapters for multiple possible API return shapes. The
     API wrapper must have one stable return contract; if it does not, fix the
     API wrapper.

Comments and block structure:

1. Preserve Fortran comments and divider blocks.
   - Keep comments in the same relative location.
   - Keep major block headers even if the target language code below them looks
     different syntactically.
   - Do not delete diagnostic, error-checking, or setup comments because the
     corresponding target code is shorter.
   - Preserve source module/routine descriptions, references, and important
     warnings in target module docstrings or routine docstrings when a literal
     comment location does not make sense in the target language.
   - Preserve comments that explain call-site constraints, boundary-level
     exclusions, parameter ordering, solver choices, or stats/diagnostic
     semantics. These comments often prevent future off-by-one or wrong-branch
     mistakes even when they do not affect executable code.

2. Remove target-only explanatory comments unless they document a real language
   adaptation.
   - Do not narrate obvious assignments.
   - Do not add comments that explain the porting process inside the code.
   - If a comment exists only because a target-language workaround is needed,
     keep it short and specific.
   - Use TODOs only for known remaining work. A TODO should not be a substitute
     for understanding the source; it should record a concrete source/target
     mismatch that could not be fixed in the current pass.

Calculations and state updates:

1. Translate calculations mechanically.
   - Preserve algebraic order where practical.
   - Preserve clipping, max/min, masking, and conditional update order.
   - Preserve temporary variables when the Fortran uses them for ordering,
     diagnostics, or BFB-sensitive behavior.
   - Do not combine expressions, precompute values, or vectorize in a way that
     changes evaluation order unless the Fortran has already done so.

2. Use the same variable names for the same physical quantities.
   - Avoid local aliases that only rename Fortran variables.
   - Avoid shape aliases such as `shzt` or `shzm`; write dimensions explicitly
     at allocation sites so the port is easy to compare with Fortran dimensions.
   - Avoid defensive dimension changes such as `max(sclr_dim, 1)` unless the
     target runtime requires them. If required, use the smallest workaround and
     document why it exists.

3. Preserve diagnostics and error checks.
   - Do not drop debug-gated checks because the target function lacks a local
     debug argument; use the Fortran-equivalent global/debug mechanism.
   - Preserve beginning-of-routine and end-of-routine validation blocks.
   - Preserve fatal-error early returns in the same positions.
   - Diagnostic routines count as behavior. Do not replace a detailed Fortran
     diagnostic routine with a stub unless the reduced diagnostic surface is
     explicitly documented as incomplete.

4. Double-check branch-specific behavior.
   - Do not assume that matching one default case exercises all branches.
   - Check passive scalar, wind-prediction, perturbed-wind, PDF-type, implicit
     versus explicit turbulent-advection, clipping, and stats branches when they
     exist in the Fortran.
   - If a branch is known inactive for the current driver, document why instead
     of silently omitting it.
   - If a driver rejects a Fortran feature, the port may leave that feature
     unsupported only when the target file or final response says so and points
     to the driver gate or configuration reason.

Preprocessor and configuration handling:

1. Resolve Fortran preprocessor conditionals before porting them.
   - Check compile, CMake, and config code before assuming a macro is active.
   - Omit branches only when they are truly inactive for this target and note the
     reason.
   - Ignore ifdef-only optional arguments when the user explicitly asks the port
     to match the normal API surface without those build-specific extras.

2. Use target-language constants for Fortran constants.
   - Do not pass Fortran compile-time order/debug constants as user-facing
     arguments if the Fortran routine does not.
   - Import or define constants near the top of the file if the Fortran treats
     them as module constants.
   - If a source constant is an index into an array and the target language uses
     a different indexing convention, define the target-native constants in a
     dedicated mirror file or constants section. Do not scatter call-site
     offsets such as `idx - 1` throughout the port unless the offset is a true
     local Fortran-to-target boundary conversion.

API and binding boundaries:

1. Target ports must not import or call raw generated Fortran/F2PY modules
   directly.
   - In Python/JAX ports, imports such as `import clubb_f2py` are defects unless
     the file is itself an API wrapper whose purpose is to expose F2PY safely.
   - If a target port needs a lower-level Fortran-backed operation, call the
     public API wrapper, such as `clubb_api.band_solve`,
     `clubb_api.clip_covar`, or `clubb_api.stats_update`.
   - If the public API lacks a needed argument, output, or option, update the API
     wrapper and its contract rather than adding raw F2PY calls to the port.

2. Dtype and memory-layout conversion belongs at API boundaries.
   - Do not add per-port helpers such as `_f64`, `_np64`, or local
     `asfortranarray` wrappers unless the source Fortran has an equivalent
     routine or the helper is a documented, unavoidable target-language
     adaptation.
   - Prefer API wrappers that convert inputs exactly once before entering a
     Fortran-backed routine.
   - Avoid bulk-converting every routine input at entry. It hides which calls
     truly require NumPy/Fortran buffers and makes future JAX transformations
     harder.

3. API return contracts should be strict.
   - A port should not accept several possible return shapes from the same API
     call.
   - Compatibility code such as `if isinstance(result, tuple)` for an API result
     is a sign that the API wrapper needs to be fixed or documented more tightly.
   - The port should call the API as if its contract is stable, then fail clearly
     if that contract is broken.

JAX-specific porting guidance:

Use this section when the target file is a JAX port or a Python file intended to
become JAX/JIT-capable. The Fortran still controls routine inventory, call order,
comments, names, and behavior. These rules only describe how to express the
Fortran logic in JAX without creating future JIT blockers.

1. Keep internal physics math in JAX.
   - Use `jax.numpy` (`jnp`) for array math, allocation, masking, reductions, and
     elementwise functions inside ported physics routines.
   - Do not import or use normal `numpy` in a JAX physics port except in explicit
     non-physics I/O, test, driver, or transitional API-boundary code.
   - Do not add local NumPy conversion helpers such as `_f64`, `_np64`,
     `_zeros`, or `_as_fortran` inside JAX physics files. If a Fortran-backed API
     call needs NumPy or Fortran-contiguous buffers, the shared API wrapper should
     own that conversion.
   - Avoid bulk-converting every input with either `np.asarray` or
     `jnp.asarray` as a substitute for understanding the data boundary. Normalize
     values where needed, and keep conversion intent clear.

2. Prefer array expressions over Python translations of Fortran mutation.
   - Do not port vertical loops as Python `for k in range(...)` loops that update
     one level at a time with `.at[:, k].set(...)`.
   - Prefer whole-slice expressions such as `arr.at[:, 1:nzm-1].set(...)`,
     strided interleaved updates, `jnp.stack`, `jnp.concatenate`, `jnp.where`,
     `jnp.take`, or explicit index arrays.
   - Coarse `.at[..., slice].set(...)` updates are acceptable when they build
     Fortran-like banded matrices, interleaved solve vectors, boundary values, or
     stats arrays. The smell is scalar/vertical repeated updates, not every
     `.at.set`.
   - Use `jax.lax.scan` for true recurrences where each level depends on the
     previous level's computed state. Do not hide a recurrence inside a Python
     loop and do not vectorize it incorrectly just to remove the loop.

3. Keep Python control flow static or replace it with JAX control flow.
   - Python `if` statements are fine for compile-time/static configuration such
     as solver choice, enabled physics options, or fixed dimensions.
   - Do not branch in Python on JAX array values: avoid patterns such as
     `if bool(jnp.any(mask))`, `if jnp.any(mask)`, `if array_value`, or
     `if np.any(jnp_array)`.
   - Use array masks, `jnp.where`, `jax.lax.cond`, `jax.lax.select`, or carry an
     active/fatal mask through the computation when the branch condition depends
     on model data.
   - Avoid data-dependent `break`, `continue`, and early `return` inside code
     intended to become a pure JAX kernel. Carry status through the computation
     or use structured JAX control flow.

4. Do not extract Python scalars from traced values.
   - Avoid `float(array_value)`, `int(array_value)`, `.item()`, and
     `bool(array_value)` when the value may be a JAX array or tracer.
   - If a value is a static configuration quantity, keep it outside the pure
     kernel or pass it as a static argument.
   - If a value is model data, keep it as a JAX scalar and use it directly in
     JAX expressions.

5. Keep shapes and object state JIT-friendly.
   - Dimensions used in `jnp.zeros`, `reshape`, `stack`, and solver allocation
     must be static for JIT. Treat `nzm`, `nzt`, `ngrdcol`, scalar counts, and
     similar shape controls as static kernel metadata unless the design says
     otherwise.
   - Avoid dynamic Python lists, `.append`, and variable-length tuple assembly
     in pure kernels. Prefer fixed tuples, stacked arrays, or static loops over
     known field/scalar counts.
   - Python objects such as grid, error-info, PDF-parameter, and vertical-grid
     derived types should not be mutated inside pure JAX kernels. Either treat
     metadata fields as static or pass array fields explicitly. If pytrees are
     introduced, make the static versus dynamic fields explicit.

6. Keep API and stats details from shaping the pure JAX design.
   - Do not add JAX-only optional arguments, local stats flags, or derived-type
     setter calls to make API calls work.
   - If Fortran global state such as `stats%l_sample` is needed, expose it
     through the shared API, for example a small public helper, instead of adding
     a target-only optional argument.
   - Keep stats collection outside pure compiled physics where practical, or
     represent diagnostics as explicit returned arrays.
   - Do not let API boundary needs justify NumPy math throughout the port.
     Physics code should remain JAX-first, with conversion isolated to explicit
     boundaries.

7. Prefer JAX kernel style over eager Python mirror style.
   - Prefer vectorized internal kernels, static-shape arrays, whole-slice
     updates, masks for conditional updates, and `lax.scan` for recurrences.
   - Avoid eager Python/NumPy mirror code, mutable object state,
     data-dependent Python branching, local compatibility adapters, and
     one-level-at-a-time `.at` updates.

Recommended porting workflow:

1. Routine inventory pass.
   - Ensure no Fortran routine is missing and no unjustified target-only helper
     was added.

2. Signature pass.
   - Make function/subroutine arguments match names, order, grouping, and
     optional status.
   - Update call sites after signature changes.

3. Block-order pass.
   - Walk the Fortran file from top to bottom and compare each block with the
     target file.
   - Move target blocks until the order matches.

4. Comment pass.
   - Restore missing Fortran comments and divider blocks.
   - Remove unnecessary target-only commentary.

5. Call-format pass.
   - Convert target calls to positional arguments where possible.
   - Match the Fortran line grouping at each call site.
   - Mark only unavoidable keyword calls with a TODO.

6. Calculation pass.
   - Compare formulas, temporaries, initialization, array dimensions, clipping,
     and conditional updates.
   - Remove renamed aliases and unnecessary shape helpers.

7. Error/diagnostic pass.
   - Restore Fortran debug gates, parameterization checks, stats calls, and fatal
     return behavior.

8. Wrapper-removal pass.
   - Search the target file for old top-level wrapper calls such as
     `f2py_advance_<routine>` or equivalent API shortcuts.
   - Keep lower-level calls only when they are deliberate adaptation boundaries.
   - Confirm imports do not include unused compatibility layers or fallbacks that
     hide a broken target runtime.
   - Search for raw generated bindings such as `clubb_f2py`; target physics
     ports should not import them. Fix the public API wrapper instead.
   - Search for local helpers that only wrap API calls or dtype/layout
     conversion. Delete them or move genuinely needed behavior to the shared API.

9. Branch-coverage audit pass.
   - Identify branch-specific behavior that could be missed by a straight-line
     read.
   - Pay special attention to passive scalar, wind-prediction, perturbed-wind,
     PDF-type, implicit versus explicit turbulent-advection, clipping, stats,
     and diagnostic branches when they exist in the Fortran.
   - Do not assume that matching the common path proves the less common branches
     were ported correctly.

10. Compliance audit pass.
   - Compare the final target routine order against the Fortran outline.
   - Review every target-only helper and classify it as removed, unavoidable, or
     still needing cleanup.
   - For every retained target-only helper or unresolved conflict, add a local
     in-code TODO at the relevant import, helper, branch, call site, or
     conversion block. The TODO should briefly say why it remains and what would
     remove it.
   - Search explicitly for helper and adapter smells such as `def _`, `_f64`,
     `_np64`, `_stats`, `_clip`, `_fatal`, `compat`, `adapter`, and `wrapper`.
     Classify each occurrence as mirrored from Fortran, unavoidable and
     documented, useful decomposition with clear reasoning, or a defect to
     remove.
   - Ensure every unsupported or omitted Fortran block has a nearby in-code
     TODO with its reason.
   - Check diagnostic/error routines and signature drift explicitly.
   - Check API-boundary drift explicitly: no raw F2PY imports in target ports, no
     one-line local wrappers around public API calls, and no defensive
     multi-contract adapters for API return values.
   - Re-read non-default branches even if the comparison cases are unlikely to
     exercise them. Include debug fatal paths, per-column versus all-column error
     updates, optional solver paths, scalar loops, perturbed-wind paths,
     single-LHS versus multiple-LHS paths, stats-enabled paths, and clipping
     paths when they exist.
   - Do this before claiming the port is complete; do not rely on tests to catch
     inactive branches or structural drift.

11. Validation expectation.
   - Testing is mandatory, but detailed test selection belongs in the user's
     current request, repo-local test guidance, or task-specific debugging plan.
   - Passing tests are necessary but not sufficient. After tests pass, still do
     the structural and helper hygiene audit; tests can miss inactive branches
     and target-language convenience artifacts.
   - After porting, validate the changed target code enough to support the
     behavioral claim being made: syntax/import health, focused behavior, and
     Fortran-vs-target equivalence are different levels of evidence.
   - When porting code that affects a Python or JAX driver path, include
     whole-result comparisons where possible, not only unit tests. In this repo,
     examples include `run_python_vs_fortran_cases.py` and
     `run_jax_vs_fortran_cases.py`.
   - For broad timestep or shared-driver changes, focused validation alone is
     not enough; broader comparison validation should follow once focused issues
     are fixed.
   - Do not add target-runtime fallbacks that hide missing dependencies or broken
     execution. A failed target runtime should fail clearly.
   - In the final response, state what validation was run, or explicitly state
     why validation was not run.

Final review checklist:

1. The target file has the same routine inventory as the Fortran portion.
2. Routines appear in the same order as the Fortran unless a documented
   target-language constraint prevents that.
3. No target-only helper exists unless it is documented as unavoidable.
4. Argument names, order, grouping, optional status, and returns match the Fortran
   or documented target API surface.
5. Calls appear in the same order as the Fortran.
6. Call arguments are positional and grouped like the Fortran where practical.
7. Comments and divider blocks match the Fortran.
8. Calculations, temporaries, array dimensions, and conditionals are mechanically
   equivalent to the Fortran.
9. Debug checks, validation routines, stats updates, and fatal returns are not
   missing.
10. Diagnostic stubs, unsupported features, and signature drift are documented
    as incomplete or unavoidable; none are silent.
11. No alias-only variables, shape tuple shortcuts, or compatibility kwargs were
   introduced.
12. Appropriate validation was run, or any missing validation is explained.
13. No top-level Fortran/F2PY wrapper remains in place of the requested port.
14. No raw generated Fortran/F2PY binding is imported or called from a target
    physics port; any needed low-level operation goes through the public API.
15. No one-line local helper merely wraps `clubb_api`, dtype conversion, array
    layout conversion, or blank-name filtering.
16. API calls rely on one stable return contract; no local defensive adapter
    accepts multiple shapes for the same API result.
17. Field-specific call arguments were compared against the Fortran call sites,
    not inferred from similar fields.
18. Broad shared behavior was validated beyond only the easiest/common path when
    the port touched shared driver behavior.
19. Every retained target-only helper, reduced diagnostic, unsupported block,
    reordered routine, or delegated kernel has a nearby in-code TODO with a
    short reason and expected removal path.
20. Passing tests were followed by a structural mirror audit; tests alone were
    not treated as port completion.
