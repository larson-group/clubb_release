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
- Do not introduce abstractions, helpers, optionals, aliases, or reordered logic
  just because they are convenient in the target language.

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

3. Decide what counts as a necessary language adaptation.
   - Index-base conversion is necessary, for example Fortran 1-based parameter
     indices mapped to Python 0-based array indices.
   - Array-layout conversion may be necessary.
   - Different syntax for inout/out returns may be necessary.
   - A target-language API wrapper may require a small adapter.
   - Document any adaptation that prevents a literal mirror.

Routine inventory:

1. The target file should have a 1-to-1 routine inventory with the Fortran file
   for the portion being ported.
   - Do not omit Fortran routines because they are private helpers.
   - Do not create new helpers that do not exist in the Fortran.
   - If the Fortran has a helper routine, port that helper using the same name and
     relative position where possible.
   - If a helper truly cannot be ported because it belongs in a lower-level API,
     state that explicitly and keep the target call site as close as possible.

2. Preserve routine names.
   - Use the Fortran routine name unless the target language has a hard naming
     restriction.
   - Avoid "compat", "adapter", "impl", or underscored aliases unless they are
     already established in the Fortran or required by the target binding layer.
   - Do not use temporary renamed variables such as `_name` only to unpack and
     then assign back to `name`. Unpack or assign directly to the Fortran name.

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

Comments and block structure:

1. Preserve Fortran comments and divider blocks.
   - Keep comments in the same relative location.
   - Keep major block headers even if the target language code below them looks
     different syntactically.
   - Do not delete diagnostic, error-checking, or setup comments because the
     corresponding target code is shorter.

2. Remove target-only explanatory comments unless they document a real language
   adaptation.
   - Do not narrate obvious assignments.
   - Do not add comments that explain the porting process inside the code.
   - If a comment exists only because a target-language workaround is needed,
     keep it short and specific.

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

Recommended iterative workflow:

1. Routine inventory pass.
   - Ensure no Fortran routine is missing and no target-only helper was added.

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

8. Validation pass.
   - Run syntax checks for the target language.
   - Run focused tests for touched routines.
   - Run available argument-contract audits if wrappers are involved.
   - Run Fortran-vs-port comparison tests for representative cases.

Validation guidance for CLUBB Python/JAX ports:

1. For Python files, run at least:

   ```bash
   python -m py_compile <changed-python-files>
   ```

2. If the Python API was touched, run:

   ```bash
   python -m pytest clubb_python_api/tests
   ```

3. For Python driver behavior, start with short Fortran-vs-port comparisons
   before the full suite:

   ```bash
   ./run_scripts/run_python_vs_fortran_cases.py --cases bomex --jobs 1
   ./run_scripts/run_python_vs_fortran_cases.py --cases atex --jobs 1
   ```

4. If a JAX mirror exists and was changed, run the analogous JAX checks.

5. Always run:

   ```bash
   git diff --check
   ```

Final review checklist:

1. The target file has the same routine inventory as the Fortran portion.
2. No target-only helper exists unless it is documented as unavoidable.
3. Argument names, order, grouping, optional status, and returns match the Fortran
   or documented target API surface.
4. Calls appear in the same order as the Fortran.
5. Call arguments are positional and grouped like the Fortran where practical.
6. Comments and divider blocks match the Fortran.
7. Calculations, temporaries, array dimensions, and conditionals are mechanically
   equivalent to the Fortran.
8. Debug checks, validation routines, stats updates, and fatal returns are not
   missing.
9. No alias-only variables, shape tuple shortcuts, or compatibility kwargs were
   introduced.
10. Tests and comparisons were run, or any missing validation is explained.
