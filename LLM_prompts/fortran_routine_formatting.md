Apply CLUBB Fortran routine formatting to new or changed subroutines/functions.

Use this prompt when adding, moving, extracting, or materially changing Fortran routines in `src/CLUBB_core` or nearby CLUBB Fortran code. The goal is to keep routine interfaces readable, reviewable, and consistent with the repo's local style.

Before editing:

1. Inspect the surrounding module and at least two nearby routines.
   - Match the local indentation, comment divider style, `use` ordering, and declaration grouping where it is already consistent.
   - Prefer the existing module's vocabulary for comments such as `thermodynamic levels`, `momentum levels`, `CLUBB tunable parameter`, `PDF parameters`, and similar units/descriptions.
   - If the nearby code conflicts with this prompt, preserve behavior first, then improve only the routine(s) you are touching.

2. Identify every changed routine signature and every call site.
   - Include direct call sites, wrapper call sites, tests, Python/f2py wrapper signatures, and generated/ported driver mirrors if applicable.
   - Keep call-site argument grouping aligned with the routine definition.

Routine header and description:

1. Every new or changed routine must have a `Description:` section immediately below the subroutine/function statement and before `use` statements.
2. Include `References:` when the surrounding module normally includes it, or when the routine implements a named paper, scheme, or published formula.
3. Keep the description factual and short. Explain what the routine computes and any important ordering/cache reason for the routine existing.

Argument order:

1. Order routine arguments strictly by intent group:
   - `intent(in)`
   - `intent(inout)`
   - `intent(out)`
   - `intent(in), optional`
   - `intent(inout), optional`
   - `intent(out), optional`
2. Do not mix intent groups on the same argument-list continuation line.
3. Within each intent group, keep related variables together:
   - dimensions and grid/control flags first (`nzm`, `nzt`, `ngrdcol`, `gr`, logical flags)
   - physical state/input fields next
   - tunable parameters/configuration inputs after fields
   - `stats` and `err_info` with the appropriate `InOut` group
   - output fields grouped by grid/location or physical role
4. Preserve BFB-sensitive ordering unless the user explicitly accepts behavior changes.

Call-site formatting:

1. Format call sites to mirror the routine definition grouping.
2. Each continuation line at the call site must contain arguments from only one intent group.
3. Add an end-of-line comment for each call-site continuation line using only:
   - `! In`
   - `! InOut`
   - `! Out`
   - `! Optional In`
   - `! Optional InOut`
   - `! Optional Out`
4. Do not use longer comments such as `! Intent(in)` at call sites for new formatting.
5. If preprocessor conditionals split a call argument list, keep the conditional argument in the correct intent group and comment it consistently.

Declaration formatting inside routines:

1. Declare variables in the same order they appear in the argument list, grouped by intent.
2. Use explicit section dividers before each declaration group. Prefer this shape unless the surrounding module has a stronger local convention:

   ```fortran
    !--------------------------- Input Variables ---------------------------
   ```

   Use analogous labels:
   - `Input Variables`
   - `Input/Output Variables`
   - `Output Variables`
   - `Optional Input Variables`
   - `Optional Input/Output Variables`
   - `Optional Output Variables`
   - `Local Variables`
3. End the declaration section with:

   ```fortran
    !--------------------------- Begin Code ---------------------------
   ```

4. Every declared variable should have a short description, preferably with units in square brackets.
   - Put comments inline for compact declarations.
   - Use aligned trailing comments for multi-line declarations where it improves readability.
   - Use `[-]` for dimensionless values.
   - Use `[units vary]` only when a variable genuinely carries case/scalar-dependent units.
5. Do not leave undocumented temporary arrays. If a temporary exists only to avoid OpenACC/function side effects or preserve BFB behavior, say that briefly.
6. Keep intrinsic declarations, parameters, and local constants near the relevant declaration group and document them too.

Style details observed in this repo:

1. Use `real( kind = core_rknd )`, not shortened kind syntax, unless the surrounding routine already consistently uses another style.
2. Use `type(grid)` or `type (grid)` consistently with nearby code; do not churn existing style outside the touched routine.
3. Keep `implicit none` after `use` statements and before declarations.
4. Use `! In`, `! InOut`, and `! Out` in call sites, with capitalization exactly as shown.
5. Prefer local helper routines for diagnostic/cache calculations only when the helper owns the consistency point and all dependent outputs are returned or updated there.
6. For OpenACC data regions, keep local scratch arrays created/deleted inside the routine that owns them. Do not leave scratch arrays in a parent routine after extracting the calculation.
7. If a routine moves a stats update inside itself, remove duplicate stats updates from the caller and keep the stats update at the same logical point as the calculation.
8. When a diagnostic value must match mutable prognostic values, either compute it after the mutation or make the routine that mutates the prognostic values also update the dependent diagnostic. Do not introduce stale diagnostic caches.

Review checklist before finishing:

1. Routine signature grouping and declaration grouping match.
2. Every call site mirrors the new argument grouping and has intent comments.
3. Every argument and local variable has a useful description and units where possible.
4. `stats` and `err_info` are in the correct intent group at both definition and call sites.
5. Optional arguments are declared and passed after required arguments.
6. OpenACC create/delete lists match moved local variables.
7. No removed argument remains in wrappers, tests, Python API files, or driver mirrors.
8. Run at least:

   ```bash
   git diff --check
   ```

9. If the build/test environment is available, compile the touched Fortran path or run the smallest relevant test. If not available, state the environment blocker clearly.
