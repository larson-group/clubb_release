Update the Python API and Python/JAX drivers so they match the current Fortran CLUBB core.

The broad goal is to make the Python API, Python driver, and JAX driver structurally match the current Fortran implementation after Fortran refactors. This usually means updating f2py-facing wrapper signatures, fixing changed cross-module public routine argument lists, moving Python-side calculations to match newly internalized Fortran helper behavior, and keeping the Python and JAX driver copies in sync.

Do not automatically apply this prompt after every Fortran change that might affect Python. If you are doing other Fortran work and notice that argument lists for cross-module public routines have changed, remind the user that the Python API and Python/JAX drivers may need follow-up updates, then ask whether they want that fix-up done now. Apply this prompt only when the user has asked for Python API/driver work, has agreed to the follow-up, or the current task explicitly includes keeping Python wrappers and drivers in sync.

Do the work in this order:

1. Inspect the recent Fortran changes in `src/CLUBB_core`.
   - Focus on public routines whose argument lists changed.
   - Look for routines that moved between modules.
   - Look for helper routines that now compute values internally that Python may still be recomputing.
   - Treat any changed public Fortran routine as a possible Python API impact, including changed arguments, changed intent, added optional arguments, removed outputs, renamed routines, moved routines, deleted routines, and newly added routines.
   - If a whole new Fortran source file or module was added and it exposes functionality needed by the Python driver, expect that a new Python API wrapper or Python-side port may be needed.
   - Do not limit the scan to routines used by `advance_clubb_core`; cross-module public routines anywhere in `src/CLUBB_core` can affect the Python API or drivers if they are wrapped, imported, or mirrored in Python.
   - Also read the relevant Python API and driver README files so the current API/driver split is clear before editing. Check for files such as:

     ```text
     clubb_python_api/README*
     clubb_python_driver/README*
     clubb_jax/README*
     ```

2. Update the Python API wrappers first.
   - Fix changed argument lists in `clubb_python_api/`.
   - Add wrappers for new public routines if needed.
   - Remove or update wrappers for routines that were deleted, renamed, moved, or had outputs made local.
   - Keep wrapper argument order aligned with the Fortran routine order.
   - Treat the current Fortran public routine definition as the authority for the top-level Python API signature. For every wrapped public routine, the exported `clubb_api.<routine>` signature and the F2PY-facing wrapper must match the current Fortran argument set, order, and ownership unless a documented adapter exception is truly required.
   - Do not expose stale adapter arguments in `clubb_api` merely because older wrappers or drivers still pass them. If Fortran now computes a value internally, such as a value that used to be supplied by Python, remove that value from the public Python API and from driver call sites unless another current Fortran argument still requires it.
   - Be especially suspicious of wrapper-only arguments that are not passed through to the underlying Fortran routine. These are usually signature drift and should be removed or explicitly documented as adapter-only reconstruction inputs.
   - If a routine gains a Fortran derived-type argument or return value, add or update the corresponding Python-to-Fortran type conversion path. Check the existing `derived_types` converters and `derived_type_storage` patterns before creating anything new.

3. Compile the Python API.
   - First make sure the compiler and NetCDF Fortran environment are available.
   - If the current shell already has `module` or `ml`, load the needed modules directly:

     ```bash
     module load gcc netcdf-fortran
     ./compile.py -python
     ```

   - If `module` is not available, try the site-specific environment setup for the machine you are on. Common options include:

     ```bash
     source /etc/profile.d/larson-group.sh
     source /etc/profile
     source /usr/share/lmod/lmod/init/bash
     source /usr/local/lmod/lmod/init/bash
     source /usr/local/spack/share/spack/setup-env.sh
     ```

   - After sourcing an environment setup file, check whether `module` or `ml` is available, then load the closest available compiler and NetCDF Fortran modules. On some systems the module names may include versions, for example `gcc/13.1.0` instead of `gcc`.
   - If modules are not used on the current machine, inspect the repo's existing build documentation or compiler config files and use the local equivalent environment.
   - Fix compile errors before moving on.

4. Run the Python API tests.
   - Use:

     ```bash
     python -m pytest clubb_python_api/tests
     ```

   - Some negative-path tests may print expected Fortran fatal-error text. Trust the pytest result, not just the log noise.

5. Update the Python and JAX drivers.
   - Update `clubb_python_driver/advance_clubb_core.py`.
   - Mirror equivalent changes in `clubb_jax/advance_clubb_core.py`.
   - Use both driver paths in `advance_clubb_to_end.py` as validation tools:

     ```text
     clubb_python_driver/advance_clubb_to_end.py
     clubb_jax/advance_clubb_to_end.py
     ```

   - `_advance_clubb_core_api` exercises the Python API binding to the Fortran `advance_clubb_core` path. This is useful for checking API functionality and wrapper compatibility.
   - `_advance_clubb_core_python` exercises the Python port of the core logic and checks much more of the Python-side implementation. This path is usually more important for validating that the Python/JAX drivers are up to date with the Fortran refactor.
   - When debugging, swap between `_advance_clubb_core_api` and `_advance_clubb_core_python` intentionally so both paths are known-good. Keep track of which path each comparison is testing.

6. Run small SCM smoke comparisons before running the full suites.
   - Start with `bomex` or `atex` and a small iteration count:

     ```bash
     python run_scripts/run_scm.py -max_iters 5 -out_dir /tmp/clubb_bomex_fortran bomex
     python run_scripts/run_scm.py -max_iters 5 -python -out_dir /tmp/clubb_bomex_python bomex
     python run_scripts/run_bindiff_all.py -v 2 -case bomex -t 1e-12 -pt 1e-12 /tmp/clubb_bomex_python /tmp/clubb_bomex_fortran
     ```

   - Keep the iteration count low while debugging so the first diverging field and timestep are easier to isolate.

7. Run targeted comparison harness checks.
   - Use serial jobs first:

     ```bash
     ./tests/run_python_vs_fortran_cases.py --cases bomex --jobs 1
     ./tests/run_python_vs_fortran_cases.py --cases atex --jobs 1
     ./tests/run_jax_vs_fortran_cases.py --cases bomex --jobs 1
     ./tests/run_jax_vs_fortran_cases.py --cases atex --jobs 1
     ```

   - `--jobs 1` makes logs easier to inspect and avoids multiprocessing-related noise.

8. Final success criterion: run the full comparison suites.
   - The work is complete only when these pass, or any remaining differences are understood, documented, and explicitly accepted:

     ```bash
     ./tests/run_python_vs_fortran_cases.py --jobs 1
     ./tests/run_jax_vs_fortran_cases.py --jobs 1
     ```

Debugging techniques and common failure modes:

- If many prognostic fields diverge, look for duplicated calculations in Python that Fortran routines now do internally. Common examples are interpolation calls like `zm2zt`, `zt2zm`, or `zm2zt2zm` that may have moved inside a Fortran helper.

- If only stats fields differ, look for stale explicit Python-side `stats_update` calls. These may now belong inside Fortran helper wrappers or a centralized stats routine.

- Check logs for stats oversampling warnings. They often identify duplicate stats updates directly.

- Watch for BFB-sensitive ordering. If Fortran moved a computation inside a helper, the Python driver should generally call the same helper and avoid recreating intermediate arrays in a slightly different order.

- Confirm index-base expectations when calling Fortran wrappers from Python. Python grid bounds may be zero-based, while Fortran routines or wrappers may expect one-based bounds.

- If a Fortran output was made local inside a routine, remove it from Python call sites unless another downstream Python routine still needs the same value. If it is still needed, compute it at the same point in the workflow as Fortran does.

- Passing tests are not enough to prove API correctness. Existing wrapper/order audits may only compare Python wrappers against the generated F2PY interface; they do not necessarily prove that either layer still matches the current Fortran source routine. When wrapper drift is suspected, compare directly against the Fortran subroutine definition.

- If a Fortran stats call moved inside a helper, make sure the Python path either calls the helper that performs the stats update or performs an equivalent update at the same logical time.

- If the Python driver and JAX driver share the same structure, fix both together. Do not leave the JAX driver stale after fixing the Python driver.

- Use NetCDF field comparisons with tight thresholds after each small fix. This keeps the first changed field and earliest changed timestep visible.

- Treat `bomex` and `atex` as smoke tests only. They are useful early checks, but they are not the final target.

Expected final state:

- `./compile.py -python` succeeds.
- `python -m pytest clubb_python_api/tests` passes.
- Short `run_scm.py` Python-vs-Fortran bindiffs pass for representative cases.
- `./tests/run_python_vs_fortran_cases.py --jobs 1` passes.
- `./tests/run_jax_vs_fortran_cases.py --jobs 1` passes.
