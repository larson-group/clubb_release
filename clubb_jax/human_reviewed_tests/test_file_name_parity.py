"""Check that the CLUBB-JAX source tree mirrors the CLUBB Fortran tree file for file.

The JAX port is meant to stay easy to navigate for someone who knows the Fortran code: the JAX version of
src/<path>/<name>.F90 should be clubb_jax/src/<path>/<name>.py. Fortran forbids circular `use` between
modules, so a file-for-file mirror cannot create a circular import in Python. There is no technical reason to
rename, merge or split files during the port.

This module compares the two trees and fails on three kinds of structural drift:

  * rename -- a Fortran <name>.F90 ported to a JAX file with a different name
  * fold   -- a Fortran file's code merged into some other JAX file
  * split  -- part of one Fortran file moved into a new JAX file of its own

Usage
-----
Run from the repository root (the clubb_jax Jenkins job's "File-name parity" stage runs the same commands):

    python3 clubb_jax/run_jax.py -profile=cpu -init_env
    .venv-jax/bin/python -m pytest -v clubb_jax/human_reviewed_tests/test_file_name_parity.py

The check uses only the Python standard library and reads file names, so it needs no build and runs in well
under a second. The JAX environment is used only because it provides pytest.

Understanding the lists
-----------------------
Every difference between the trees must be fixed or listed in one of the tables below with its reason.
Renames, folds and splits are never listed: the failure output names each one, so it doubles as a work list,
and the test keeps failing until they are fixed by renaming, splitting or merging the JAX files.

  * FORTRAN_ONLY_ALLOWED_DIRS -- whole Fortran subtrees outside the port; their files are not compared.
  * FORTRAN_ONLY_ALLOWED -- Fortran files whose code is genuinely NOT in JAX (out-of-scope subsystems,
    features the JAX driver rejects at initialization, Fortran-only build and host machinery).
  * JAX_ONLY_ALLOWED -- JAX files with no Fortran source to mirror at all.

test_lists_are_still_live keeps the tables honest: an entry that no longer describes a real difference (the file
was ported, renamed or deleted) fails, so an exception cannot outlive its reason.
"""
import os

import pytest

from utilities.output_paths import REPO_ROOT

_FORT_ROOT = os.path.join(REPO_ROOT, "src")
_JAX_ROOT = os.path.join(REPO_ROOT, "clubb_jax", "src")

# Every Fortran source form in this tree. Fixed-form .F/.f and the include files are real ported source too
# (e.g. BUGSrad/gases_ckd_data.h has a JAX twin), so they are all compared.
_FORTRAN_EXTS = (".F90", ".f90", ".F", ".f", ".h", ".inc")

# Fortran subtrees (path relative to src/) with no clubb_jax/src/ counterpart at all -- whole subsystems
# outside the supported port. Files under these are not compared. If a JAX counterpart directory ever
# appears, test_lists_are_still_live fails so the exclusion gets revisited.
FORTRAN_ONLY_ALLOWED_DIRS = {
    "Lapack": "Vendored LAPACK source (85 .f files). JAX uses jnp.linalg / jax.scipy.linalg.",
    "G_unit_test_types": "Fortran's own unit-test framework; used only by src/G_unit_tests.F90. JAX has its "
                         "own pytests under clubb_jax/auto_llm_generated_tests/pytests/.",
    "Numerical_recipes": "Vendored Numerical Recipes routines (simulated annealing, quicksort, RNG state). "
                         "Used only by the Fortran tuner (clubb_tuner.F90, error.F90, tuner_tests.F90).",
    "readiopdata_core": "IOP (SCAM netCDF) reader core; used only by src/readiopdata_module.F90, itself not "
                        "ported -- see that entry.",
    "Microphys/COAMPS_microphys": "COAMPS microphysics scheme (fixed-form .F). JAX keeps only disabled "
                                  "interfaces in Microphys/coamps_microphys_driver_module.py.",
    "Microphys/SCM_Activation": "GFDL aerosol-activation code. Initialization rejects this feature; "
                                "Microphys/gfdl_activation.py holds only disabled signature dummies.",
}

# Fortran files (path relative to src/, without extension) whose code is genuinely NOT ported.
# NOT for renames/folds/splits -- see this file's docstring.
FORTRAN_ONLY_ALLOWED_FILES = {
    # ---- src/ root: Fortran test programs, the legacy tuner and host machinery ----
    "clubb_loss_driver_test": "Fortran test program for the loss driver. JAX loss coverage is "
                              "clubb_jax/auto_llm_generated_tests/integration_tests/test_loss_driver.py.",
    "clubb_thread_test":"Fortran test program (thread-safety); JAX has no OpenMP equivalent to test.",
    "G_unit_tests": "Fortran unit-test driver program; JAX's equivalent is the pytest suite.",
    "generalized_grid_test": "Grid-generalization self-test invoked from clubb_driver.F90. The JAX driver "
                             "rejects it (clubb_case_initalization.py: l_test_grid_generalization).",
    "clubb_tuner": "Legacy Fortran tuner program. JAX takes part in tuning through the Python tuner "
                   "(tuner/, loss through clubb_loss_driver.py), not through a port of this program.",
    "error": "Despite the name, this module is the legacy Fortran tuner (tuner_init, min_les_clubb_diff, "
             "output_nml_tuner), not error handling. Not ported, as with clubb_tuner.",
    "readiopdata_module": "IOP (SCAM netCDF) sounding/surface reader, used when l_readiopdata is true. That "
                          "input path is not ported.",
    "text_writer": "write_text/write_date formatting helpers for Fortran's text log output; host I/O.",

    # ---- CLUBB_core ----
    "CLUBB_core/lapack_wrap": "LAPACK band/tridiagonal solver wrappers -- JAX uses jnp.linalg.solve "
                              "(matrix_solver_wrapper.py) instead of calling LAPACK.",
    "CLUBB_core/lapack_interfaces": "LAPACK single/double precision interface declarations; same as "
                                    "lapack_wrap.",
    "CLUBB_core/clubb_api_module": "Doorway for Fortran host models (E3SM, CAM, WRF): re-exports plus *_api "
                                   "wrappers that mostly copy single-column arguments into one-column arrays. "
                                   "JAX has no host models; its driver calls the underlying routines directly, "
                                   "and *_api routines owned by other modules are ported in those modules. "
                                   "Revisit if clubb_driver.F90, which uses this module, is mirrored.",
    "CLUBB_core/endian": "Byte-order detection for Fortran binary I/O. No callers anywhere in src/.",
    "CLUBB_core/grid_adaptation_module": "Adaptive grid remapping. The JAX driver rejects it "
                                         "(clubb_case_initalization.py: grid_adapt_in_time_method > 0).",
    "CLUBB_core/penta_bicgstab_solver": "Experimental iterative pentadiagonal solver; its own header says it "
                                        "is slower than the direct LU solver even on GPU, and is not "
                                        "GPU-ized. Only penta_lu_solver is ported.",

    # ---- Benchmark_cases ----
    "Benchmark_cases/clex9_nov02": "Deliberately excluded per prescribe_forcings.py's docstring: the JAX "
                                   "driver 'intentionally omits the Fortran CLEX November case ... Its "
                                   "case-specific input reader has not been ported.'",

    # ---- Input_fields: LES input fields ----
    "Input_fields/extrapolation": "lin_ext_zm_bottom/lin_ext_zt_bottom, used by input_fields.F90 only in the "
                                  "LES input-fields branch, which the JAX port of input_fields.F90 leaves "
                                  "gated (the JAX driver rejects l_input_fields).",

    # ---- Microphys ----
    "Microphys/microphys_interface": "Include file declaring the Fortran `interface` of the microphysics "
                                     "subroutine argument passed to estimate_scm_microphys_module and "
                                     "lh_microphys_driver_module. Python passes the function directly and "
                                     "needs no interface declaration.",
    "Microphys/KK_microphys/AiryFunction": "Airy-function routines that the KK parabolic cylinder orders do "
                                           "not need. Parabolic.f90 uses only xpowy (x**y), which "
                                           "parabolic_cylinder.py writes as `**`.",
    "Microphys/Morrison_microphys/microphysics": "M2005 host driver, commented out of the Fortran build "
                                                 "itself (Morrison_microphys/CMakeLists.txt: 'current "
                                                 "unused').",

    # ---- Radiation/BUGSrad: variants the CLUBB build compiles out ----
    "Radiation/BUGSrad/kinds": "Fortran KIND parameters (int_kind, dbl_kind, ...). No Python equivalent; JAX "
                               "arrays carry their dtype.",
    "Radiation/BUGSrad/driver_read": "Standalone BUGSrad test program; nothing in CLUBB calls it.",
    "Radiation/BUGSrad/two_rt_lw_iter": "Cloud-overlap solver variant. CLUBB compiles with -Dnooverlap "
                                        "(cmake/utils/CPPDefinitions.cmake), so only two_rt_lw is built.",
    "Radiation/BUGSrad/two_rt_lw_sel": "Overlap solver variant compiled out by -Dnooverlap.",
    "Radiation/BUGSrad/two_rt_lw_ocastrndm": "Overlap solver variant, additionally behind "
                                             "USE_BUGSrad_ocast_random, which CLUBB does not define.",
    "Radiation/BUGSrad/two_rt_sw_bs": "Overlap solver variant compiled out by -Dnooverlap.",
    "Radiation/BUGSrad/two_rt_sw_iter": "Overlap solver variant; no callers at all.",
    "Radiation/BUGSrad/two_rt_sw_ocastrndm": "Overlap solver variant, additionally behind "
                                             "USE_BUGSrad_ocast_random, which CLUBB does not define.",

    # ---- SILHS ----
    "SILHS/mt95": "Mersenne Twister random number generator. SILHS/generate_uniform_sample_module.py draws "
                  "with jax.random keys seeded by Fortran's lh_seed * itime convention, so random streams "
                  "differ from Fortran on purpose. l_lh_deterministic_test gives both codes the same "
                  "hard-coded sampling inputs for comparisons (see clubb_jax/SILHS_PORT_NOTES.md).",
}

# JAX files (path relative to clubb_jax/src/, without extension) with no Fortran source to mirror.
# NOT for renames/folds/splits -- see this file's docstring.
JAX_ONLY_ALLOWED_FILES = {
}

def _twin(path, old_exts, new_root, new_ext):
    """Return the repo-relative path that the mirrored counterpart of `path` would have.

    `path` is repo-relative (e.g. src/CLUBB_core/foo.F90); its first directory is replaced by `new_root` and
    its extension, one of `old_exts`, by `new_ext`. The directory casing of `path` is kept.
    """
    rel = path.split(os.sep, 1)[1] if os.sep in path else path
    for ext in old_exts:
        if rel.endswith(ext):
            rel = rel[: -len(ext)]
            break
    return f"{new_root}/{rel}{new_ext}"


def _keys(table):
    """Return the lowercased keys of a table; table keys keep the real casing, comparisons ignore it."""
    return {key.lower() for key in table}


def _rel_key(path, root, exts):
    """Return the comparison key for `path`: relative to `root`, extension removed, '/'-separated, lowercase."""
    rel = os.path.relpath(path, root)
    for ext in exts:
        if rel.endswith(ext):
            rel = rel[: -len(ext)]
            break
    return rel.replace(os.sep, "/").lower()


def _is_excluded(rel_dir):
    """Return whether a directory relative to src/ lies in one of FORTRAN_ONLY_ALLOWED_DIRS."""
    rel_dir = rel_dir.replace(os.sep, "/")
    return any(rel_dir == d or rel_dir.startswith(d + "/") for d in FORTRAN_ONLY_ALLOWED_DIRS)


def _fortran_files():
    """Return {comparison key: repo-relative path} for every Fortran source outside the excluded subtrees."""
    found = {}
    for dirpath, dirnames, filenames in os.walk(_FORT_ROOT):
        rel_dir = os.path.relpath(dirpath, _FORT_ROOT)
        rel_dir = "" if rel_dir == "." else rel_dir
        if rel_dir and _is_excluded(rel_dir):
            dirnames[:] = []
            continue
        for name in filenames:
            if name.endswith(_FORTRAN_EXTS):
                full = os.path.join(dirpath, name)
                found[_rel_key(full, _FORT_ROOT, _FORTRAN_EXTS)] = os.path.relpath(full, REPO_ROOT)
    return found


def _jax_files():
    """Return {comparison key: repo-relative path} for every JAX source module (package __init__ excluded)."""
    found = {}
    for dirpath, dirnames, filenames in os.walk(_JAX_ROOT):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for name in filenames:
            if name.endswith(".py") and not name.startswith("__"):
                full = os.path.join(dirpath, name)
                found[_rel_key(full, _JAX_ROOT, (".py",))] = os.path.relpath(full, REPO_ROOT)
    return found


def test_no_unexplained_fortran_only_files():
    """Every compared Fortran source has a same-named JAX file or an entry explaining why not.

    Fails with the list of Fortran files that have neither, each with the JAX path it should have.
    """
    fort, jax = _fortran_files(), _jax_files()
    listed = _keys(FORTRAN_ONLY_ALLOWED_FILES)
    missing = sorted(k for k in fort.keys() - jax.keys() if k not in listed)
    if missing:
        # pytrace=False reports only this message: the file list, without pytest's source/traceback dump.
        pytest.fail(
            "Fortran source with no same-named JAX file (port it, or -- if its code IS in JAX under another "
            "name/file -- rename or split the JAX file so the names line up; only add a FORTRAN_ONLY_ALLOWED "
            "entry if the code is genuinely not ported):\n  "
            + "\n  ".join(f"{fort[k]}  ->  expected {_twin(fort[k], _FORTRAN_EXTS, 'clubb_jax/src', '.py')}"
                          for k in missing),
            pytrace=False)


def test_no_unexplained_jax_only_files():
    """Every JAX source has a same-named Fortran file or an entry explaining why not.

    Fails with the list of JAX files that have neither, each with the Fortran path it would mirror.
    """
    fort, jax = _fortran_files(), _jax_files()
    listed = _keys(JAX_ONLY_ALLOWED_FILES)
    extra = sorted(k for k in jax.keys() - fort.keys() if k not in listed)
    if extra:
        pytest.fail(
            "JAX source with no same-named Fortran file (name it after the Fortran file whose code it holds, or "
            "merge it back into the JAX file mirroring that Fortran file; only add a JAX_ONLY_ALLOWED entry if "
            "there is no Fortran source for it at all):\n  "
            + "\n  ".join(f"{jax[k]}  ->  no {_twin(jax[k][len('clubb_jax/'):], ('.py',), 'src', '.F90')}"
                          for k in extra),
            pytrace=False)


def test_lists_are_still_live():
    """Every table entry still describes a real difference between the trees.

    A Fortran-side entry must name an existing Fortran file without a JAX twin, a JAX-side entry an existing
    JAX file without a Fortran twin, and an excluded directory must exist in src/ but not in clubb_jax/src/.
    Fails with each entry to remove, e.g. after its file was ported or deleted.
    """
    fort, jax = _fortran_files(), _jax_files()
    stale = []
    for key in sorted(FORTRAN_ONLY_ALLOWED_FILES):
        if key.lower() not in fort:
            stale.append(f"FORTRAN_ONLY_ALLOWED_FILES[{key!r}]: no such Fortran file any more -- remove the entry")
        elif key.lower() in jax:
            stale.append(f"FORTRAN_ONLY_ALLOWED_FILES[{key!r}]: a JAX file exists now -- remove the entry")
    for key in sorted(JAX_ONLY_ALLOWED_FILES):
        if key.lower() not in jax:
            stale.append(f"JAX_ONLY_ALLOWED_FILES[{key!r}]: no such JAX file any more -- remove the entry")
        elif key.lower() in fort:
            stale.append(f"JAX_ONLY_ALLOWED_FILES[{key!r}]: a Fortran file exists now -- remove the entry")
    for excluded in sorted(FORTRAN_ONLY_ALLOWED_DIRS):
        if not os.path.isdir(os.path.join(_FORT_ROOT, excluded)):
            stale.append(f"FORTRAN_ONLY_ALLOWED_DIRS[{excluded!r}]: no such Fortran directory -- remove it")
        elif os.path.isdir(os.path.join(_JAX_ROOT, excluded)):
            stale.append(f"FORTRAN_ONLY_ALLOWED_DIRS[{excluded!r}]: clubb_jax/src/{excluded} exists now, so the "
                         f"subtree is being ported -- drop the exclusion and compare it")
    if stale:
        pytest.fail("Table entries that no longer describe a real difference:\n  " + "\n  ".join(stale),
                    pytrace=False)
