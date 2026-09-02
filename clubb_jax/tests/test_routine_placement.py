#!/usr/bin/env python3
"""test_routine_placement.py — routine-level structural guard that every JAX routine lives in the file that
mirrors its Fortran home module.

`test_mirror_audit.py` already asserts the audit's MISPLACED=0, but that check is *file-level* and consults the
audit's documented-rename / fold allowlist machinery. This test is an INDEPENDENT, source-parsed, routine-level
cross-check: it parses every `subroutine`/`function` name straight out of this checkout's Fortran `.F90` sources, then for
every JAX `def` whose bare name (or `_jax`-stripped base name — the JAX differentiable-variant suffix) matches a
Fortran routine name, asserts the JAX file stem is among the
Fortran home-file stems — with a small EXPLICIT allowlist of the handful of deliberate file-renames (the
PDF-parameter organization + the `advance_clubb_to_end` split). If anyone adds a routine to the wrong
file (its name matches a Fortran routine that lives in a different module), this fails with no dependence on the
audit's fold rules. Source-grounded; the Fortran tree is required. (iter 528)
"""
import os
import re
import glob
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "../.."))

# Deliberate, documented file-renames: (jax_file_stem, fortran_home_stem) pairs whose name-exact routines are
# intentionally NOT in an identically-named file. Each is an accepted architectural choice (see DESIGN.md /
# mirror_audit.py rename map), NOT a drift:
#   * CLUBB_core/pdf_params.py           ↔ pdf_parameter_module.F90  (the pdf_parameter type + its alloc/zero/
#     pack/unpack/print helpers use a shorter JAX filename)
#   * advance_clubb_to_end.py            ↔ clubb_driver.F90          (the timestep loop split into its own file)
#   * clubb_case_initalization.py        ↔ clubb_driver.F90          (the host initialization/finalization split)
_DOCUMENTED_RENAMES = {
    ("pdf_params", "pdf_parameter_module"),
    ("advance_clubb_to_end", "clubb_driver"),
    ("clubb_case_initalization", "clubb_driver"),
}


def _fortran_routine_homes():
    fort_dir = os.path.join(_ROOT, "src")
    homes = {}
    for f in glob.glob(fort_dir + "/**/*.F90", recursive=True):
        stem = os.path.basename(f)[:-4].lower()
        for line in open(f, errors="ignore"):
            m = re.match(r"\s*(?:pure |elemental |recursive )*(?:subroutine|function)\s+([A-Za-z]\w*)", line)
            if m:
                homes.setdefault(m.group(1).lower(), set()).add(stem)
    return homes


def test_routines_live_in_their_fortran_home_file():
    homes = _fortran_routine_homes()
    assert homes, f"No Fortran routines found under {os.path.join(_ROOT, 'src')}"
    offenders = []
    for jf in glob.glob(_ROOT + "/clubb_jax/src/**/*.py", recursive=True):
        jstem = os.path.basename(jf)[:-3].lower()
        for line in open(jf, errors="ignore"):
            m = re.match(r"def ([A-Za-z]\w*)", line)
            if not m:
                continue
            name = m.group(1).lower()
            # Match the bare name AND the `_jax`-stripped base (the JAX-only differentiable-variant suffix), so the
            # large class of `<fortran_name>_jax` routines is placement-checked too, not just name-exact ones.
            base = name[:-4] if name.endswith("_jax") else name
            fhomes = homes.get(name) or homes.get(base)
            if not fhomes or jstem in fhomes:
                continue
            # match in a different file: allowed only if a documented rename covers it
            if any((jstem == js and fh in fhomes) for (js, fh) in _DOCUMENTED_RENAMES):
                continue
            offenders.append(f"{os.path.relpath(jf, _ROOT)}: def {m.group(1)} — Fortran home(s) {sorted(fhomes)}")
    assert not offenders, (
        "JAX routine(s) defined in a file that does NOT mirror their Fortran home module (move them, or add a "
        "documented rename to _DOCUMENTED_RENAMES if the relocation is intentional):\n  " + "\n  ".join(offenders))
    print(f"  routine placement: every name-exact AND `_jax`-suffixed JAX routine is in its Fortran home file "
          f"({len(_DOCUMENTED_RENAMES)} documented renames excepted)  PASS")


def test_documented_renames_are_all_live():
    """Guard the allowlist itself: every documented-rename pair must still correspond to a real (jax_file, fortran
    routine-in-that-home) situation, so a stale exception can't silently mask a future genuine misplacement."""
    homes = _fortran_routine_homes()
    assert homes, f"No Fortran routines found under {os.path.join(_ROOT, 'src')}"
    stale = []
    for jstem, fhome in _DOCUMENTED_RENAMES:
        jf = glob.glob(_ROOT + f"/clubb_jax/src/**/{jstem}.py", recursive=True)
        if not jf:
            stale.append(f"({jstem}, {fhome}): no JAX file {jstem}.py")
            continue
        defs = {re.match(r"def ([A-Za-z]\w*)", l).group(1).lower()
                for l in open(jf[0], errors="ignore") if re.match(r"def ([A-Za-z]\w*)", l)}
        if not any(name in defs and fhome in homes.get(name, set()) for name in defs):
            stale.append(f"({jstem}, {fhome}): no name-exact routine from {fhome} now lives in {jstem}.py")
    assert not stale, "stale entries in _DOCUMENTED_RENAMES (remove them):\n  " + "\n  ".join(stale)
    print(f"  documented-rename liveness: all {len(_DOCUMENTED_RENAMES)} rename exceptions still active  PASS")


def main():
    print("test_routine_placement:")
    test_routines_live_in_their_fortran_home_file()
    test_documented_renames_are_all_live()
    print("All routine-placement checks PASSED")


if __name__ == "__main__":
    main()
