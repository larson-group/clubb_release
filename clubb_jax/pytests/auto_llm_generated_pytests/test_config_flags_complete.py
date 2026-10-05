"""Standing guard: the JAX `ConfigFlags` covers every Fortran *configurable* model flag.

The Fortran `clubb_config_flags` derived type carries all the namelist-settable model flags. If the JAX
`ConfigFlags` were missing one, a case_setup that sets that flag would be **silently ignored** (the JAX never
loads it) — a footgun. This test extracts the configurable flag names from model_flags.F90 (every
`clubb_config_flags%<field>` reference) and asserts each is a JAX `ConfigFlags` field. Verified complete
(67 of 67) iter 375.

NB this guards COVERAGE (the flag is loadable), not WIRING — whether each flag's behavior is implemented or
fail-loud guarded is enforced separately by `test_unsupported_config_guards.py`.

Pure-Python (reads the F90 source + the JAX namedtuple); requires the checked-in Fortran source.
"""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT
import os
import re

_ROOT = str(_REPO_ROOT)
# The AUTHORITATIVE list of case-settable flags is the namelist a case's *_model.in is read into:
# clubb_driver.F90 `namelist /configurable_clubb_flags_nl/` (iter 377 — was model_flags.F90's
# `clubb_config_flags%X` proxy, a superset that also counts non-case-settable internal references).
_DRIVER = os.path.abspath(os.path.join(
    _ROOT, "src", "clubb_driver.F90"))


def _fortran_configurable_flags():
    """Names in the `configurable_clubb_flags_nl` namelist — exactly the flags a case's model.in can set."""
    txt = open(_DRIVER, errors="ignore").read()
    m = re.search(r"namelist\s*/\s*configurable_clubb_flags_nl\s*/(.*?)(?:\n\s*\n|\n\s*namelist\s*/)",
                  txt, re.S | re.I)
    if not m:
        return set()
    body = m.group(1)
    names = set()
    for line in body.splitlines():
        line = line.split("!", 1)[0].replace("&", " ")
        for n in re.findall(r"[a-zA-Z_][a-zA-Z_0-9]*", line):
            names.add(n.lower())
    return names


def test_config_flags_cover_all_fortran_configurable_flags():
    if not os.path.exists(_DRIVER):
        raise AssertionError("  clubb_driver.F90 oracle source not present")
    from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags
    jax_fields = {f.lower() for f in get_default_config_flags()._fields}
    fort_flags = _fortran_configurable_flags()
    assert fort_flags, "extracted 0 namelist flags — the configurable_clubb_flags_nl extraction broke"
    missing = sorted(f for f in fort_flags if f not in jax_fields)
    assert not missing, (
        "JAX ConfigFlags is MISSING these case-settable flags (a model.in setting one would be silently "
        f"ignored): {missing}. Add them to CLUBB_core/config_flags.py + model_flags.py defaults.")
    print(f"  ConfigFlags covers all {len(fort_flags)} case-settable namelist flags  PASS")
