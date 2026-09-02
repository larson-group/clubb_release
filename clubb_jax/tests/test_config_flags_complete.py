"""Standing guard: the JAX `ConfigFlags` covers every Fortran *configurable* model flag.

The Fortran `clubb_config_flags` derived type carries all the namelist-settable model flags. If the JAX
`ConfigFlags` were missing one, a case_setup that sets that flag would be **silently ignored** (the JAX never
loads it) — a footgun. This test extracts the configurable flag names from model_flags.F90 (every
`clubb_config_flags%<field>` reference) and asserts each is a JAX `ConfigFlags` field. Verified complete
(67 of 67) iter 375.

NB this guards COVERAGE (the flag is loadable), not WIRING — whether each flag's behavior is implemented or
fail-loud guarded is enforced separately by `test_unsupported_config_guards.py`.

Pure-Python (reads the F90 source + the JAX namedtuple); SKIPs cleanly if the Fortran oracle source is absent.
"""
import os
import re
import sys

_HERE = os.path.dirname(__file__)
_ROOT = os.path.normpath(os.path.join(os.path.abspath(_HERE), "../.."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
# The AUTHORITATIVE list of case-settable flags is the namelist a case's *_model.in is read into:
# clubb_driver.F90 `namelist /configurable_clubb_flags_nl/` (iter 377 — was model_flags.F90's
# `clubb_config_flags%X` proxy, a superset that also counts non-case-settable internal references).
_DRIVER = os.path.abspath(os.path.join(
    _HERE, "..", "..", "src", "clubb_driver.F90"))


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
        print("  clubb_driver.F90 oracle source not present — SKIP")
        return
    from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags
    jax_fields = {f.lower() for f in get_default_config_flags()._fields}
    fort_flags = _fortran_configurable_flags()
    assert fort_flags, "extracted 0 namelist flags — the configurable_clubb_flags_nl extraction broke"
    missing = sorted(f for f in fort_flags if f not in jax_fields)
    assert not missing, (
        "JAX ConfigFlags is MISSING these case-settable flags (a model.in setting one would be silently "
        f"ignored): {missing}. Add them to CLUBB_core/config_flags.py + model_flags.py defaults.")
    print(f"  ConfigFlags covers all {len(fort_flags)} case-settable namelist flags  PASS")


def test_default_values_match_fortran():
    """The JAX default ConfigFlag VALUES match the Fortran defaults, field-for-field. Coverage (above) guards that a
    case CAN set a flag; this guards that a case NOT setting a flag gets the SAME default behavior as the Fortran — a
    diverged default would silently change results for every case that leaves the flag unset. Compares
    `get_default_config_flags()` against the Fortran `clubb_api.get_default_config_flags()`. SKIPs if clubb_python is
    unbuilt. (iter 466)"""
    for _p in (_ROOT, _ROOT + "/clubb_python_api"):
        if _p not in sys.path:
            sys.path.append(_p)
    try:
        from clubb_python import clubb_api
    except Exception as e:
        print(f"  config-flag defaults vs Fortran: SKIP ({type(e).__name__})")
        return
    from clubb_jax.src.CLUBB_core.model_flags import get_default_config_flags
    jflags = get_default_config_flags()
    fflags = clubb_api.get_default_config_flags()
    mism = []
    spurious = []
    checked = 0
    for name in jflags._fields:
        if not hasattr(fflags, name):
            # Reverse-direction guard: a JAX ConfigFlags field with no Fortran `clubb_config_flags` counterpart is a
            # misnamed/spurious flag the JAX reads but the oracle has no notion of — surface it, don't silently skip.
            spurious.append(name)
            continue
        jv = getattr(jflags, name)
        fv = getattr(fflags, name)
        checked += 1
        if isinstance(jv, bool):
            if bool(fv) != jv:
                mism.append(f"{name}: JAX {jv} vs Fortran {bool(fv)}")
        else:
            if fv != jv and str(fv) != str(jv):
                mism.append(f"{name}: JAX {jv} vs Fortran {fv}")
    assert checked > 0, "no overlapping flags compared — clubb_api default flags object unexpected"
    assert not spurious, ("JAX ConfigFlags has field(s) absent from the Fortran clubb_config_flags type "
                          f"(misnamed/spurious): {spurious}")
    assert not mism, ("JAX default ConfigFlags VALUES diverge from the Fortran defaults (cases leaving these unset "
                      "would behave differently):\n  " + "\n  ".join(mism))
    print(f"  default ConfigFlag values match Fortran field-for-field, no spurious JAX flag ({checked} flags)  PASS")


if __name__ == "__main__":
    print("test_config_flags_complete:")
    test_config_flags_cover_all_fortran_configurable_flags()
    test_default_values_match_fortran()
    print("Done.")
