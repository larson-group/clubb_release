"""Check that JAX default ConfigFlags values match the Fortran defaults."""


def test_default_values_match_fortran():
    """The JAX default ConfigFlag VALUES match the Fortran defaults, field-for-field, so a case NOT setting a flag
    gets the SAME default behavior as the Fortran — a diverged default would silently change results for every case
    that leaves the flag unset. Compares `get_default_config_flags()` against the Fortran
    `clubb_api.get_default_config_flags()`."""
    from clubb_python import clubb_api
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
