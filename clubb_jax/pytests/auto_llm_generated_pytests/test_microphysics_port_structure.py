"""Structural audit of the source-shaped microphysics interface (cores exempt)."""
from utilities.output_paths import REPO_ROOT as _REPO_ROOT

ROOT = _REPO_ROOT
MODULES = sorted((ROOT / 'src/Microphys').glob('*.F90'))


def test_runtime_interface_has_no_fortran_fallback():
    for source in MODULES:
        code=(ROOT/'clubb_jax'/source.relative_to(ROOT).with_suffix('.py')).read_text()
        assert 'import clubb_api' not in code
        assert 'pure_callback' not in code
