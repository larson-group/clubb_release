"""Fortran routine/signature contract, with explicit native-JAX adaptations."""

import ast
from utilities.output_paths import REPO_ROOT as ROOT
import pytest
from clubb_jax.pytests.port_structure_helpers import (
    fortran_routines,
    python_routines,
)

SOURCES = sorted(p for p in (ROOT / "src/SILHS").glob("*.F90") if p.stem != "mt95")

# Random operations require explicit keys; persistent permutations are a pytree.
ADAPTATIONS = {
    ("generate_uniform_sample_module", "rand_uniform_real"): ("key", "shape"),
    ("generate_uniform_sample_module", "generate_uniform_lh_sample"): (
        "key",
        "sampling_state",
    ),
    ("generate_uniform_sample_module", "choose_permuted_random"): ("key",),
    ("generate_uniform_sample_module", "permute_height_time"): ("key",),
    ("generate_uniform_sample_module", "rand_permute"): ("key",),
    ("latin_hypercube_driver_module", "generate_silhs_sample"): ("sampling_state",),
    ("latin_hypercube_driver_module", "generate_all_uniform_samples"): (
        "key",
        "sampling_state",
    ),
    ("math_utilities", "rand_integer_in_range"): ("key",),
    ("silhs_api_module", "generate_silhs_sample_api"): ("sampling_state",),
    ("silhs_importance_sample_module", "importance_sampling_driver"): ("key",),
    ("silhs_importance_sample_module", "cloud_weighted_sampling_driver"): ("key",),
    ("silhs_importance_sample_module", "generate_strat_uniform_variate"): ("key",),
    ("silhs_importance_sample_module", "choose_x_u_scaled"): ("key",),
}


def test_silhs_inventory_excludes_mt95():
    target = ROOT / "clubb_jax/src/SILHS"
    assert {p.stem for p in target.glob("*.py") if p.stem != "__init__"} == {
        p.stem for p in SOURCES
    }


@pytest.mark.parametrize("source", SOURCES, ids=lambda p: p.stem)
def test_silhs_source_routine_and_argument_order(source):
    target = ROOT / "clubb_jax" / source.relative_to(ROOT).with_suffix(".py")
    expected = fortran_routines(source)
    # Nested JAX scan/cond bodies express source DO/IF blocks, not new APIs.
    actual = {
        k: v for k, v in python_routines(ast.parse(target.read_text())).items() if len(k) == 1
    }
    assert list(actual) == list(expected)
    for key, function in actual.items():
        args = [a.arg.lower() for a in function.args.args]
        adaptations = ADAPTATIONS.get((source.stem, key[0]), ())
        assert tuple(a for a in args if a in adaptations) == adaptations
        args = [a for a in args if a not in adaptations]
        source_args, intents = expected[key]
        assert args == [a for a in source_args if a in args]
        assert not [a for a in source_args if intents.get(a) in ("in", "inout") and a not in args]
        assert not function.args.kwonlyargs
        assert function.args.vararg is None and function.args.kwarg is None
        if key != ("rand_uniform_real",):
            assert not function.args.defaults
