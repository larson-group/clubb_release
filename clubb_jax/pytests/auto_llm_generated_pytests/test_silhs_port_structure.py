"""Fortran routine/signature contract, with explicit native-JAX adaptations."""

from utilities.output_paths import REPO_ROOT as ROOT

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
