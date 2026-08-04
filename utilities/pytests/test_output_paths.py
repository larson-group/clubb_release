from pathlib import Path

import pytest

from utilities.output_paths import OUTPUT_ROOT, REPO_ROOT, resolve_output_dir


def test_resolve_output_dir_roots_bare_names_under_output():
    assert resolve_output_dir("test1") == OUTPUT_ROOT / "test1"
    assert resolve_output_dir("experiments/example") == OUTPUT_ROOT / "experiments" / "example"


def test_resolve_output_dir_preserves_explicit_output_relative_path():
    assert resolve_output_dir("output/test1") == REPO_ROOT / "output" / "test1"


def test_resolve_output_dir_preserves_absolute_path():
    absolute = Path("/tmp/clubb-output-test")
    assert resolve_output_dir(absolute) == absolute


def test_resolve_output_dir_defaults_to_output_root():
    assert resolve_output_dir(None) == OUTPUT_ROOT
    assert resolve_output_dir("") == OUTPUT_ROOT


def test_resolve_output_dir_rejects_relative_traversal():
    with pytest.raises(ValueError, match="may not contain"):
        resolve_output_dir("../outside")
