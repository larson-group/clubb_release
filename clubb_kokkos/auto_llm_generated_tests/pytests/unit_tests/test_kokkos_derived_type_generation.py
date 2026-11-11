"""Keep handwritten grid APIs when regenerating the native type-field mirror.

The real generator reads the current Fortran field schema. These checks change
that schema in a temporary output tree and verify that only the marked type
region changes; a missing marker must fail without overwriting the public API.
"""

import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "clubb_kokkos/src/f2k"))
import generate_derived_types as generator


def test_kokkos_grid_api_survives_type_regeneration(tmp_path, monkeypatch):
    """Update a native field while preserving every handwritten API byte and mtime on a no-op."""
    source = "CLUBB_core/grid_class"
    target = tmp_path / (source + ".hpp")
    target.parent.mkdir(parents=True)
    original = (ROOT / "clubb_kokkos/src" / (source + ".hpp")).read_text()
    target.write_text(original)
    fields = generator.parse_type("grid", source)
    fields[0] = replace(fields[0], name="new_grid_dimension")
    monkeypatch.setattr(generator, "KOKKOS", tmp_path)
    monkeypatch.setattr(generator, "SOURCES", {"grid": source})
    generator.emit_headers({"grid": fields})
    actual = target.read_text()
    begin = "// BEGIN GENERATED GRID TYPE"
    end = "// END GENERATED GRID TYPE"
    assert actual.split(begin)[0] == original.split(begin)[0]
    assert actual.split(end)[1] == original.split(end)[1]
    assert "int new_grid_dimension{};" in actual
    assert "int nzm{};" not in actual
    before = target.stat().st_mtime_ns
    generator.emit_headers({"grid": fields})
    assert target.read_text() == actual
    assert target.stat().st_mtime_ns == before


@pytest.mark.parametrize("marker", ["BEGIN", "END"])
def test_kokkos_grid_missing_marker_preserves_header(tmp_path, monkeypatch, marker):
    """Reject a damaged generated-region marker without deleting any handwritten API."""
    source = "CLUBB_core/grid_class"
    target = tmp_path / (source + ".hpp")
    target.parent.mkdir(parents=True)
    original = (ROOT / "clubb_kokkos/src" / (source + ".hpp")).read_text()
    original = original.replace(f"// {marker} GENERATED GRID TYPE", "// missing marker")
    target.write_text(original)
    monkeypatch.setattr(generator, "KOKKOS", tmp_path)
    monkeypatch.setattr(generator, "SOURCES", {"grid": source})
    with pytest.raises(ValueError, match="marker pair"):
        generator.emit_headers({"grid": generator.parse_type("grid", source)})
    assert target.read_text() == original
