"""Preserve the loss source's routine inventory and JAX math boundaries."""
import ast
from pathlib import Path
import re

from clubb_jax.src import clubb_loss_driver as loss

ROOT = Path(loss.__file__).resolve().parents[2]


def test_loss_routine_inventory_and_order_match_source():
    source = (ROOT / 'src/clubb_loss_driver.F90').read_text()
    names = re.findall(r'^\s*(?:logical function|subroutine)\s+(\w+)\(', source, re.MULTILINE)
    target = ast.parse((ROOT / 'clubb_jax/src/clubb_loss_driver.py').read_text())
    assert [node.name for node in target.body if isinstance(node, ast.FunctionDef)] == names


def test_loss_profile_math_uses_jax():
    target = ast.parse((ROOT / 'clubb_jax/src/clubb_loss_driver.py').read_text())
    for function in target.body:
        if isinstance(function, ast.FunctionDef) and function.name in ('calculate_taylor_metrics', 'calculate_field_loss'):
            assert not any(isinstance(node, ast.Name) and node.id == 'np' for node in ast.walk(function))
