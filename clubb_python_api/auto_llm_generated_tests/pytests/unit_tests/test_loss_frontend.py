"""Loss CLI help must not initialize the compiled model.

Moved from the JAX loss checks: this contract belongs to the Fortran Python
API frontend, while JAX keeps an independent equivalent without F2PY imports.
"""

import sys

import pytest

from clubb_python import clubb_standalone_loss


@pytest.mark.parametrize("option", ["-h", "-help", "--help"])
def test_loss_entry_point_accepts_help_without_initializing_model(monkeypatch, option):
    def unexpected_run(*args):
        pytest.fail("A help request must not initialize CLUBB")

    monkeypatch.setattr(clubb_standalone_loss.clubb_api, "clubb_get_loss", unexpected_run)
    monkeypatch.setattr(sys, "argv", ["loss_driver", option])
    assert clubb_standalone_loss.main() == 0
