"""Apply compiler policy before loading the CLUBB JAX model.

Importing any model submodule configures the CPU workaround through
clubb_jax.backends.cpu.configure_environment, also used by the launcher.
The compiler defect is documented at https://github.com/jax-ml/jax/issues/42185.
CLUBB_JAX_ACCELERATOR selects the backend; CPU is the default. Import the model
before initializing a JAX backend, or restart an already initialized process.
"""

import os

from clubb_jax.backends import cpu

if os.environ.get("CLUBB_JAX_ACCELERATOR", "cpu").strip().lower() == "cpu":
    cpu.configure_environment(os.environ)
