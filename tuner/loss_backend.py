"""Select a loss runtime without importing the other backend."""


def load_loss_backend(backend):
    """Import only the backend selected by the tuning request."""
    if backend == "jax":
        from clubb_jax.src import clubb_loss_driver

        return clubb_loss_driver
    if backend == "fortran":
        from clubb_python import clubb_api

        return clubb_api
    raise ValueError(f"Unknown tuner backend: {backend}")


def parameter_metadata(backend):
    """Read the selected runtime's parameter names and hard bounds."""
    if backend == "jax":
        import numpy as np
        from clubb_jax.src.CLUBB_core.parameters_tunable import (
            get_param_names, get_parameter_hard_bounds,
        )

        names = get_param_names()
        bounds = get_parameter_hard_bounds()
        no_bound = np.finfo(bounds.dtype).max
        return names, [
            {
                "name": name,
                "min": float(bounds[0, i]) if bounds[0, i] > -no_bound else None,
                "max": float(bounds[1, i]) if bounds[1, i] < no_bound else None,
            }
            for i, name in enumerate(names)
        ]
    api = load_loss_backend(backend)
    names = list(api.get_param_names())
    return names, list(api.get_parameter_hard_bounds(len(names)))
