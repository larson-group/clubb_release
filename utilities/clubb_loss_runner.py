"""Reusable CLUBB loss sessions for persistent backend workers.

clubb_loss_runner initializes a case and its benchmark once, evaluates successive
parameter batches, and closes the loss driver on context exit. The tuner and
explorer retain this object inside their backend processes; process transport
belongs to their existing worker owners. Importing this module loads no model.

Use init() before reading field_names, param_names, default_params or bounds.
evaluate(params, grad=False, gradient_parameters=None, gradient_mode="auto",
gradient_target="loss", profiles=False) returns
four named metrics, loss[candidate], field_loss[field,candidate] and optional
gradients/profiles. Metrics have axes (window,field,candidate). The default
gradient["loss"] has axes (candidate,parameter); explicit gradient_target="metrics" adds
a parameter axis to each diagnostic metric. Profile dictionaries
contain heights, units, benchmark (window, level) and model (window, candidate,
level). close() releases the session. One active loss session is allowed per
process because the underlying drivers own process-global state.
Reverse mode checkpoints each numerical timestep and replays it during
backpropagation. Forward/value-only evaluations use normal advancement.
"""
from pathlib import Path

GRADIENT_BACKENDS = frozenset({'jax'})
GRADIENT_MODES = ('auto', 'forward', 'reverse')
GRADIENT_TARGETS = ('loss', 'metrics')


def validate_gradient_mode(mode):
    """Require an explicit supported differentiation mode or automatic selection."""
    if mode not in GRADIENT_MODES:
        raise ValueError('gradient_mode must be auto, forward or reverse')
    return mode


def validate_gradient_target(target):
    """Choose the scalar case objective or the explicit explorer metric Jacobians."""
    if target not in GRADIENT_TARGETS:
        raise ValueError('gradient_target must be loss or metrics')
    return target


def resolve_gradient_mode(mode, num_parameters, num_outputs):
    """Compare parameters with one case output, or field groups for metric maps."""
    validate_gradient_mode(mode)
    return ('reverse' if num_parameters > num_outputs else 'forward') if mode == 'auto' else mode


def validate_gradient_backend(backend, grad):
    """Reject gradient requests for backends outside the explicit whitelist."""
    if not isinstance(grad, bool):
        raise ValueError('grad must be a boolean')
    if grad and backend not in GRADIENT_BACKENDS:
        raise ValueError('Loss gradients require the jax backend')


class clubb_loss_runner:
    """Own an initialized loss driver throughout a backend worker's lifetime."""

    def __init__(self, backend, case, fields, work_dir, *, batch_size=1,
                 case_defaults=None, loss_weights=None, field_weights=None,
                 config='default', override=None):
        from utilities.loss_metrics import normalize_loss_weights, normalize_named_weights
        self.backend, self.case, self.fields = backend, case, list(fields)
        self.loss_weights = normalize_loss_weights(loss_weights)
        self.field_weights = normalize_named_weights(field_weights, self.fields, 'field')
        self.work_dir = Path(work_dir)
        self.batch_size = batch_size
        self.case_defaults, self.config, self.override = case_defaults, config, override
        self.api = None

    def init(self):
        """Prepare the benchmark and model once; return this session."""
        if self.api is not None:
            return self
        from input.scripts.create_case_namelist import create_loss_case_namelist
        from tuner.loss_backend import load_loss_backend, parameter_metadata
        import numpy as np

        runfile, _, _ = create_loss_case_namelist(
            self.case, self.work_dir, self.fields, batch_size=self.batch_size,
            num_time_windows=(self.case_defaults or {}).get('num_time_windows', 1),
            case_defaults=self.case_defaults, config=self.config, override=self.override,
            loss_weights=self.loss_weights, field_weights=self.field_weights,
            duplicate_params_for_batch=True, disable_stats_storage=True)
        api = load_loss_backend(self.backend)
        self.field_names, params = api.init_clubb_loss(str(runfile), return_default_params=True)
        self.api = api
        self.default_params = np.asarray(params, dtype=float).copy()
        self.param_names, self.bounds = parameter_metadata(self.backend)
        return self

    def evaluate(self, params, *, grad=False, gradient_parameters=None, gradient_mode='auto', gradient_target='loss', profiles=False, progress=None):
        """Evaluate candidates from initial conditions and return metrics and requested data.

        Reverse gradients use one-timestep checkpoint replay. Forward gradients
        and value-only evaluations use ordinary advancement.
        """
        validate_gradient_backend(self.backend, grad)
        validate_gradient_mode(gradient_mode)
        validate_gradient_target(gradient_target)
        self.init()
        import numpy as np
        from utilities.loss_metrics import LOSS_METRIC_NAMES

        if self.backend == 'jax':
            result = self.api.evaluate_loss(params, grad=grad,
                                           gradient_parameters=gradient_parameters,
                                           gradient_mode=gradient_mode, gradient_target=gradient_target,
                                           progress=progress)
            if profiles:
                request = self.api.active_request
                writer = request.state['stats_writer']
                grids = dict(zip(('zt', 'zm'), writer.get_source_grid()))
                result['profiles'] = {}
                for field in request.fields:
                    grid, units, _ = writer.registry[field.clubb_var_name]
                    batches = [np.stack([np.asarray(self.api._field_profile(field, stats))[:count]
                                         for stats in windows])
                               for count, windows in request.state['_loss_stats_snapshots']]
                    result['profiles'][field.clubb_var_name] = {
                        'heights_m': np.asarray(grids[grid])[0, field.k_min:field.k_max + 1],
                        'units': units, 'benchmark': np.asarray(field.truth_profile),
                        'model': np.concatenate(batches, axis=1)}
        else:
            if profiles:
                raise ValueError('In-memory loss profiles currently require the jax backend')
            values = self.api.clubb_get_loss_for_params(params)
            result = {'metrics': dict(zip(LOSS_METRIC_NAMES, values[:4])), 'loss': values[4]}
        result['metrics'] = {name: np.asarray(values) for name, values in result['metrics'].items()}
        result['loss'] = np.asarray(result['loss'])
        from utilities.loss_metrics import root_mean_square
        result['field_loss'] = root_mean_square(result['metrics']['rmse'], axis=0)
        result['loss_weights'] = dict(self.loss_weights)
        result['field_weights'] = dict(self.field_weights)
        if 'gradient' in result:
            result['gradient'] = {name: np.asarray(values) for name, values in result['gradient'].items()}
        return result

    def close(self):
        """Release the active model and benchmark session."""
        if self.api is not None:
            self.api.finalize_clubb_loss()
            self.api = None

    def __enter__(self):
        return self.init()

    def __exit__(self, *_):
        self.close()
