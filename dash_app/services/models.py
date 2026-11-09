"""Strict public request models for CLUBB domain operations."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from clubb_jax.run_jax import runtime_selection
from dash_app.shared.stats import validate_stats_selection


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)


class JaxDeviceOptions(StrictModel):
    """Saved device preferences, validated by the JAX launcher."""

    jax_gpu: str = Field(
        default="", max_length=128,
        description="Device identifier returned by JAX runtime discovery; empty inherits the environment.",
    )
    jax_xla_prealloc: bool | None = Field(
        default=None, strict=True,
        description="GPU memory preallocation override; null preserves the launcher default/environment.",
    )

    @field_validator("jax_gpu")
    @classmethod
    def _device_selection(cls, value: str) -> str:
        return runtime_selection("gpu", device=value)["device"]


class TimeWindow(StrictModel):
    """A physical averaging interval in seconds from the case time origin."""

    start_seconds: float = Field(ge=0.0, description="Window start in case seconds.")
    duration_seconds: float = Field(gt=0.0, description="Averaging duration in seconds.")


class CompileRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    debug: bool = True
    python_bindings: bool = False
    fresh: bool = False
    gptl: bool = False


class ScmRunOptions(StrictModel):
    """Bounded diagnostic controls accepted by the normal SCM runner."""

    max_iters: int | None = Field(default=None, ge=1, le=1_000_000)
    dt_main: float | None = Field(default=None, gt=0.0)
    dt_rad: float | None = Field(default=None, gt=0.0)
    tout: float | None = Field(default=None, gt=0.0)


class ScmRunRequest(JaxDeviceOptions):
    request_id: str = Field(min_length=8, max_length=128)
    case: str = Field(pattern=r"^[A-Za-z0-9_]+$")
    implementation: Literal["fortran", "python", "jax"] = "fortran"
    jax_profile: Literal["cpu", "gpu"] = "cpu"
    stats_file: str = "standard"
    config: str = Field(default="default", pattern=r"^[A-Za-z0-9_.-]+$")
    overrides: dict[str, str | int | float | bool] = Field(default_factory=dict)
    run_options: ScmRunOptions = Field(default_factory=ScmRunOptions)
    out_dir: str | None = Field(default=None, max_length=512)

    @field_validator("stats_file")
    @classmethod
    def _simple_stats_name(cls, value: str) -> str:
        return validate_stats_selection(value)


class ScmRunBatchRequest(JaxDeviceOptions):
    """One idempotent request for a shared-output group of SCM cases."""

    request_id: str = Field(min_length=8, max_length=128)
    cases: list[str] = Field(min_length=1, max_length=64)
    implementation: Literal["fortran", "python", "jax"] = "fortran"
    jax_profile: Literal["cpu", "gpu"] = "cpu"
    stats_file: str = "standard"
    config: str = Field(default="default", pattern=r"^[A-Za-z0-9_.-]+$")
    overrides: dict[str, str | int | float | bool] = Field(default_factory=dict)
    run_options: ScmRunOptions = Field(default_factory=ScmRunOptions)
    max_workers: int | None = Field(default=None, ge=1, le=64)
    out_dir: str | None = Field(default=None, max_length=512)

    @field_validator("cases")
    @classmethod
    def _unique_cases(cls, values: list[str]) -> list[str]:
        cleaned = [str(value).strip() for value in values if str(value).strip()]
        if len(cleaned) != len(set(cleaned)):
            raise ValueError("cases must not contain duplicates")
        return cleaned

    @field_validator("stats_file")
    @classmethod
    def _simple_stats_name(cls, value: str) -> str:
        return validate_stats_selection(value)


class TuneCaseRequest(StrictModel):
    case_name: str = Field(pattern=r"^[A-Za-z0-9_]+$")
    time_average_range: tuple[float, float] | None = None
    altitude_comparison_range: tuple[float, float] | None = None
    average_time_seconds: int | None = Field(default=None, gt=0)
    num_time_windows: int | None = Field(default=None, ge=1, le=128)


class ParameterRange(StrictModel):
    name: str = Field(pattern=r"^[A-Za-z]\w*$")
    # ``targets`` lets one sampled coordinate drive a required equal set of
    # physical namelist parameters (for example, C6rt = C6thl).  Omit it for
    # the ordinary one-name range; the broker normalizes that to ``[name]``.
    targets: list[str] | None = Field(default=None, min_length=1, max_length=16)
    min: float
    max: float


class TuneRequest(JaxDeviceOptions):
    grad: bool = False
    gradient_mode: Literal['auto', 'forward', 'reverse'] = 'auto'
    request_id: str = Field(min_length=8, max_length=128)
    backend: Literal["fortran", "jax"] = "fortran"
    jax_options: str = Field(default="cpu", max_length=512)
    cases: list[TuneCaseRequest | str] = Field(min_length=1)
    parameter_ranges: list[ParameterRange] = Field(min_length=1, max_length=24)
    fields: list[str] = Field(default_factory=list, max_length=128)
    config: str = Field(default="default", pattern=r"^[A-Za-z0-9_.-]+$")
    strategy: str = Field(default="random", pattern=r"^[a-z][a-z0-9_]*$")
    strategy_options: dict[str, Any] = Field(default_factory=dict)
    batch_size: int = Field(default=1, ge=1, le=64)
    max_workers: int = Field(default=1, ge=1)
    loss_weights: dict[str, float] = Field(default_factory=lambda: {"mean": 1.0, "std": 1.0, "corr": 1.0})
    aggregation_mode: str | None = None
    aggregation_weights: list[float] = Field(default_factory=lambda: [0.1, 0.4, 0.4, 0.1], min_length=4, max_length=4)
    overrides: dict[str, str | int | float | bool] = Field(default_factory=dict)

    @field_validator("loss_weights")
    @classmethod
    def _loss_weights(cls, value):
        from utilities.loss_metrics import normalize_loss_weights
        return normalize_loss_weights(value)

    @field_validator("jax_options")
    @classmethod
    def _runtime_options(cls, value: str) -> str:
        runtime_selection(value)
        return value

    @model_validator(mode='after')
    def _gradient_backend(self):
        """Validate declared options and mandatory gradient/backend requirements."""
        from tuner.algorithm_registry import apply_requirements, normalize_options
        from utilities.clubb_loss_runner import validate_gradient_backend
        self.strategy_options = normalize_options(self.strategy, self.strategy_options)
        request = {"strategy": {"name": self.strategy, "options": self.strategy_options},
                   "backend": self.backend, "grad": self.grad}
        apply_requirements(request)
        self.grad = request["grad"]
        validate_gradient_backend(self.backend, self.grad)
        return self


class ProfileArtifactRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    run_id: str = Field(min_length=8)
    variables: list[str] = Field(min_length=1)
    time_window: TimeWindow
    vertical_coordinate: Literal["zt_m"] = "zt_m"


class LeaderboardRerunRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    workspace_id: str | None = Field(default=None, pattern=r"^[A-Za-z0-9_.-]+$")
    revision_id: str | None = Field(default=None, pattern=r"^[A-Za-z0-9_.-]+$")
    mode: Literal["window", "complete"] = "window"
    max_results: int = Field(default=16, ge=1, le=16)
