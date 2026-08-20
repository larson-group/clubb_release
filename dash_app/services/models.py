"""Strict public request models for CLUBB domain operations."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)


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


class ScmRunRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    case: str = Field(pattern=r"^[A-Za-z0-9_]+$")
    stats_file: str = "standard_stats.in"
    config: str = Field(default="default", pattern=r"^[A-Za-z0-9_.-]+$")
    overrides: dict[str, str | int | float | bool] = Field(default_factory=dict)
    run_options: ScmRunOptions = Field(default_factory=ScmRunOptions)
    out_dir: str | None = Field(default=None, max_length=512)

    @field_validator("stats_file")
    @classmethod
    def _simple_stats_name(cls, value: str) -> str:
        if value == "none":
            return value
        if "/" in value or "\\" in value or value in {"", ".", ".."}:
            raise ValueError("stats_file must be a checked-in stats filename or 'none'")
        return value


class ScmRunBatchRequest(StrictModel):
    """One idempotent request for a shared-output group of SCM cases."""

    request_id: str = Field(min_length=8, max_length=128)
    cases: list[str] = Field(min_length=1, max_length=64)
    stats_file: str = "standard_stats.in"
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
        if value == "none":
            return value
        if "/" in value or "\\" in value or value in {"", ".", ".."}:
            raise ValueError("stats_file must be a checked-in stats filename or 'none'")
        return value


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


class TuneRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    cases: list[TuneCaseRequest | str] = Field(min_length=1)
    parameter_ranges: list[ParameterRange] = Field(min_length=1, max_length=24)
    fields: list[str] = Field(default_factory=list, max_length=128)
    config: str = Field(default="default", pattern=r"^[A-Za-z0-9_.-]+$")
    strategy: Literal["random", "resolve", "simann", "adam"] = "random"
    max_samples: int = Field(default=12, ge=1, le=1000)
    resolve_spacing: float = Field(default=0.1, gt=0.0)
    simann_max_iters: int = Field(default=200, ge=1, le=5000)
    simann_initial_temp: float = Field(default=1.0, gt=0.0)
    simann_final_temp: float = Field(default=1.0e-12, gt=0.0)
    adam_max_updates: int = Field(default=100, ge=1, le=5000)
    adam_learning_rate: float = Field(default=0.01, gt=0.0)
    adam_perturbation: float = Field(default=0.05, gt=0.0, le=0.5)
    adam_spsa_pairs: int = Field(default=2, ge=1)
    batch_size: int = Field(default=1, ge=1, le=64)
    max_workers: int = Field(default=1, ge=1)
    loss_mode: str | None = None
    aggregation_mode: str | None = None
    aggregation_weights: list[float] = Field(default_factory=lambda: [0.1, 0.4, 0.4, 0.1], min_length=4, max_length=4)
    time_window_aggregation_scope: Literal["overall", "by_case"] = "overall"
    overrides: dict[str, str | int | float | bool] = Field(default_factory=dict)


class ProfileArtifactRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    run_id: str = Field(min_length=8)
    variables: list[str] = Field(min_length=1)
    time_window: TimeWindow
    vertical_coordinate: Literal["zt_m"] = "zt_m"


class LeaderboardRerunRequest(StrictModel):
    request_id: str = Field(min_length=8, max_length=128)
    mode: Literal["window", "complete"] = "window"
    max_results: int = Field(default=16, ge=1, le=16)
