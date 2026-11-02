"""CLUBB services used by Dash callbacks, MCP tools and numerical workers.

Public request and job types load on demand so importing a numerical service
in a model worker does not load dashboard dependencies.
"""

from importlib import import_module

__all__ = [
    "ArtifactStore",
    "CompileRequest",
    "JobStore",
    "LeaderboardRerunRequest",
    "ParameterRange",
    "ProfileArtifactRequest",
    "ScmRunBatchRequest",
    "ScmRunOptions",
    "ScmRunRequest",
    "SubmissionConflict",
    "TimeWindow",
    "TuneCaseRequest",
    "TuneRequest",
]


def __getattr__(name):
    """Load a public job or request type only when its caller requests it."""
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    owner = ".jobs" if name in {"ArtifactStore", "JobStore", "SubmissionConflict"} else ".models"
    value = getattr(import_module(owner, __name__), name)
    globals()[name] = value
    return value
