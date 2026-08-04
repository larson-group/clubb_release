"""Small declarative registry for dashboard semantic tab operations.

The broker, MCP adapter, and persistent bridge only know two generic verbs:
``inspect_dashboard`` and ``invoke_dashboard``.  Individual tabs register the
operations they support here, including a concise JSON-schema-like description
for an agent to inspect before acting.  This intentionally does *not* turn the
browser into an unrestricted component/DOM automation surface.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable


OperationHandler = Callable[[dict[str, Any]], dict[str, Any]]


@dataclass(frozen=True)
class OperationSpec:
    """One semantic operation exposed by a dashboard tab."""

    tab: str
    name: str
    description: str
    schema: dict[str, Any]
    handler: OperationHandler

    def describe(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "input_schema": dict(self.schema),
        }


@dataclass(frozen=True)
class TabSpec:
    """Metadata needed to present one top-level dashboard tab to an agent."""

    name: str
    title: str
    description: str

    def describe(self, operations: list[OperationSpec]) -> dict[str, Any]:
        return {
            "tab": self.name,
            "title": self.title,
            "description": self.description,
            "operations": [operation.describe() for operation in operations],
        }


_TABS: dict[str, TabSpec] = {}
_OPERATIONS: dict[tuple[str, str], OperationSpec] = {}


def register_tab(name: str, title: str, description: str) -> None:
    """Register or refresh a top-level tab's static agent-facing metadata."""
    value = str(name or "").strip()
    if not value:
        raise ValueError("tab name is required")
    _TABS[value] = TabSpec(value, str(title or value), str(description or ""))


def register_operation(
    tab: str,
    name: str,
    description: str,
    schema: dict[str, Any],
    handler: OperationHandler,
) -> None:
    """Register or refresh one safe semantic operation for a known tab."""
    tab_name = str(tab or "").strip()
    operation_name = str(name or "").strip()
    if tab_name not in _TABS:
        raise ValueError(f"register tab '{tab_name}' before registering operations")
    if not operation_name:
        raise ValueError("operation name is required")
    _OPERATIONS[(tab_name, operation_name)] = OperationSpec(
        tab_name,
        operation_name,
        str(description or ""),
        dict(schema or {}),
        handler,
    )


def describe_tabs(tab: str | None = None) -> dict[str, Any]:
    """Return all known tabs, or one tab's declared operation surface."""
    requested = str(tab or "").strip()
    if requested and requested not in _TABS:
        raise ValueError(f"unknown dashboard tab: {requested}")
    names = [requested] if requested else sorted(_TABS)
    return {
        "tabs": [
            _TABS[name].describe(
                sorted(
                    [operation for (operation_tab, _), operation in _OPERATIONS.items() if operation_tab == name],
                    key=lambda operation: operation.name,
                )
            )
            for name in names
        ]
    }


def invoke(tab: str, operation: str, arguments: dict[str, Any] | None = None) -> dict[str, Any]:
    """Invoke a registered operation after a small common payload check."""
    tab_name = str(tab or "").strip()
    operation_name = str(operation or "").strip()
    spec = _OPERATIONS.get((tab_name, operation_name))
    if spec is None:
        available = sorted(name for registered_tab, name in _OPERATIONS if registered_tab == tab_name)
        if tab_name not in _TABS:
            raise ValueError(f"unknown dashboard tab: {tab_name}")
        raise ValueError(
            f"unknown operation '{operation_name}' for {tab_name}; available: "
            + (", ".join(available) or "none")
        )
    payload = dict(arguments or {})
    required = set(spec.schema.get("required") or [])
    missing = sorted(name for name in required if name not in payload or payload[name] in (None, ""))
    if missing:
        raise ValueError("missing required argument(s): " + ", ".join(missing))
    return spec.handler(payload)
