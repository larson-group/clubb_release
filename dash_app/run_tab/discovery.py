"""Filesystem and config discovery helpers for the run tab."""

from input.scripts.create_case_namelist import available_case_names


def list_cases():
    """Discover catalog and saved cases through the shared case resolver."""
    return available_case_names()


def load_available_cases():
    """Public wrapper for available run cases."""
    return list_cases()


def load_case_groups(available_cases):
    """Resolve catalog run lists for the Run tab's established selection buttons."""
    from input.scripts.input_catalog import resolve_selection
    mapping = {"all": "all_cases", "standard": "standard_cases", "priority": "priority_cases",
               "minimum": "min_cases", "short": "short_cases"}
    return {label: [case for case in resolve_selection("cases", name) if case in available_cases]
            for label, name in mapping.items()}
