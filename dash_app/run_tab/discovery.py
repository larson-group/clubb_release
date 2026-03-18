"""Filesystem and config discovery helpers for the run tab."""

import ast
import os

from .state import CASE_DIR, RUN_SCM_ALL, STATS_DIR


def list_cases():
    """Return available SCM case names from input/case_setups."""
    cases = []
    if not os.path.isdir(CASE_DIR):
        return cases
    for entry in os.listdir(CASE_DIR):
        if entry.endswith("_model.in"):
            cases.append(entry[: -len("_model.in")])
    return sorted(cases)


def load_available_cases():
    """Public wrapper for available run cases."""
    return list_cases()


def list_stats_files():
    """Return available stats namelist filenames from input/stats."""
    files = []
    if not os.path.isdir(STATS_DIR):
        return files
    for entry in os.listdir(STATS_DIR):
        if entry.endswith(".in"):
            files.append(entry)
    return sorted(files)


def load_stats_choices():
    """Public wrapper for available stats-file choices."""
    return list_stats_files()


def load_case_groups(available_cases):
    """Parse run_scm_all.py case groups and keep only cases present in this repo."""
    groups = {
        "all": [],
        "standard": [],
        "priority": [],
        "minimum": [],
        "short": [],
    }
    if not os.path.isfile(RUN_SCM_ALL):
        return groups
    try:
        with open(RUN_SCM_ALL, "r", encoding="utf-8") as handle:
            tree = ast.parse(handle.read(), filename=RUN_SCM_ALL)
    except Exception:
        return groups
    mapping = {
        "ALL_CASES": "all",
        "STANDARD_CASES": "standard",
        "PRIORITY_CASES": "priority",
        "MIN_CASES": "minimum",
        "SHORT_CASES": "short",
    }
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in mapping:
                try:
                    value = ast.literal_eval(node.value)
                except Exception:
                    continue
                if isinstance(value, list):
                    groups[mapping[target.id]] = [case for case in value if case in available_cases]
    return groups
