"""Standing code-hygiene guard: no orphaned (dead) imports in clubb_jax/src.

Dead imports silently accumulated (~58) from the iter-160 whole-driver extractions before the iters
336-341 sweep removed them. This guard prevents recurrence: it flags any module-level import whose bound
name is never loaded (precise `ast.Name`-Load analysis — not grep), excluding `from __future__`, names
re-exported via `__all__`, and a small explicit ALLOWLIST of deliberate keeps.

Pure-Python (stdlib `ast` only), so it never SKIPs and runs fast.
"""
import ast
import os
import glob

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))

# Deliberate keeps (NOT cruft) — verified iter 341/342:
#   mpace_a `grav`                          : carries `# noqa: F401` (intentional)
#   adg1 `calc_comp_corrs_binormal`         : documented re-export mirroring the Fortran `use pdf_utilities`
# (iter 392: dropped the grid_class one/zero/one_half entries — those re-exported constants had no
#  consumer and were removed from grid_class.py, so they no longer need an allowance.)
_ALLOW = {
    ("mpace_a.py", "grav"),
    ("adg1_adg2_3d_luhar_pdf.py", "calc_comp_corrs_binormal"),
}


def _dead_imports(path):
    tree = ast.parse(open(path, errors="ignore").read())
    loads = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}
    strings = {n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)}
    imported = {}
    for n in ast.walk(tree):
        if isinstance(n, ast.ImportFrom):
            if n.module == "__future__":
                continue
            for a in n.names:
                if a.name != "*":
                    imported[a.asname or a.name] = n.lineno
        elif isinstance(n, ast.Import):
            for a in n.names:
                imported[(a.asname or a.name).split(".")[0]] = n.lineno
    return [k for k in imported if k not in loads and k not in strings]


def test_no_dead_imports():
    offenders = []
    for f in sorted(glob.glob(os.path.join(_SRC, "**", "*.py"), recursive=True)):
        base = os.path.basename(f)
        if base.startswith("__"):
            continue
        for name in _dead_imports(f):
            if (base, name) in _ALLOW:
                continue
            offenders.append(f"{f.split('/src/')[-1]}: {name}")
    assert not offenders, (
        "Dead (never-loaded) imports found — remove them, or if intentional (re-export / API-mirror / noqa) "
        "add to _ALLOW in this test:\n  " + "\n  ".join(offenders))
    print(f"  no dead imports in clubb_jax/src ({len(_ALLOW)} deliberate keeps allowlisted)  PASS")


def _fortran_runtime_refs(path):
    """AST refs to the compiled Fortran oracle module `clubb_python` in a src file — live `import clubb_python` /
    `from clubb_python ...` or any `clubb_python` name use. AST-based, so docstrings/comments mentioning it (the
    historical 'fallback removed iter 388' notes) are ignored — only EXECUTABLE references count."""
    tree = ast.parse(open(path, errors="ignore").read())
    refs = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                if a.name.split(".")[0] == "clubb_python":
                    refs.append(n.lineno)
        elif isinstance(n, ast.ImportFrom):
            if (n.module or "").split(".")[0] == "clubb_python":
                refs.append(n.lineno)
        elif isinstance(n, ast.Name) and n.id == "clubb_python":
            refs.append(n.lineno)
    return refs


def test_src_has_no_fortran_runtime_import():
    """Standing guard for the port's core 'runs 100% in JAX, zero Fortran calls per timestep' property: the compiled
    Fortran oracle `clubb_python` must NOT be imported or referenced anywhere in `clubb_jax/src` (it is the comparison
    ORACLE — used only by f2py tests, which SKIP when it is unbuilt). The old lazy `clubb_python.clubb_api` fallbacks
    were removed iters 388/389 and replaced with fail-loud raises. This asserts no EXECUTABLE `clubb_python` reference
    has crept back into the runtime (AST-based, so the historical mentions in comments/docstrings are ignored). A new
    src-side `clubb_python` import fails loudly — the port would no longer be Fortran-free. (iter 590)"""
    offenders = []
    for f in sorted(glob.glob(os.path.join(_SRC, "**", "*.py"), recursive=True)):
        if os.path.basename(f).startswith("__"):
            continue
        for ln in _fortran_runtime_refs(f):
            offenders.append(f"{f.split('/src/')[-1]}:{ln}")
    assert not offenders, (
        "clubb_jax/src now has an EXECUTABLE reference to the Fortran oracle `clubb_python` — the runtime is no longer "
        "100% JAX / zero-Fortran-calls. Remove the import/call (or replace with a fail-loud raise as in iter 389):\n  "
        + "\n  ".join(offenders))
    print("  clubb_jax/src is free of executable clubb_python (Fortran-oracle) references — runtime is 100% JAX  PASS")


if __name__ == "__main__":
    print("test_no_dead_imports:")
    test_no_dead_imports()
    test_src_has_no_fortran_runtime_import()
    print("Done.")
