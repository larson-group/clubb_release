"""Standing guard: no DEAD top-level public functions in clubb_jax/src.

Complements test_no_dead_imports (dead imports) by catching dead DEFINITIONS — a top-level `def` that is:
  * public (no leading underscore),
  * NOT a mirror of a Fortran routine (its base name, sans `_jax`, is not a Fortran subroutine/function — those
    are deliberately kept for the name-mirror even when unused, e.g. pack_pdf_params_api),
  * NOT exported via a module `__all__` (public API surface),
  * and has ZERO call/reference sites anywhere in src + tests (only its own def line).

Such a function is dead code with no Fortran equivalent — exactly the "remove jax-only routines that don't have a
direct Fortran equivalent" class. This guard caught (and now prevents the re-accumulation of) the f2py-wrapper-era
pack/unpack code removed iter 383. Found 0 at iter 394.

Pure Python; parses the current checkout's Fortran routine names directly.
"""
import os
import re
import glob

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
# Functions intentionally defined-but-uncalled-in-repo (e.g. an external entry point). Empty at iter 394.
_ALLOW = set()

_DEF = re.compile(r"^def ([a-zA-Z][a-zA-Z_0-9]*)", re.M)


def _fortran_routines():
    """Return names of routines declared in this checkout's Fortran sources."""
    routines = set()
    pattern = re.compile(
        r"\s*(?:pure |elemental |recursive )*(?:subroutine|function)\s+([A-Za-z]\w*)",
        re.IGNORECASE,
    )
    for path in glob.glob(os.path.join(_ROOT, "src", "**", "*.F90"), recursive=True):
        with open(path, errors="ignore") as source:
            for line in source:
                match = pattern.match(line)
                if match:
                    routines.add(match.group(1).lower())
    return routines


def _dead_functions():
    fort = _fortran_routines()
    corpus = ""
    for d in ("clubb_jax/src", "clubb_jax/tests"):
        for f in glob.glob(os.path.join(_ROOT, d, "**", "*.py"), recursive=True):
            corpus += open(f, errors="ignore").read()
    dead = []
    for f in glob.glob(os.path.join(_ROOT, "clubb_jax", "src", "**", "*.py"), recursive=True):
        if "__" in os.path.basename(f):
            continue
        txt = open(f, errors="ignore").read()
        allblock = re.search(r"__all__\s*=\s*\[(.*?)\]", txt, re.S)
        exported = set(re.findall(r"['\"]([a-zA-Z_0-9]+)['\"]", allblock.group(1))) if allblock else set()
        for name in _DEF.findall(txt):
            if name.startswith("_"):
                continue
            base = name[:-4] if name.endswith("_jax") else name
            if base.lower() in fort or name in exported:
                continue
            if (os.path.basename(f), name) in _ALLOW:
                continue
            if len(re.findall(r"\b" + re.escape(name) + r"\b", corpus)) <= 1:
                dead.append(f"{os.path.relpath(f, _ROOT)}: {name}")
    return dead


def _dead_private_helpers():
    """Module-level PRIVATE (`_`-prefixed) defs that are jax-only internals (never mirror Fortran) and have ZERO
    references anywhere in src+tests beyond their own def line. `test_no_dead_functions` deliberately
    scopes to PUBLIC defs (API surface); this complements it for the file-local helpers, which are exactly the
    "jax-only routines with no direct Fortran equivalent" the instruction targets for removal once dead.

    Excludes: dunder defs; and **decorator-registered** defs (a `def` whose immediately-preceding line begins with
    `@` — e.g. a `@x.defjvp` custom-derivative rule like `_safe_div_jvp`, registered by the decorator rather than
    called by name). Indented (nested) helpers are out of scope — only module-level `^def _…`."""
    corpus = ""
    for d in ("clubb_jax/src", "clubb_jax/tests"):
        for f in glob.glob(os.path.join(_ROOT, d, "**", "*.py"), recursive=True):
            corpus += open(f, errors="ignore").read()
    dead = []
    for f in glob.glob(os.path.join(_ROOT, "clubb_jax", "src", "**", "*.py"), recursive=True):
        if "__" in os.path.basename(f):
            continue
        lines = open(f, errors="ignore").read().split("\n")
        for i, ln in enumerate(lines):
            mt = re.match(r"def (_[a-zA-Z][a-zA-Z_0-9]*)\b", ln)   # module-level (unindented) private def
            if not mt:
                continue
            name = mt.group(1)
            if name.startswith("__"):
                continue
            if i > 0 and lines[i - 1].lstrip().startswith("@"):   # decorator-registered (jvp/vjp/etc.)
                continue
            if (os.path.basename(f), name) in _ALLOW:
                continue
            if len(re.findall(r"\b" + re.escape(name) + r"\b", corpus)) <= 1:
                dead.append(f"{os.path.relpath(f, _ROOT)}: {name}")
    return dead


def test_no_dead_functions():
    fort_src = os.path.join(_ROOT, "src", "CLUBB_core")
    assert os.path.isdir(fort_src), f"Fortran source directory is missing: {fort_src}"
    dead = _dead_functions()
    assert not dead, (
        "Dead (never-called) public functions with no Fortran mirror found — remove them, or if intentional "
        "(external entry point) add to _ALLOW in this test:\n  " + "\n  ".join(sorted(dead)))
    print("  no dead public functions in clubb_jax/src (mirror routines + __all__ exports excluded)  PASS")


def test_no_dead_private_helpers():
    dead = _dead_private_helpers()
    assert not dead, (
        "Dead (never-referenced) module-level private helpers found — remove them, or if intentional add to "
        "_ALLOW in this test:\n  " + "\n  ".join(sorted(dead)))
    print("  no dead private helpers in clubb_jax/src (decorator-registered rules excluded)  PASS")


if __name__ == "__main__":
    print("test_no_dead_functions:")
    test_no_dead_functions()
    test_no_dead_private_helpers()
    print("Done.")
