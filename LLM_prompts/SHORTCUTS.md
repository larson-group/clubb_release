# LLM prompt shortcuts

Choose the matching workflow from the current request and conversation. Read
that workflow before substantial work; use only the applicable parts and state
material exclusions. Selection does not authorize edits, broader validation,
external writes, a goal or a token budget. Ask only when an unresolved choice
changes the outcome. Explicit user instructions take priority.

| Request | Guidance |
| --- | --- |
| **Update Host Models After CLUBB Changes**: a host-consumed Fortran/C API signature, semantics, public type/constant, generated wrapper, flag or configuration changed. E3SM is excluded unless explicitly requested. | [Host compatibility](update_host_models_after_clubb_changes.md) |
| **Fix Python API**: repair F2PY/wrappers, update Python drivers after Fortran changes, make Python/JAX drivers match Fortran, or repair `run_python_vs_fortran_cases.py` / `run_jax_vs_fortran_cases.py`. | [API and drivers](update_python_api_and_drivers.md) |
| **Port Underlying Fortran**: port, re-port, audit or “similarize” a target file; align routines, calls, comments, argument lists; remove target-only helpers/aliases/optionals/reordered logic. | [Faithful port](port_underlying_fortran_to_other_languages.md) |
| **JAXize CLUBB Core File**: translate `src/CLUBB_core` into `clubb_jax/src`, replace a `clubb_api` call, wire `advance_clubb_core_module.py`, or continue a file port such as `advance_xm_wpxp`. | [JAX workflow](jax_porting_workflow.md), with the faithful-port standards |
| **Format Fortran Routines**: format interfaces/calls, extract a helper, document arguments/local variables, or apply routine formatting. | [Fortran formatting](fortran_routine_formatting.md) |
| **Refactor Column Loops**: move computation/column loops into routines, change multi-column interfaces or stats dimensions. | [Column refactor](fortran_column_refactor.md), with routine formatting |
| **Change Script Arguments**: options, forwarding, worker defaults, multicolumn settings, or Jenkins/documentation callers. | [CLI workflow](script_cli_conventions.md) |
| **Diagnose Numerical Comparisons**: JAX/Fortran, compiler, BFB or multi-column mismatches; distinguish regression, roundoff and untested behavior. | [Numerical diagnosis](numerical_regression_workflow.md) |
| **Use or Change Jenkins Tests**: launch named jobs, inspect failures/hangs, audit coverage, combine/remove/rename jobs or edit pipelines. | [Jenkins workflow](jenkins_workflow.md) |
| **Work Through CLUBB Dash**: connect/use/open/control Dash; compile/run, show profiles/contours/plots/console, tune, navigate pages, or create/update/view an investigation report. | [Dash workflow](dash_app_workflow.md) |

## Scope boundaries

- Host compatibility is host-owned work. Never modify, copy, synchronize,
  reformat or include vendored CLUBB/SILHS source in a host PR. Source sync is
  separate. Dash/MCP and Python-only changes do not trigger host work unless
  the host directly consumes that interface. Without provided host repos,
  audit/report exact host call-site changes; do not implicitly clone/edit them.
- API/driver work may target the Python API, Python driver, JAX driver or a
  combination. Distinguish compile/focused/smoke validation from a full suite.
- A faithful port uses current Fortran as authority. An explicitly requested
  idiomatic rewrite or relaxed numerical core is a different scope; it does
  not relax unrelated interfaces. Determine whether the request is discussion,
  an audit or edits, and whether it covers one file or multiple mirrors.
- For JAX work, the JAX workflow is primary and the faithful-port standards
  supply source matching. Existing/external JAX files are secondary references;
  reading this index does not create a file-level goal or require cloning.
- Fortran formatting may be formatting only, a behavior-preserving refactor,
  or a wider caller/wrapper update. Column refactors must distinguish BFB work
  from separately authorized behavioral fixes.
- CLI spelling applies to CLUBB commands, not external flags or structured APIs.
- Dash needs one unambiguous live browser view and matching broker. Read its
  discovery/permission/fallback rules; do not silently start a replacement or
  guess between instances. A named Jenkins request uses Jenkins.

[Maintaining guidance](README.md) describes how to record approved lessons
without duplicating requirements. Optional repo skills `clubb-jax-port` and
`clubb-script-cli` point to the same canonical workflows; they aid discovery
for agents with skill support, while this index supports other agents.
