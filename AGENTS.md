# Agent Instructions

Before substantial work, select the relevant repo skills below and read their
`SKILL.md` files. Apply only the workflows needed for the task. Determine scope
from the current request and conversation; ask when an unresolved choice changes
the outcome. Skill selection does not itself authorize source edits, external writes, goals, or
broader validation. Explicit user instructions take priority.

Use only the applicable parts and state material exclusions. Selecting a skill
does not authorize a token budget.

Before drafting documentation or change messages, including small edits, read
[CLUBB documentation](.agents/skills/clubb-documentation/SKILL.md) and the
relevant writing guidance. Inspect the implementation or diff, including
helpers, callers, or argument parsers needed to support the explanation. Do
not infer behavior solely from names or existing comments.

## Skills

Related skills may be nested under an entry point. Search recursively for
`SKILL.md` when listing or validating them. Follow explicit links for shared
requirements; selecting a child does not automatically load its parent.

| Request | Skill |
| --- | --- |
| **Write Documentation or Change Messages**: create or revise READMEs, file and procedure descriptions, code comments, PR titles/descriptions, or commit messages. | [Documentation](.agents/skills/clubb-documentation/SKILL.md) |
| **Update Host Models After CLUBB Changes**: a host-consumed Fortran/C API signature, semantics, public type/constant, generated wrapper, flag or configuration changed. E3SM is excluded unless explicitly requested. | [Host compatibility](.agents/skills/clubb-host-compatibility/SKILL.md) |
| **Fix Python API**: repair F2PY/wrappers, update Python drivers after Fortran changes, make Python/JAX drivers match Fortran, or repair `run_python_vs_fortran_cases.py` / `run_jax_vs_fortran_cases.py`. | [API and drivers](.agents/skills/clubb-python-api/SKILL.md) |
| **Port Underlying Fortran**: port, re-port, audit or “similarize” a target file; align routines, calls, comments, argument lists; remove target-only helpers/aliases/optionals/reordered logic. | [Faithful port](.agents/skills/clubb-fortran-port/SKILL.md) |
| **JAXize CLUBB Core File**: translate `src/CLUBB_core` into `clubb_jax/src`, replace a `clubb_api` call, wire `advance_clubb_core_module.py`, or continue a file port such as `advance_xm_wpxp`. | [JAX workflow](.agents/skills/clubb-jax-port/SKILL.md), with the faithful-port standards |
| **Format Fortran Routines**: format interfaces/calls, extract a helper, document arguments/local variables, or apply routine formatting. | [Fortran formatting](.agents/skills/clubb-fortran-format/SKILL.md) |
| **Refactor Column Loops**: move computation/column loops into routines, change multi-column interfaces or stats dimensions. | [Column refactor](.agents/skills/clubb-column-refactor/SKILL.md), with routine formatting |
| **Change Script Arguments**: options, forwarding, worker defaults, multicolumn settings, or Jenkins/documentation callers. | [CLI workflow](.agents/skills/clubb-script-cli/SKILL.md) |
| **Diagnose Numerical Comparisons**: JAX/Fortran, compiler, BFB or multi-column mismatches; distinguish regression, roundoff and untested behavior. | [Numerical diagnosis](.agents/skills/clubb-numerical-diagnosis/SKILL.md) |
| **Run or Inspect Jenkins Tests**: launch named jobs, monitor builds, inspect failures/hangs or query live jobs. | [Jenkins access](.agents/skills/clubb-jenkins/SKILL.md) |
| **Create or Change Jenkins Tests**: add/edit pipelines, audit coverage, combine/remove/rename tests or coordinate their live jobs. | [Jenkins test authoring](.agents/skills/clubb-jenkins/test-authoring/SKILL.md) |
| **Declare Bit-Changing Commits**: decide whether to mark a Fortran answer change, identify affected BFB tests/options, or split preserving and changing edits. | [Bit-change policy](.agents/skills/clubb-bit-changing/SKILL.md) |
| **Choose Validation or Report Test Results**: plan relevant checks or summarize evidence in a PR, commit or response. | [Testing and reporting](.agents/skills/clubb-testing-and-reporting/SKILL.md) |
| **Finalize a Branch for Merge**: audit all branch changes, confirm current validation, stop for substantive unaccepted issues, then commit, push and create/update a draft PR. | [Finalization](.agents/skills/clubb-finalization/SKILL.md) |
| **Add or Review Pytests**: focused checks, generated-test admission, test pruning, or pytest CI coverage. | [Pytest workflow](.agents/skills/clubb-pytests/SKILL.md) |
| **Learn From a Session**: review recent corrections/decisions and improve their AGENTS, skill, shortcut or documentation owner. | [Session retrospective](.agents/skills/clubb-session-retrospective/SKILL.md) |
| **Work Through CLUBB Dash**: connect/use/open/control Dash; compile/run, show profiles/contours/plots/console, tune, navigate pages, or create/update/view an investigation report. | [Dash workflow](.agents/skills/clubb-dash/SKILL.md) |

## Repository conventions

- A request to investigate a problem authorizes diagnosis and a proposed remedy,
  not a substantive fix. Explain the findings and obtain acceptance before
  fixing, unless the user explicitly authorized fixing that problem or working
  through issues.
- Use [CLUBB documentation](.agents/skills/clubb-documentation/SKILL.md) for
  READMEs, change messages, file and procedure descriptions, and code comments,
  including documentation needed by new or changed code. It selects the
  relevant writing and language guidance.
- Fortran is the authority for faithful ports. Passing tests does not replace
  the routine/branch/signature/comment/formatting audit in the porting workflow.
- For script/CLI work, use
  [clubb-script-cli](.agents/skills/clubb-script-cli/SKILL.md): single-dash
  descriptive options, `output_*`, and `-workers` for process concurrency with
  a half-logical-CPU default. Keep worker counts distinct from batch width and
  preserve external-tool flags and documented GPU/plotting exceptions.
- Follow [clubb-pytests](.agents/skills/clubb-pytests/SKILL.md): fast specific
  checks belong in `pytests/`; new agent coverage goes in
  `auto_llm_generated_pytests/` pending human review. Real cases/full workflows
  belong in their component's `tests/`; reserve root `tests/` for shared
  workflows and keep Jenkins additions within the requested scope.
- Reuse the existing utility/parser/state owner. Keep UI/runner adapters thin;
  justify a new abstraction by actual shared behavior. Remove obsolete paths
  in the authorized change instead of adding silent fallback/retry machinery.
- Write user-facing documentation around purpose, a usable quick start, then
  needed detail. Update relevant sections of the existing document; keep implementation detail
  in its owner and define scientific symbols/units when they matter.
- Use [testing and reporting](.agents/skills/clubb-testing-and-reporting/SKILL.md)
  when choosing validation or summarizing results. Prefer relevant Jenkins build
  links in PRs; include local checks when they add coverage or rigor, and explain
  what they tested rather than relying on counts.
- Compare failures with the starting revision before attributing them to a
  change. Report supported, tested and untested behavior accurately.
- Agent pull requests are always drafts; do not automatically mark them ready
  for review. End PR comments posted through a user's account with an agent
  signature on its own line, such as `-codex`, `-claude` or `-kimi`. See
  [PRs and commits](.agents/skills/clubb-documentation/change-messages/SKILL.md).
- Preserve unrelated work and machine-local artifacts. Inspect the diff before
  staging; do not include ignored worklogs, outputs or local symlinks in commits.

## Session Worklog

Use repository-root `worklog.md` as the durable record for substantial work.
This is a local-only record excluded from normal Git staging by `.gitignore`;
create the file if it is not already present.
Read recent entries before continuing an investigation. Append a dated milestone
with the objective, decisions/assumptions, files, verification/results and remaining work;
never rewrite earlier history. Keep it concise and exclude secrets and large
raw output. This file is ignored: updating it does not mean committing it.

For an approved durable convention discovered during work, update its narrow
canonical guidance as part of the authorized change, or record an explicit
follow-up if documentation changes are outside scope. Keep task-specific choices
in the worklog; do not promote a question or one-off exception to a repo rule.
