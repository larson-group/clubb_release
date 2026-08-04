# Update host models after CLUBB interface changes

Use this workflow when a CLUBB API, configuration, or Fortran interface change has merged into the relevant CLUBB branch and host-model repositories need host-side compatibility updates.

This workflow updates host-owned code only. It must never update the CLUBB or SILHS source bundled inside a host model.

## Scope

The caller should provide the host repositories and the CLUBB change to port. Treat each repository as an independent change. Do not copy unrelated work from a more advanced checkout. A host-model PR must not modify, copy, merge, subtree-sync, reformat, or regenerate any vendored CLUBB/SILHS path. Examples include `SRC/CLUBB`, `WRF/phys/clubb`, `src/physics/clubb/src/CLUBB_core`, SILHS directories, and every destination listed as synchronized by `sys_admin`.

If the host's vendored CLUBB/SILHS source also needs synchronization, report that as a separate prerequisite or separate workflow. Do not perform that synchronization, or include its files, in this host-compatibility PR.

## Workflow

1. Read the repository instructions and inspect `LLM_prompts/SHORTCUTS.md` before making substantial changes.
2. Resolve the exact target branch for each host before touching files. This is usually `master`, but may be an integration branch such as CAM's `clubb_silhs_devel`; use the branch that the PR is meant to target, not an assumed default.
3. Use a fresh dedicated clone when possible. If an existing checkout is used, inspect `git status -sb` and its remote, but treat it as potentially stale even when clean. Do not assume a checkout from a previous task is current.
4. Verify that `origin` points to the intended canonical host-model repository, not another local checkout. Fetch the actual remote target branch before comparison or branch creation: `git fetch origin --prune <target>`. Verify the fetched `origin/<target>` SHA against `git ls-remote origin refs/heads/<target>`, and create the work branch from `origin/<target>`, never from a possibly stale local `<target>` branch. Record the base SHA.
5. Identify the CLUBB change's removed symbols, changed routine signatures, generated wrappers, configuration fields, and tests. Search each host repository for all old names before editing.
6. Inspect `sys_admin`'s CLUBB synchronization configuration only to identify forbidden synchronized destinations. In particular, check `gitUpdateScripts/run_updates.py` and the sync mappings for each host. Do not run the synchronization or edit/copy any mapped CLUBB/SILHS destination as part of this workflow.
7. Update only host-owned compatibility code. Preserve host-specific behavior and formatting. Prefer the host's interface layer, especially routines and wrappers with `_api` suffixes. Edit configuration or namelist files only when `sys_admin` confirms they are host-owned rather than synchronized. Do not reintroduce a removed CLUBB flag merely to preserve an obsolete configuration path; if CLUBB behavior is now unconditional, remove the host-side flag from the owning host interface.
8. For WRF, check whether the interface update adds any host-owned source files, modules, or new dependencies. If so, update the relevant WRF build-system inputs and test that the new dependency is discoverable. Never add or update vendored CLUBB/SILHS source through this step. If no host files or dependencies are added, record that no build-system change is needed.
9. Run focused source scans and the most relevant available build, interface, or unit tests for each host repository. Record missing dependencies or unavailable full builds explicitly.
10. Before staging, compare the changed-file list against the fetched remote target branch and verify that every changed file is intentionally host-owned. Fail closed if any CLUBB/SILHS vendored or `sys_admin`-mapped path appears. Also run `git diff --check`.
11. For each repository, create a branch named `agent/<short-description>` from the freshly fetched `origin/<target>` branch. Stage only the intended files and commit with a concise message.
12. Before opening a PR, run `git diff --name-status origin/<target>...HEAD` and verify the exact file list, then confirm the PR base is the same remote target branch. Do not open the PR if the list contains an unexpected file or if the base advanced; fetch and rebase the work branch first.
13. Open a draft pull request targeting the exact fetched target branch. The PR body should state:
   - the CLUBB change being synchronized;
   - the host-side files and interfaces updated;
   - behavior preserved or intentionally changed;
   - validation performed and any limitations.
14. If a PR branch must be refreshed after its target branch advances, fetch the target branch, rebase onto the latest remote target branch, rerun focused validation, and push with `--force-with-lease`. Use this only with explicit authorization because it rewrites the remote branch history.
15. Report each repository's branch, base SHA, commit, PR URL, changed files, validation result, and synchronized files intentionally left untouched. Do not merge the PRs.

## Safety

- Inspect `git status -sb` before staging every repository.
- Use the available authenticated Git remote and GitHub API/integration when publishing or inspecting PRs. If authentication is unavailable, stop before making external writes.
- Never use `git add -A` in a mixed or dirty checkout.
- A clean but stale checkout is still unsafe; verify freshness against the remote target SHA before creating a branch.
- Treat any unexpected changed file, especially anything under a vendored CLUBB/SILHS or `sys_admin`-mapped path, as a blocking scope error and stop before staging or opening the PR.
- Never force-push, reset hard, delete branches, or merge without explicit permission.
- If a repository is not writable or GitHub authentication is unavailable, stop before making external changes and report the blocker.
- If the same host change is needed in multiple repositories, implement and validate it independently; do not assume identical layouts.
