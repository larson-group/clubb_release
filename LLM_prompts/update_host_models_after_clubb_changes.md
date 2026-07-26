# Update host models after CLUBB interface changes

Use this workflow when a CLUBB API, configuration, or Fortran interface change has merged into the CLUBB default branch and host-model repositories must be updated.

## Scope

The caller should provide the host repositories and the CLUBB change to port. Treat each repository as an independent change. Do not copy unrelated work from a more advanced checkout.

## Workflow

1. Read the repository instructions and inspect `LLM_prompts/SHORTCUTS.md` before making substantial changes.
2. Verify GitHub CLI availability and authentication with `gh --version` and `gh auth status`.
3. Clone each requested host repository and `larson-group/sys_admin` into dedicated temporary or workspace directories. Never reuse a dirty checkout without inspecting it first.
4. Identify the merged CLUBB change's removed symbols, changed routine signatures, generated wrappers, configuration fields, and tests. Search each host repository for all old names before editing.
5. Inspect `sys_admin`'s CLUBB synchronization configuration before touching copied files. In particular, check `gitUpdateScripts/run_updates.py` and the sync mappings for each host. Treat every mapped destination as subtree-owned, including parameter/flag and stats directories; do not edit those destinations manually because the sync process will overwrite them and may trigger notifications.
6. Update only the host-model compatibility code. Preserve host-specific behavior and formatting. Do not edit or copy files under vendored `CLUBB_core` or `SILHS` subtree directories; those sources are synchronized separately. Prefer the host's interface layer, especially routines and wrappers with `_api` suffixes. Edit configuration or namelist files only when `sys_admin` confirms they are host-owned rather than synchronized. Do not reintroduce a removed CLUBB flag merely to preserve an obsolete configuration path; if the CLUBB behavior is now unconditional, remove the host-side flag from the owning source/interface instead.
7. For WRF, check whether the interface update adds any source files, modules, or new dependencies. If so, update the relevant WRF build-system inputs and test that the new dependency is discoverable. If no files or dependencies are added, record that no build-system change is needed.
8. Run focused source scans and the most relevant available build, interface, or unit tests for each host repository. Record missing dependencies or unavailable full builds explicitly.
9. Before staging, compare the changed-file list against the `sys_admin` mappings and verify that no synchronized subtree/configuration destination is included.
10. For each repository, fetch the latest default branch and create a branch named `agent/<short-description>` from it. Stage only the intended files and commit with a concise message.
11. Open a draft pull request targeting the repository's default branch. The PR body should state:
   - the CLUBB change being synchronized;
   - the host-side files and interfaces updated;
   - behavior preserved or intentionally changed;
   - validation performed and any limitations.
12. If a PR branch must be refreshed after its target branch advances, fetch the target branch, rebase onto the latest remote default branch, rerun focused validation, and push with `--force-with-lease`. Use this only with explicit authorization because it rewrites the remote branch history.
13. Report each repository's branch, commit, PR URL, changed files, validation result, and any synchronized files intentionally left untouched. Do not merge the PRs.

## Safety

- Inspect `git status -sb` before staging every repository.
- Never use `git add -A` in a mixed or dirty checkout.
- Never force-push, reset hard, delete branches, or merge without explicit permission.
- If a repository is not writable or GitHub authentication is unavailable, stop before making external changes and report the blocker.
- If the same host change is needed in multiple repositories, implement and validate it independently; do not assume identical layouts.
