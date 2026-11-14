---
name: clubb-archive
description: Prepare a draft pull request to preserve explicitly selected code, experiments or documentation on CLUBB's archive branch. Use for archiving requests; ordinary cleanup or deletion alone does not request archiving.
---

# Archive selected material

Use the [archive branch README](https://github.com/larson-group/clubb/blob/archive/README.md)
as the maintained description of its layout and practices. Read its current
version before changing the archive, for example with
`git show origin/archive:README.md` after fetching that branch.

## Searching the archive

Exclude the archive branch, its checkouts and archived contents from routine
repository searches and code analysis. Read or search archived code only when
the user explicitly asks to consult that material. An archiving request permits
reading archive instructions, relevant entry metadata and the selected material
needed for that request; it does not invite broad exploration of historical code.

## Selection and scope

The user decides what to preserve. Archiving is manual: do not add deletion
hooks, scheduled jobs or automatic snapshot machinery, or establish a general
requirement to archive before deleting or merging anything.

Resolve the requested selection from the actual files and diff. Distinguish
retired files from relocations, and include explicitly selected previous
implementations even when their files were replaced rather than deleted.
Archiving does not itself authorize deleting the originals or publishing
unrelated cleanup work.

Retired tools, unfinished experiments and uncommitted or unversioned material
are all eligible. A source branch and commit are useful provenance, not
requirements. Do not require material to come from master, or create a source
branch/commit merely to qualify it for archiving.

## Review the selection through a pull request

Prepare the selected material and its entry documentation on a topic branch
based on the current `origin/archive`, then open a draft PR with `archive` as
its base. The diff is the reviewable selection: describe what is included,
why it is worth preserving, deliberate exclusions, provenance and preservation
checks. A request to archive authorizes preparing and publishing this PR;
it does not authorize merging it or pushing directly to `archive`.

Use [PRs and commits](../clubb-documentation/change-messages/SKILL.md) for the
message and draft state. Keep archive payload changes separate from active
CLUBB development branches. If the archive already contains the selected entry,
propose the update through the same PR workflow.

When archiving accompanies a removal, link the archive PR from the removal PR
description or the associated removal commit message. Link back to the source
or removal change from the archive PR and entry README. Check these references
before handing off either change, and state whether the archive is proposed or
merged.

If the user explicitly requests commits without PRs, follow that publication
scope. Reference the archive commit in the removal commit message, or in a
follow-up provenance commit if the removal was already committed; identify the
source/removal commit in the archive entry or commit message. Do not create a
PR or rewrite published history solely to add the connection.

## Branch and entry layout

The branch is named `archive`. Its root contains its README and archive entry
directories, without an active CLUBB checkout. Work on the archive topic branch
in a separate checkout or repository so the active branch and uncommitted work
remain intact.
Inspect and refresh the existing archive base rather than recreating or
overwriting it.
When explicitly asked to create it for the first time, start with independent
history and only the archive contents; do not copy the active tree into its root.

Name entries `MON_YEAR_SOURCE/`: uppercase three-letter English archive month,
four-digit archive year, and preferably the first seven hexadecimal characters
of a relevant source commit. For example, `NOV_2026_743b524/`. For uncommitted
work, a base commit can identify the context if the entry README clearly says
the files include uncommitted changes. With no useful commit, use a descriptive
label such as `NOV_2026_pdf_trial/`; never invent a hash or claim an exact snapshot.

Preserve the selected versions at their original or intended repository-relative
paths when known, including contents, hidden files, executable modes and symlinks.
For committed files, export from the chosen commit rather than a modified working
tree. For experiments, copy the selected working files and record their state.
Do not collect dependencies, generated outputs or caches unless requested.
Archived material has no functionality or maintenance requirement; do not repair,
format or refactor historical code to make it runnable.

## Entry documentation, including partial documents

Use [CLUBB documentation](../clubb-documentation/SKILL.md) for authored READMEs
and change messages. Each entry README should explain what was saved, when and
why, with available origin information: original paths, repository, branch,
full source/base hashes, working-directory context or explicitly unknown origin.
Distinguish exact committed copies from uncommitted or unversioned material.
Related PRs/issues, a directory list and deliberate exclusions can help readers.
Removal from master is not required or assumed.

Keep whole retired documents at their original paths. For selected passages
removed from an otherwise retained README or other document, preserve the excerpt
and useful surrounding context in the archive entry's README. Identify the source
file/section and available source version, explain the removal or revision, and
label the text as historical rather than current instructions. Do not copy an
entire retained document just to save one paragraph. Verify quotations and source
line links against the actual original text.

An existing entry can be extended when newly selected material belongs to the
same snapshot or experiment. Keep distinct historical versions in separate
entries; context corrections need not change the saved payload. If requested,
keep a pointer in the active CLUBB README small and put archive details in the
archive README rather than duplicating them.

## Verify and publish within the request

Check the exact selected-file list and compare saved contents and modes with
the chosen commit or working sources. Check excerpt accuracy, provenance claims,
links and authored README formatting. Preserve historical file bytes rather
than fixing their existing whitespace. Archive PRs require no tests of any kind:
no builds, unit, integration, functional or Jenkins tests. Preservation checks
(contents, modes, selection and provenance) are still required; they verify the
archive, not whether historical code works. Do not add a test requirement or
run tests for an archive PR unless the user explicitly requests that separate
work.

Stage only intended payload and metadata. Nested historical `.gitignore` files
may hide selected files; explicitly include those files instead of staging local
artifacts indiscriminately. Commit and push the topic branch within the request,
then create or update its draft PR targeting `archive`. Refresh the archive base
and retain other entries and concurrent remote work; do not force-push the archive
branch or merge it into an active development branch. Verify the PR base, head,
selected payload and draft state. Report the PR link, proposed entry, origin and
preservation checks. Describe the material as proposed until that PR is merged;
keep any active-tree work separate. Merge only when explicitly requested.
