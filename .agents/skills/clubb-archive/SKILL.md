---
name: clubb-archive
description: Preserve explicitly selected code, experiments or documentation on CLUBB's archive branch. Use for archiving requests; ordinary cleanup or deletion alone does not request archiving.
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

## Branch and entry layout

The branch is named `archive`. Its root contains its README and archive entry
directories, without an active CLUBB checkout. Work in a separate checkout or
repository so the active branch and uncommitted work remain intact. Inspect and
refresh an existing archive branch rather than recreating or overwriting it.
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
than fixing their existing whitespace. Builds and functional tests are not
needed merely to preserve an archive; they are separate work when requested.

Stage only intended payload and metadata. Nested historical `.gitignore` files
may hide selected files; explicitly include those files instead of staging local
artifacts indiscriminately. Commit and push when covered by the request and its
existing authorization. Retain other archive entries and concurrent remote work;
do not force-push or merge the archive into an active development branch.
Report the published branch/entry, origin information and preservation checks,
and distinguish archive publication from any still-uncommitted active-tree work.
