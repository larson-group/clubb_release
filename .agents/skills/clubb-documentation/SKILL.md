---
name: clubb-documentation
description: Route requests to write or edit CLUBB READMEs, change messages, file descriptions, or code comments.
---

# CLUBB documentation

Read the skill for the text being written or revised:

| Text | Skill |
| --- | --- |
| README or substantial feature's usage documentation | [READMEs](../clubb-readmes/SKILL.md) |
| PR title/description or commit message | [PRs and commits](../clubb-change-messages/SKILL.md) |
| Opening file description or module docstring | [File descriptions](../clubb-file-descriptions/SKILL.md) |
| Procedure descriptions and comments within code | [Procedure and code comments](../clubb-code-comments/SKILL.md) |
| Wording and how much context to explain | [Language](../clubb-doc-language/SKILL.md) |

Each writing skill uses the language guide. Read the relevant skills, rather
than the whole set. A file-wide comment rewrite uses both file descriptions
and procedure/code comments. A request just to run or summarize code does not
require a documentation edit.

## Before drafting

For source documentation, also review the applicable language or porting
standards listed in AGENTS.md: [Fortran formatting](../clubb-fortran-format/SKILL.md),
[JAX porting](../clubb-jax-port/SKILL.md), and
[faithful ports](../clubb-fortran-port/SKILL.md) for Fortran-derived code.
Read and apply the parts relevant to the documentation task, including comment
placement, section markers, and dividers, within the requested scope.

Inspect the implementation before writing. Read the complete target routine
and the helpers or callers needed to understand its purpose, inputs, outputs,
defaults, and important boundary or error behavior. For usage documentation,
inspect the actual entry points and argument parsers. For change messages,
read the relevant diff and available test records. Existing comments and names
are clues, not evidence of behavior. If a claim cannot be verified, omit it or
state the uncertainty.

## Before finishing

Check the draft's factual claims against the implementation or recorded
results. Make a short checklist from the relevant writing guide and check each item.
For file descriptions, check the language's header form, purpose, public-routine
list, and usage/understanding details. For READMEs, check named entry points,
examples, and first-mention links. For procedures, check purpose, inputs,
outputs, and relevant references. Check named options and paths against their actual owners.
Repeat this verification after revisions; a wording or formatting change can
introduce a new factual error. Keep comment/docstring syntax valid and preserve
code during documentation-only edits.
