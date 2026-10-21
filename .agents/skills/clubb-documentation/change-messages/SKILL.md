---
name: clubb-change-messages
description: Write or revise CLUBB commit messages and PR titles/descriptions, with problem-first explanations, test evidence, and scoped BFB claims when relevant.
---

# PRs and commits

Read the [shared documentation requirements](../SKILL.md) for source
inspection and verification. Directory placement does not load them implicitly.

Use the [language guide](../language/SKILL.md). Write for a CLUBB
colleague who has not followed the chat.

Read the diff and available results. Describe a commit against its parent and
a PR against its base; tests must apply to the revision described. Separate
observations, expectations, and pending work. Never invent results or settings.
Missing records do not prove a test was never run.

Use [testing and reporting](../../clubb-testing-and-reporting/SKILL.md) to choose
and verify test evidence. In PRs, lead with relevant recorded Jenkins build
results and direct links. A local run of the same test usually need not be
listed unless it adds coverage or rigor; explain what supplemental local checks
and pytests tested rather than only giving a total. Include unexpected failures
and material gaps. Apply the same evidence selection to concise commit messages.

## Draft PRs and agent comments

Agent PRs are always drafts. Create them in draft state and do not automatically
mark them ready for review. For the full audit, validation and publication
sequence, use [branch finalization](../../clubb-finalization/SKILL.md) when the
user asks to get a branch ready to merge; writing a message alone does not
require that workflow.

End every PR comment posted through a user's account with `-AGENT` on its own
line, using the actual agent name, for example `-codex`, `-claude` or `-kimi`.
This includes conversation comments, inline review comments and replies.
Attribution does not itself authorize posting a comment.

## Structure

**Title:** the simplest specific statement of the change. Name the affected
behavior or option when it distinguishes the problem. Leave secondary details
for the body; do not force conventional-commit prefixes.

Use these PR headings:

- **Motivation and Change:** the concrete problem or unmet need, its cause when
  known, and how the PR addresses it.
- **Result:** resulting behavior and measured performance or progress when
  relevant. Discuss BFB only when the change could plausibly affect results.
- **Testing:** what ran, what it showed, and material failures or gaps.
- **Additional notes:** optional tradeoffs, follow-up work, or context that
  explains a decision.

Open Motivation and Change with the problem and why it matters, then explain
the solution. For a refactor, describe the duplication, unclear ownership or
other maintenance problem before describing the code movement. Connect the
implementation details to how they address that problem.

A commit uses a title, blank line, and compact body covering the same facts;
headings are optional. Keep commits especially simple. Avoid repeating facts
between sections, development chronology, file inventories, and exhaustive
command logs. Put useful links beside the claims they support.

## BFB and evidence

Read the canonical [bit-change policy](../../clubb-bit-changing/SKILL.md)
before deciding whether to declare a numerical answer change or how to split
its commits. Apply its test-or-condition justification to commit and PR messages.
Discuss preservation evidence when relevant to numerical answers; otherwise
omit BFB statements, including “not applicable.”

When relevant, describe the status in Result:

- **Verified BFB:** identify the bitwise comparison, baseline, and tested scope.
  Include cases, options, compiler/optimization, backend/device, duration, and
  fields when they limit the claim. Do not generalize beyond that scope.
- **Declared answer change:** follow the [bit-change policy](../../clubb-bit-changing/SKILL.md).
  Name the affected test(s) or precise untested conditions and explain the
  intentional numerical change and why it is acceptable.
- **Not established:** say what the available comparisons establish and what
  remains unknown. Say “not checked” only when that is known.

Tolerance matches, plots, passing tests, and numeric equality alone do not
establish bitwise equality. Unchanged-looking arithmetic is not proof either.
Do not call small differences scientifically harmless without evidence.

Testing must support Result. Use branch-finalization/Jenkins results when
available: list actual passes/failures and skipped or pending checks. BFB
claims must follow comparisons that ran. For bug fixes, name evidence that
exercises the affected path; claim a test catches the old bug only when
demonstrated.
Describe material failures accurately; call them expected only with evidence.
Performance claims need measured before/after results under stated conditions.
Keep missing evidence explicit and pending work pending.
