---
name: clubb-change-messages
description: Write or revise CLUBB commit messages and PR titles/descriptions, with clear purpose, test evidence, and scoped BFB claims when relevant.
---

# PRs and commits

Use the [language guide](../clubb-doc-language/SKILL.md). Write for a CLUBB
colleague who has not followed the chat.

Read the diff and available results. Describe a commit against its parent and
a PR against its base; tests must apply to the revision described. Separate
observations, expectations, and pending work. Never invent results or settings.
Missing records do not prove a test was never run.

## Structure

**Title:** the simplest specific statement of the change. Name the affected
behavior or option when it distinguishes the problem. Leave secondary details
for the body; do not force conventional-commit prefixes.

Use these PR headings:

- **Purpose:** what changes and what kind of change it is.
- **Motivation:** the concrete problem or goal, and why this change helps.
- **Result:** resulting behavior and measured performance or progress when
  relevant. Discuss BFB only when the change could plausibly affect results.
- **Testing:** what ran, what it showed, and material failures or gaps.
- **Additional notes:** optional tradeoffs, follow-up work, or context that
  explains a decision.

A small PR may combine Purpose and Motivation. A commit uses a title, blank
line, and compact body covering the same facts; headings are optional. Keep
commits especially simple. Avoid repeating facts between sections, development
chronology, file inventories, and exhaustive command logs. Put useful links
beside the claims they support.

## BFB and evidence

Discuss BFB only when a change could plausibly affect model results or
scientific output. Renames, documentation, or high-level Dash/tuner/PyPlotGen
changes usually need no BFB statement. A possible build failure or crash alone
does not qualify. Otherwise omit BFB, including “not applicable.”

When relevant, describe the status in Result:

- **Verified BFB:** identify the bitwise comparison, baseline, and tested scope.
  Include cases, options, compiler/optimization, backend/device, duration, and
  fields when they limit the claim. Do not generalize beyond that scope.
- **`BIT_CHANGING`:** use for known or expected changes to model numerical
  values or scientific output, including diagnostics or output schema. Explain
  what changes, why, and why it is acceptable. Distinguish observed changes
  from expected changes in an untested path. Console text and timers alone do
  not require this marker.
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
