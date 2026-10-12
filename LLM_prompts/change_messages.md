# CLUBB change messages

Explain what changed, why it was worth doing, and what the results actually establish. Write for a colleague who knows CLUBB but has not followed the chat. Use the categories below to organize the information; do not turn a small change into a report.

## Establish the facts

Read the final diff and available run results/worklog when drafting for current work. For a rewrite, use the provided facts and original messages. Resolve the final scope, baseline, affected options, and tests before writing. An original claim is evidence of what the author reported, not a license to invent details. Review comments can contain important later results; commit messages do not necessarily summarize the final PR.

Use a short private fact list if helpful. Distinguish what ran, what was observed, what is inferred, and what is planned. Missing records do not prove a test was never run. In historical rewrites, say “No BFB result was recorded” or “The recorded comparison does not establish this path,” rather than “This path was not tested,” unless that absence is confirmed. If evidence is missing, say so in the relevant section and finish a useful draft. Do not invent a command, baseline SHA, compiler, case, duration, benchmark, test count, or approval. Use recorded test results directly in a draft for its author. Avoid repeating “reported” or “the supplied description”; keep provenance outside the message if needed. Describe an important missing comparison in plain language, such as “The recorded comparison does not establish BFB for the Luhar path,” rather than narrating every absent detail. Drafting/revising does not authorize posting, committing, rewriting history, or rerunning costly tests.

## Structure

**Title:** A single simple, specific statement of the change. Prefer a concrete action and object. Put the important option/trigger in the title if it distinguishes the bug. Leave secondary backends and maintenance changes for the body; the title does not have to name everything. Avoid branch names, “updates,” “improvements,” overloaded titles, and implementation jargon that belongs in the body. Do not force conventional-commit prefixes or an arbitrary character limit.

For a PR, use these headings by default:

- **Purpose:** What this does and the kind of change: bug fix, feature, organization, optimization, test coverage, etc. Explain behavior rather than merely naming the category.
- **Motivation:** The concrete problem or next goal. Explain how the old behavior caused the problem or why the change helps. A useful “why the test missed it” belongs here.
- **Result:** The resulting behavior, explicit BFB status and scope, and performance/progress when supported. Separate physical results from statistics, output schema, and workflow behavior.
- **Testing:** What ran and what it showed, with enough settings to understand the Result claims. State material failures or gaps plainly.
- **Additional notes:** Only when useful: tradeoff, next step, host impact, split-off work, or a rejected method that explains the choice.

A short PR may combine Purpose and Motivation into one opening paragraph, followed by Result and Testing. Commit messages use the title, blank line, and a compact body covering the same information; headings are optional. A commit describes its own change against its parent, while a PR describes the final change against its base. Do not transfer PR-wide test results to an earlier commit unless they apply to that revision. Keep useful issue/test/plot links beside the explanation they support.

Purpose says what this does, Motivation says why, Result says what was achieved, and Testing says how we know. Avoid explaining the mechanism or repeating the same evidence in every section. Most PRs need a sentence or short paragraph per category, with a few test bullets for a larger change. As a useful default, aim for about 120–250 words for a PR body and 60–140 for a commit body; use more when the scientific qualifications need it. These are guides, not limits. A rewrite should usually be shorter than an overgrown source. Retain decisive evidence, not every available fact.

## BFB and evidence

Always address BFB in Result, including in the compact commit body:

- **Verified BFB:** Name the comparison and tested scope. Cases, flags, columns, grid direction, compiler/optimization, device, duration, fields and baseline matter when they limit the claim. Report the settings you have; do not fabricate missing ones.
- **Known or expected bit changes:** Use the CLUBB marker `BIT_CHANGING` for changes to model numerical values or scientific output (including changed diagnostic fields or schema). Console wording, timer printouts, comments, and UI changes alone do not call for that marker; discuss their behavior and model BFB status separately. Explain which cases/options/columns/fields change, the mechanism, and why that change is acceptable. Distinguish observed changes from expected changes in an untested path. Default BFB can coexist with changes under another option or at `-O2`.
- **Not established:** Say “BFB is not established by the available comparisons” or describe the limited comparison. Use “not checked” only when that is known. A refactor or unchanged-looking arithmetic is not proof.
- **Not applicable:** Use only when there is no model/output comparison to make, such as documentation-only changes; explain briefly. Test or tooling changes that affect model execution/output still need a scoped status, or “not checked.”

“Matched,” close plots, successful runs, unit-test counts, and tolerance comparisons do not by themselves prove BFB. If a record says values or array leaves are “exactly equal” without identifying a bitwise comparison or explicitly establishing BFB, report that exact-value result without promoting it to verified bitwise equality. Say “within tolerance” and give the available tolerance. BFB at `-O0` does not prove BFB at `-O2`; CPU BFB does not prove GPU BFB. Identical physical fields with changed diagnostics or schema are not blanket BFB. Do not call tiny differences scientifically harmless without evidence.

For a bug fix, prefer evidence that reaches the affected path and catches the old bug. A passing default run does not establish an unused option or columns beyond the first. Name a demonstrated before-fails/after-passes test when available. Never claim that a regression test catches the bug just because it was added. Include the test's result, not merely that it exists.

Describe a failure as expected only when the evidence establishes that; explain what changed and which check detected it. Separate infrastructure/build failures from model comparison failures. Pending work stays pending. Performance needs measured before/after evidence with conditions; an old loop's share of runtime is not the achieved speedup. Label a quick observation or proposed benefit accordingly.

## Voice and final edit

Use ordinary causal language: “This happened because…”, “The old test only checked…”, “Moving this calculation lets us…”. Prefer concrete verbs such as fixes, moves, removes, runs, compares, skips. Keep technical names where they explain the trigger or mechanism. First person is fine when the author supplied it; do not invent personal anecdotes or imitate typos.

Avoid marketing claims, abstract noun chains, exhaustive file inventories, long validation dumps, defensive lists of unchanged things, and development chronology. Keep a limitation when it changes how the reader should interpret the result. Prefer the main regression/numerical test and a concise summary of other checks over every rerun, intermediate revision, inventory count, test-ownership detail, or validation command. Do not reproduce “merged latest master,” branch bookkeeping, or file-count assertions. Omit vague incidental “Python fix,” “small UI bug,” or “prompt update” titles unless you can explain a relevant behavior. When no tests are recorded, say “No test results were recorded”; do not enumerate all the tests someone might have run. Keep scientific comparison tolerances and the affected path; keep coverage counts only when they help assess scope. Do not add a negative caveat merely because the raw packet lacks a detail unrelated to the result claim. For mixed evidence, state the tradeoff instead of flattening it into “more stable,” “faster,” or “fully validated.”

Before returning the message, check: can a colleague explain the reason and outcome after one reading; does every Result claim have evidence in Testing; are bit changes and their acceptability clear; did you preserve material gaps; can any sentence be shortened or deleted without losing that understanding? Return the requested message, not the private fact list.
