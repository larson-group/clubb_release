---
name: backlog
description: Track findings, TODOs, and out-of-scope issues discovered while working in this repo that shouldn't be fixed in the current session. Use when the user says "add to backlog", "log this for later", "don't fix now but remember", or similar -- and proactively when Claude notices a real, verified issue that is out of scope for the current task. Also use when the user asks "what's on the backlog", "check backlog", "list backlog items", or asks to work/fix a specific backlog item.
---

# Backlog

A lightweight running list of findings surfaced during work in this repo that are real but out of scope for the current task -- bugs, stale docs/comments, TODOs, cleanup opportunities. Not a replacement for a real issue tracker; just a way to not lose things discovered mid-task without derailing the task itself.

Backed by a single file: `BACKLOG.md` in this skill's directory (`.agents/skills/backlog/BACKLOG.md`; `.claude/skills` links to `.agents/skills`). It is gitignored and local to each checkout, so never commit it, and don't cite its entries from code comments or commit messages.

## Adding an entry

Trigger: user says "add to backlog" / "log this for later" / "don't fix now but remember this", OR Claude itself notices a real, independently-verified issue that is out of scope for the current task and would otherwise just get lost when the conversation ends.

Append a new `##` block to BACKLOG.md (create the file with a one-line title if it doesn't exist yet) using this template:

```
## <short title>

- Severity: <1-10 -- your judgment; 1 = trivial/cosmetic, 10 = correctness-critical or blocking>
- Found: <YYYY-MM-DD>, during <one-line context, e.g. "file-name-parity allowlist audit">
- Finding: <2-4 sentences: what's wrong or noteworthy, concrete enough to act on later>
- Verify: <one concrete instruction for re-checking this is still true later -- a grep command, a file:line to re-read, or a test to run>
```

Keep entries factual and specific -- cite file:line, quote the relevant code/docstring, don't editorialize. This is a flat log, not a priority queue; the severity field carries the priority signal, no other sorting/tagging needed. Every entry present in the file is implicitly open -- there is no "resolved" status (see "Removing an item" below): once something is done, it's gone, not marked.

## Listing / triaging

Trigger: "what's on the backlog", "check backlog", "list backlog items".

Read BACKLOG.md and summarize open entries (title + severity), sorted by severity descending, most severe first. Don't re-verify every entry just to list them -- only re-verify when the user asks to act on a specific one (see below).

## Working an item

Trigger: user asks to fix, investigate, or otherwise act on a specific backlog entry (by name, or by picking one from a list you just gave).

ALWAYS re-verify before acting, even if the entry looks recent -- code may have changed since it was logged:

1. Re-run the entry's `Verify` instruction (re-read the cited file/line, re-run the grep, re-run the test).
2. If the finding no longer holds (already fixed, code changed, file deleted, superseded) -- say so, then remove the entry (see "Removing an item"). Do not carry out the originally-described work against a premise that is no longer true.
3. If it still holds -- proceed with the fix/investigation at whatever scope the user now wants. Once it's actually done, remove the entry (see "Removing an item") -- don't leave it sitting around marked done.

## Removing an item

Trigger: an item was just fixed or confirmed stale (step 2/3 above), OR the user directly says "remove this from the backlog" / "drop that item" / "delete X from backlog".

Delete the entire `##` block for that entry from BACKLOG.md. Do not keep a "resolved"/"done"/"stale" graveyard -- the whole point of this file is a short, current list, not a growing archive. If the user wants a record of what was fixed, that belongs in a commit message or PR description, not here. If removing the last entry leaves the file empty, leave just the `# Backlog` title line.

## Notes

- Deliberately lightweight. Don't add process (extra priority dimensions, links to an external issue tracker, auto-triage rules) unless the user asks for it.
- If BACKLOG.md doesn't exist yet, create it with just a one-line title (`# Backlog`) -- the first appended entry establishes the real structure.
