# Learn from a recent session

Use this workflow when the user asks to review a recent session and improve
`AGENTS.md`, skills or shortcuts by recording lessons. Start with
[maintaining guidance](README.md) and [SHORTCUTS.md](SHORTCUTS.md). A review-only
request remains read-only; apply edits when the user requests them.

1. Read the relevant session messages, corrections, decisions and resulting
   code. Identify repeated instructions, follow-up fixes, approved conventions,
   and workflows that required unnecessary manual explanation. Distinguish a
   human-approved durable requirement from an agent suggestion, a question,
   an experiment or a task-specific exception. Use current owner code/docs to
   resolve stale examples; do not infer standards solely from what an agent did.
2. For each durable lesson, locate the existing owner. Decide whether the rule
   is missing, contradictory, hard to discover, or already present but not
   checked before completion. Prefer an edit to the narrow canonical workflow
   or owner README. Keep root `AGENTS.md` short, use the shortcut index for
   routing, and make a repo skill a thin entry point when discovery benefits.
   Keep one authoritative copy of each detailed requirement.
3. Describe the practical rule, its trigger and reason, applicable exceptions,
   and observable completion evidence. When agents ignored an existing rule,
   improve the finishing procedure rather than repeating that rule everywhere.
   Preserve all still-applicable requirements when simplifying or reorganizing;
   map moved or superseded rules to their new owner or the later user decision.
4. Check links, current command examples, scope and portability. Keep credentials,
   machine paths, host ports, temporary environment details and fixed model choices
   out of reusable guidance. Put one-off decisions in the ignored worklog.
5. When authorized, make a small coherent guidance change and summarize the
   lessons, exact owners edited, verification and limits. If the change is
   substantial, try a realistic ambiguous request or a small implementation
   trial with expectations set in advance; inspect actual evidence. Use
   independent subagent trials only when delegation is authorized. Do not
   claim a small prompt trial proves reliable execution on a full scientific port.

Do not automatically create an automation, goal, new chat or broad research
project. Identify uncertain lessons and proposals explicitly; do not encode
unapproved preferences as settled repository requirements.
