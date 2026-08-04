# Agent Instructions

Before beginning substantial work, check `LLM_prompts/SHORTCUTS.md` for reusable prompt guidance that may match the request.

If a shortcut appears relevant:

1. Read the referenced prompt.
2. Confirm that its goal and constraints match the user's current request.
3. Use it as guidance only if it fits the current task.
4. If it partially fits, apply only the relevant parts and state what does not apply.

These prompts are guidance, not automatic commands. The user's current request takes priority.

## Session Worklog

Use the repository-root `worklog.md` as the durable session record for
substantial work.

- Read the recent entries before continuing an existing investigation.
- Append a dated entry after each meaningful completed change or diagnostic
  milestone; do not rewrite or discard earlier session history.
- Record the objective, important decisions and assumptions, files changed,
  verification commands/results, and any remaining follow-up work.
- Keep entries concise enough to scan but specific enough that a later agent
  can resume without reconstructing the session from chat history.
- Do not place secrets, credentials, or large generated outputs in the worklog.
