# Maintaining agent guidance

Start with [SHORTCUTS.md](SHORTCUTS.md). It routes a request to the relevant
workflow; it is not a plan to execute every workflow. Root `AGENTS.md` contains
only common decisions and completion habits. Keep detailed requirements in a
narrow workflow, and keep current interface/runtime contracts in the owner
README and code. Repo skills are thin discovery adapters to those same files.

## Learning from a correction

When the user approves a durable convention, check whether it is missing,
contradicted, hard to discover, or present but not checked before completion.
Update the existing owner accordingly. If it was already present, strengthen
its completion evidence instead of copying the rule into several prompts.
Keep the reason and verification in the ignored worklog. A question, temporary
experiment or task-specific exception is not a new repository-wide requirement.

Prefer short descriptions with a clear trigger, an authority, important
exceptions and observable completion evidence. Preserve detailed source-faithful
port requirements; reduce repeated routing text instead. Avoid hard-coded host
ports, credentials, checkout paths, model choices and machine-local environment
names in reusable repo guidance.

## Reviewing a guidance change

Map existing requirements to retained, moved or explicitly superseded guidance.
Explain superseded requirements using the later user decision or current owner
contract. Do not silently delete an inconvenient rule. Check links and current
CLI examples against their owners; optional ignored files must remain optional.

Try the changed guidance on a realistic ambiguous request and, when useful,
a small independent implementation task. Decide expected outcomes and scope
boundaries before observing the responses. Check actual edits and evidence,
not just whether an agent says it followed the instructions. Record misses and
limits; a toy evaluation is not proof of numerical correctness for a real port.
Do not add a permanent benchmark framework for every wording edit.

For a user-requested review of a recent session, use the
[session retrospective](session_retrospective.md). For test-level decisions,
provisional admission and CI ownership, use the [pytest workflow](pytest_workflow.md).
