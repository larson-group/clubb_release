---
name: clubb-doc-language
description: Choose clear wording and reader context for CLUBB documentation, comments, PRs, and commits being written or revised.
---

# Language for CLUBB documentation

Read the [shared documentation requirements](../SKILL.md) for source
inspection and verification. Directory placement does not load them implicitly.

Match the explanation to where someone encounters it:

| Text | Language and context |
| --- | --- |
| README | Start in very simple, ordinary language: what this is for and what someone can do with it. Introduce implementation terms only when needed. |
| Top-of-file description | Start simply, then introduce and define the file's terms, using a concrete example when helpful. |
| Procedure and code comments | Reuse terms introduced in the file description. Give enough local context to understand the purpose; leave calculation details to the code when appropriate. |
| PR | Begin with the concrete problem or unmet need, then explain how the change addresses it; give technical details and evidence below. |
| Commit | Use a simple, specific title and a short explanation. Include technical detail when it is important to explain the change or qualify its result. |

Explain what something means or does before describing how to use or compute
it. When order matters, make the starting point, direction, or sequence clear.

For example, introduce **truth profiles** as “truth profiles (reference
profiles, such as results from large-eddy simulations)” before using the term
on its own. Define an abbreviation where it is first needed.

Use a few extra familiar words instead of a shorter term that needs extra
technical background. Say what is
read, computed, changed, or returned, and name the program, routine, or data
involved. Prefer “reads the case settings” to “ingests configuration state.”
Avoid stacks of abstract nouns and vague words such as “handled” when a
specific action would explain the work.

Use established CLUBB names and scientific terms consistently. Before inventing
a name for part of the code, look for its existing name. If a new term is
necessary, define it in the relevant README or file description and connect it
to the actual code and familiar concepts. Do not give an ordinary operation a
new technical-sounding label.

Put the main point first and keep each sentence focused. Replace unclear “it,”
“this,” or “that” with the thing being discussed. Preserve scientific meaning,
conditions, units, and qualifications while simplifying. If someone cannot
follow the purpose without more context, supply the missing definition or
connection rather than assuming they followed the implementation work.
