---
name: clubb-code-comments
description: Write or revise descriptions of CLUBB procedures and comments within code, including driver code and tests.
---

# Procedure and code comments

Read the [shared documentation requirements](../SKILL.md) for source
inspection and verification. Directory placement does not load them implicitly.

Use the [language guide](../language/SKILL.md).

## Procedure descriptions

Every procedure in the requested scope needs a description directly below its
definition. Explain its purpose first, then its inputs and outputs. Include
units and array meanings when needed to understand those inputs or outputs.
Outputs include changes to existing data and files written, not just return
values. Include a short usage example when helpful, especially for independent
helpers. For a test, describe the setup and what its assertions check.
Include references where appropriate for the method, equation, or algorithm;
preserve accurate existing citations.

In Python, put the docstring as the first statement inside the definition,
after the complete signature. In Fortran, put `Description:` immediately below
the subroutine/function statement, before `use` statements. Do not place a
standalone string before a Python definition.

Fortran arguments and local variables also need brief declaration comments
explaining their purpose, with units where meaningful. Keep declaration layout
unchanged.

For a file-wide rewrite, review every procedure, including helpers. For a
narrow request, stay within the requested procedures or comment region.

## Comments within code

Use concise comments to explain the purpose of a small section of code, a
choice that needs explaining, or why operations occur in that order. Avoid
comments that merely repeat an assignment or function name. When an equation
helps explain a scientific calculation, introduce what it computes and define
its symbols.

## Driver and other high-level code

Code that coordinates calls and stages needs more detailed commentary. Give
each meaningful stage a concise explanation of what happens, why, and what is
passed to the next stage. Explain important branches, repeated runs, and
cleanup when present.

Read just the comments in sequence: a reader should be able to follow the
high-level workflow without tracing every call. They may still need the code
to understand individual calculations. Add missing connections between stages
rather than narrating every statement.

Keep documentation edits to comments/docstrings; preserve executable code,
its spacing, and directive comments. Check descriptions against what the
implementation and assertions actually do.
