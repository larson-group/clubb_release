---
name: clubb-readmes
description: Write or revise CLUBB READMEs, including usage documentation for a substantial new feature.
---

# READMEs

Use the [language guide](../clubb-doc-language/SKILL.md).

Write READMEs in Markdown as `README.md`. Link named scripts, files,
directories, routines, and other codebase items at least on their first
mention in prose, where practical. Use paths relative to the README; a
routine can link to its defining file.

Create a README when the user requests one or a substantial new feature or
source area clearly needs its own explanation. Update the existing README
when it already covers the work. Do not create one for every script or helper.

Start with a concise explanation of what this README covers and what it helps
someone do. Then choose the structure for that area:

- **Broad script directories, such as `run_scripts/` or `tests/`:** give a few
  sentences of overview, then describe each user entry point with its uses
  and examples.
- **A coherent feature, such as the tuner or Dash:** explain its purpose and
  main use cases, give a usable starting example, then explain the feature in
  more detail.
- **Intermediate code, such as namelist utilities or code in `input/`:**
  explain which programs use it, what they supply, what they receive, and how
  they call or include it. Show a calling example where it helps.

Within the README's stated scope, list **every intended user-run script,
command, or executable**. Describe its purpose, required inputs, and outputs.
If it accepts optional arguments, give at least two concrete invocations: a
basic use and a use with optional arguments. Add examples for distinct use
cases and explain the options they use. Do not document only a favorite subset
of the runnable tools.

Check the actual entry points and argument parsers before describing them.
Show where commands are run and any required setup. Keep commands and example
values real; distinguish examples checked by inspection from commands actually
run. Explain purpose and use before internal implementation details. Keep detailed
test tolerances and internal setup explanations in the relevant file header;
the README should explain what the checks are for and how to run them.
