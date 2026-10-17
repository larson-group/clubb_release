---
name: clubb-file-descriptions
description: Write or revise a CLUBB file's opening description or module docstring.
---

# Top-of-file descriptions

Read the [shared documentation requirements](../SKILL.md) for source
inspection and verification. Directory placement does not load them implicitly.

Use the [language guide](../language/SKILL.md). Treat the opening
description as a small README for this file. Replace a weak header with a
complete description using the structure below; a one-line original does not
limit the rewrite to one line.

Use opening `!` comments for Fortran and a module docstring for Python. Keep
shebangs, encoding declarations, and executable code outside the replacement.

Read the opening headers of `src/clubb_loss_driver.F90` and
`src/CLUBB_core/stats_netcdf.F90` in the CLUBB checkout. Follow their useful
structure and level of explanation; do not copy their file-specific claims.

Write in this order:

1. **Purpose:** a simple one-line summary, then a short explanation of what the
   file provides, when to use it, and the important terms the reader will encounter.
2. **Usage:** list every public routine with a brief description of its purpose
   and when to call it. Group routines by use or call order, as in the reference
   headers: initialize/finalize, each timestep, each variable, and optional
   operations where those groups apply. Keep the full public list, including
   less common entry points.
3. **Understanding:** explain the details needed to use the file correctly.
   Include applicable input formats/examples, output layouts and units, call
   ordering, data kept between calls, repeated runs, and important restrictions.
   For checking scripts, explain important defaults and tolerances and what
   passing establishes. Use descriptive subsections for the file's actual behavior.
4. **References:** retain existing citations and add relevant source documents
   or scientific references when the implementation warrants them.

Prefer at least one concise usage example when it helps explain the file's
overall purpose. Existing callers can supply examples; describe the file's
interface independently of its current callers. For files collecting many
independent helpers, useful examples belong with individual procedure
descriptions. Examples are optional.

Determine public routines from exports and actual callers. Do not list imported
names or private helpers as the file's public routines. Confirm descriptions against
the implementation, including non-default behavior being described.

Preserve accurate detail and citations when reorganizing a good header. Change documentation only during a documentation
request; leave code formatting and behavior alone.
