# Port Fortran declarations to JAX docstrings

[port_fortran_declarations_to_jax_docstrings.py](port_fortran_declarations_to_jax_docstrings.py)
keeps the docstrings of the JAX ports in
[clubb_jax/src/CLUBB_core](../../clubb_jax/src/CLUBB_core) in step with the
Fortran routines in [src/CLUBB_core](../../src/CLUBB_core) they were ported
from. Each ported function's docstring ends with the routine's Input,
Input/Output and Output declarations, copied verbatim from the Fortran source,
after a description copied from the routine's Description and References
comments. The script checks the declarations, rewrites them after the Fortran
changes, and regenerates the description from the Fortran comments plus
recorded edits in [custom_ports_of_fortran_descriptions.jsonc](custom_ports_of_fortran_descriptions.jsonc).

## Quick start

Run from the repository root; only the Python standard library is needed.

```bash
# Report docstrings whose Fortran declaration section is missing or outdated
python3 utilities/port_fortran_declarations_to_jax_docstrings/port_fortran_declarations_to_jax_docstrings.py

# Rewrite those sections in two modules
python3 utilities/port_fortran_declarations_to_jax_docstrings/port_fortran_declarations_to_jax_docstrings.py -write_declaration -modules new_pdf,advance_windm_edsclrm_module

# Rewrite the description of one module from its Fortran comments and custom_ports_of_fortran_descriptions.jsonc
python3 utilities/port_fortran_declarations_to_jax_docstrings/port_fortran_declarations_to_jax_docstrings.py -write_description -modules new_pdf

# Export the matches, Fortran comments and sections as JSON
python3 utilities/port_fortran_declarations_to_jax_docstrings/port_fortran_declarations_to_jax_docstrings.py -manifest -output_file /tmp/docstring_manifest.json
```

## Modes

- **`-check`** (default): lists every matched function whose section is
  missing or differs from the Fortran source, and exits with 1 if there are
  any. Notes about uncertain matches are printed first. The "Fortran
  declaration match" stage of the [clubb_jax Jenkins job](../../jenkins_tests/clubb_jax/Jenkinsfile)
  runs it; the script itself is tested in
  [test_port_fortran_declarations_to_jax_docstrings.py](../auto_llm_generated_tests/pytests/unit_tests/test_port_fortran_declarations_to_jax_docstrings.py).
- **`-write_declaration`**: inserts or replaces the section in every matched
  function. Text above the section is kept. Before saving a file, the script
  confirms that its code, with docstrings removed, is unchanged.
- **`-write_description`**: rewrites the text above the section from the
  Fortran comments and the edits in `-edits` (default
  [custom_ports_of_fortran_descriptions.jsonc](custom_ports_of_fortran_descriptions.jsonc)),
  then writes the section as `-write_declaration` does. Functions without a Fortran description keep
  their text.
- **`-manifest`**: writes one JSON record per matched function: the JAX file,
  function and line, the Fortran file, routine and line, the Description and
  References text without `!`, the section, the current docstring, and
  matching problems. Without `-output_file` the JSON goes to standard output.

`-modules` limits any mode to the given JAX module names (file stems).
`-fortran_dir` and `-jax_dir` point at other source directories, e.g. a copy
used to preview `-write_declaration`.

## Docstring layout

```python
def calc_setter_var_params(xm, xp2, Skx, sgn_wpxp, F_x, zeta_x):
    """Calculates the PDF component means, the PDF component standard deviations,
    and the mixture fraction for the variable that sets the PDF.

    References:
    Griffin and Larson (2018)

    Returns:
        (mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, mixt_frac, coef_sigma_x_1_sqd,
          coef_sigma_x_2_sqd)

    Fortran declarations (calc_setter_var_params in src/CLUBB_core/new_pdf.F90):

        !--------------------------- Input Variables ---------------------------
        integer, intent(in) :: &
          nz
        ...
    """
```

`-write_declaration` owns everything from the `Fortran declarations (` line
to the end of the docstring. `-write_description` owns the text above it. That
text is the full Fortran comment, however long;
`custom_ports_of_fortran_descriptions.jsonc` holds only the exceptions, for
the functions whose Fortran text does not apply to the JAX port. A function without an entry gets the Fortran text unchanged. The
text is built from:

- The Fortran description and references, with `!` removed.
- A `Returns:` block read from the function's `return` statements: the
  returned names, the keys of a returned dict, or the Fortran result name for
  a function that returns one expression. It is left out when the returns
  differ.
- The edits for the function, keyed by module and function name. The file is
  JSON that may contain `//` and `/* */` comments.

  | Key | Effect |
  | --- | --- |
  | `replace` | `[old, new]` pairs replaced in the Fortran text, for Fortran-only names such as `gr%nz` or LAPACK |
  | `drop` | lines containing these strings are removed, e.g. notes about Fortran optimization |
  | `append` | a paragraph added after the references, for real JAX differences such as unported branches |

  Edits that no longer match the Fortran text, and unknown keys, stop the run,
  so stale edits are noticed when the Fortran comments change. When a Fortran
  comment does not fit this layout, make the Fortran comment more standard
  instead of replacing its text here.

## Matching and limitations

- A JAX module is paired with the Fortran file of the same name. Functions
  match routines of the same name, ignoring case. A generic interface
  resolves to the specific routine whose arguments best match the JAX
  arguments; `SPECIFIC_ROUTINE` in the script settles remaining ties.
- Private helpers (`_name`) and functions with no Fortran counterpart are
  skipped. Routines without arguments get no section.
- The copied block is labelled by the Fortran section dividers
  (`!---- Input Variables ----`, see the
  [Fortran formatting guide](../../.agents/skills/clubb-fortran-format/SKILL.md)).
  A routine without dividers is reported, and its declarations are copied
  without headings.
