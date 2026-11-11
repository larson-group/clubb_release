# PDF moment draft

This contains four public moment routines from
[pdf_closure_module.F90](../../../src/CLUBB_core/pdf_closure_module.F90), in
native order: `calc_wp4_pdf`, `calc_wp2xp_pdf`, `calc_wpxp2_pdf` and
`calc_wpxpyp_pdf`. Full source review and cleanup are recorded in
[source_audit.json](source_audit.json). The four routines are promoted and pass four CPU cases with multicolumn and all stats. A following correlation-fusion cleanup also passes four multicolumn and two cloudy all-stats cases. New GPU checks are pending.

The draft replaces each native column/level loop with one Kokkos kernel and
keeps the complete formulas, arguments, comments and physical units. Integer
squares are products; fourth powers are squared squares. No new floor, helper,
state, FMA choice or other BFB workaround was added.

For integration, merge the declarations into the existing public header before
`pdf_closure_driver`. Keep native definition order in the existing C++ file:
`calc_wp4_pdf` follows private `transform_pdf_chi_eta_component`; private
`calc_wp2xp2_pdf` follows it; the remaining three public functions follow that
helper and precede the private liquid-cloud helper. Close/reopen the anonymous
namespace at those public/private boundaries. No new compilation unit is needed.

Replace the reached calls without changing their arguments, remove the four
whole-routine entries from the existing bridge generator, regenerate, build and
run the maintained CPU comparison before promoting the next part. Standard/all
stats are needed occasionally to exercise the conditionally sampled moments;
use short multicolumn stats for the regular loop. GPU validation remains a
separate milestone check.
