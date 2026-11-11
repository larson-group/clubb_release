# LY93 PDF draft

[`LY93_pdf.hpp`](LY93_pdf.hpp) implements all three public routines from
[`LY93_pdf.F90`](../../../src/CLUBB_core/LY93_pdf.F90), in native order:
`LY93_driver` (23 arguments), `calc_mixt_frac_LY93` (2 arguments and an owned
returned profile), and `calc_params_LY93` (9 arguments). There are no optional
arguments, omitted options, native statistics calls or error objects.

The header is isolated from the live build. After source review, promote it to
the existing core directory, include it in `pdf_closure_module.cpp`, replace the
existing per-column `f2k::LY93_driver` call with the same C++ call, and remove the
obsolete whole-driver bridge through its existing generator owner. The public
profile signatures already accept the current strided column subviews. No new
lower bridge, model constant, polynomial-root dependency or CMake source file
is required; all routines are templates and existing constant owners suffice.

The native `nz` dimension remains a rank-one profile. Scratch and returned
mixtures use LayoutLeft in the input View's memory space, and each policy uses
that input's execution space. Caller fields retain their existing column-first
layout. All fields must be conformable and accessible from the chosen execution
space, as the existing callers provide. Inputs and output buffers retain their
native roles; no global storage or host staging is introduced.

The mixture bisection stays serial within each independent level. It retains
the `Sk_max > 0.84` branch, `[0.5,1]` starting interval, absolute residual
tolerance `1e-4`, bracket update order and uncapped native termination. The
otherwise branch returns `0.75`. No JAX scan limit or alternate convergence
criterion is imported. Extreme or nonfinite inputs that do not converge under
the native rule are not repaired by this port.

The parameter kernel keeps source formula and field order while lowering the
level-local `sgn_Skx` and `B_x` array intermediates to distinct work-item scalars.
It uses the native direct square root and real one-third power, including the
`Skx >= 0` sign choice. JAX's gradient-safe square root/cube-root helpers are
language-specific adaptations and are not needed here. Ordinary integer squares
use products; the sixth power remains a direct Kokkos power expression. There
is no explicit FMA, forced rounding or mismatch-driven arithmetic adjustment.
The copied component-variance unit descriptions are corrected to `(kg/kg)^2`
for total water and `K^2` for temperature; no native science file is changed.

[SOURCE_OUTLINE.json](SOURCE_OUTLINE.json) records the full native/JAX read before
implementation. The separate final [source audit](SOURCE_AUDIT.json) and
[expression audit](expression_audit.json) cover all signatures, 16 ordered
assignment trees, all three field-specific parameter calls, bisection and sign
branches, ownership, comments and formatting. They establish source structure.
Compilation, numerical agreement, alternate-case coverage and performance are
pending. Root should use the maintained CPU Fortran/CPU Kokkos comparison at zero
tolerance after integration and record any deliberate BFB adjustment in the
existing [adaptation log](../../BFB_ADAPTATIONS.md).
