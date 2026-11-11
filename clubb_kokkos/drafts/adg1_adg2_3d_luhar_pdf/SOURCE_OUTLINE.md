# Full ADG module draft source outline

Authority: all1537 lines of `src/CLUBB_core/adg1_adg2_3d_luhar_pdf.F90`, read before expansion. Actual JAX is scope context, with native equations/branches/comments retained. The compact parent+responder draft remains recoverable in `incremental_adg1/` and `/tmp/clubb_adg1_next/`.

Native routine order and coverage:

1. ADG1_pdf_driver32–240:58 required arguments; w closure, debug0 all-fatal gate, rt/thl/u/v responders (native thl_tol for winds), optional scalar loop. Preserve existing compact body except calling new C++ w closure.
2. ADG2_pdf_driver243–412:41 required arguments; per-column calc_Luhar_params→close_Luhar_pdf, overwrite sigma_sqd_w with min(1/(1+m²),0.99), rt/thl responders, scalar flag/species loop. No additional fatal gate.
3. Luhar_3D_pdf_driver415–625:25 required profile arguments; source level loop; w wins magnitude ties, thl second with strict>w and>=rt, rt fallback. Each branch setter calc→close, then two responder backsolve→close sequences. Preserve all field-specific covariance/tolerance mappings and source order.
4. ADG1_w_closure628–782:16 required arguments; debug2 variance ANY all-fatal return; abs(Skw)<=1e-5 half-mixture branch; native mixture clamp; normalized/physical means; equal component variances. Golaz et al.2002a and Larson/Golaz2005 references.
5. calc_Luhar_params786–950:8 required profile arguments; xp2>tol varying branch, fractional power1/3 with m floor0.05; m²/M/covariance sign/mixture. Nonvarying returns mixture0.5,m=M=0. Larson et al.2002 Appendix e derivation.
6. close_Luhar_pdf953–1069:15 required profile arguments; varying/nonvarying outputs in source order. Definedness adaptation: initialize sign from wpxp for both branches, matching JAX; native initializes it only for varying cells but reads it in nonvarying normalized means. No native BFB claim for those undefined normalized outputs.
7. ADG1_ADG2_responder_params1072–1233:19 required arguments; compact GPU implementation retained, debug2 all-column tolerance failure before output writes, source formulas/clamp, inactive old branch/comment/PR1276.
8. backsolve_Luhar_params1236–1437:6 scalar arguments; compile-time cubic true. Near-half m=.05/M; zero-skewness m0/M=HUGE; otherwise alpha/coefficient setup→max_cubic_root→sqrt(max(root,.05²)). Preserve inactive quadratic approximation branch, all clamps/discriminants and finalM.
9. max_cubic_root1440–1537:4 singleton arrays→scalar; thresholds .001*coefficient maxima; cubic public solve with conditional3-real max, quadratic fallback, linear fallback. Preserve native c=0 division (no JAX safe denominator).

All argument groups/intents come from current source declarations. No optionals. Public routines are first six; private responders/backsolve/root retain native privacy. No stats or persistent state in any routine.

Execution boundary: ADG1 and responder/w kernels remain Kokkos execution-space loops. Luhar profile orchestration/math is host-side C++ while native complex root solvers remain unported; public profile adapters stage caller Views explicitly. This is a draft performance boundary to revisit with root solvers, not a whole-native wrapper.

Missing lower bridge contracts: public calc_roots::cubic_solve(nz,a_coef,b_coef,c_coef,d_coef) returns complex(nz,3), quadratic_solve(nz,a_coef,b_coef,c_coef) returns complex(nz,2). Native metadata debug service already exists. No production parser/generator is modified; local draft bridge declarations document this missing integration.

BFB baseline: no runtime evidence for the expanded draft. No speculative FMA. Ordinary power/sqrt and source expressions first; native MIN/MAX operand order retained. Sign initialization is source-definedness repair from JAX, distinct from a BFB workaround. JAX's clipped roots, infinity replacement for HUGE, c=0 guard and gradient-safe masked arithmetic are not copied. Both native cubic and inactive quadratic backsolve branches will be implemented even though JAX implements only the active cubic choice.
