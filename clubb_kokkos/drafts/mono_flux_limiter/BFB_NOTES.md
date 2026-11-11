# Arithmetic and source-definedness record

No explicit FMA, forced-rounding boundary, alternate math implementation or
new numerical guard was introduced. Native association, strict comparisons,
signed spacing and sequential recurrence/integral order remain. Source
squares expand as multiplication without changing unary-minus precedence.

`l_force_descending_solves` is an existing Native option, preserved at the
actual solver boundary. Its source comment explains grid-direction bit
differences; retaining it is not a newly verified C++ arithmetic workaround.
There is no numerical evidence for this draft yet.

The fast screen and sequential flux sweep remain separate. The screen's
logical reduction uses an integer count/nonzero check, which affects only
dispatch and performs no floating-point reduction. The spike integral's
array SUM is a direct sequential expansion in Native physical level order.
The Native sampled min-limit sentinel and success/fatal flip ordering remain.

No new source-definedness deviation was identified or introduced. Lower
solver outputs on fatal return and source assumptions about grid dimensions,
indices and valid density/variance inputs remain the original contracts.
Native square-root/exp/erf behavior is used without JAX clipping or safe
divisors. CPU whole-model comparison remains pending; GPU is a later
milestone under the accepted validation policy.
