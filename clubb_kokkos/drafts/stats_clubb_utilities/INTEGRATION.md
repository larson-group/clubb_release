# Statistics accumulation source review

All three Native routines are promoted to
[stats_clubb_utilities.cpp](../../src/CLUBB_core/stats_clubb_utilities.cpp) and its
[public interface](../../src/CLUBB_core/stats_clubb_utilities.hpp). Their argument
counts are 66, 7 and 17. This directory keeps the pre-port outline, source,
cleanup, compile and independent review records; no implementation here builds.
`PROMOTION.json` maps their original draft hashes to the permanent owners.

The Native routine order, 102 ordered field updates, 15 registry queries and
29 scientific lower calls were reviewed separately from numerical acceptance.
Host accumulation and gated staging use the permanent
[statistics types](../../src/CLUBB_core/stats_netcdf.hpp). Grid metadata remains
one-based; explicit stats column/level arguments are zero-based. No new model
globals, Native compute fallback or arithmetic workaround was introduced.

See the [driver milestones](../../DRIVER_CONVERSION_PLAN.md) for integrated
comparison results and their tested scope. The original draft snapshots and
compile commands are local evidence under `build/validation/stats_cpp/`.
