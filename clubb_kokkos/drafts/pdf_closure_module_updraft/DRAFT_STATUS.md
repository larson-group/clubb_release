# Source draft status

`calc_w_up_in_cloud` is implemented and separately source-audited against
`src/CLUBB_core/pdf_closure_module.F90:3143–3340`. All 13 required arguments,
both components' three branches, all four output formulas and meaningful
computation comments are present. There are no pending lower physics contracts.

The supported computation remains **uncompiled and numerically unverified**.
Only isolated draft files were written; production sources, CMake, bridge
generators and tests were left untouched. See `INTEGRATION.md` for the merge
boundary, Native comment correction, arithmetic notes and the separate private
interpolation-cluster scope finding.
