"""JAX port of the WRF/SAM Morrison 2-moment microphysics (module_mp_graupel.F90).

Started Iter190 with the foundational, independently-verifiable special functions
(POLYSVP saturation vapor pressure, GAMMA, DERF1) — the KK-port playbook: port the
self-contained transcendental/utility functions first (validatable vs scipy / known
physical values), then the process rates (validatable via a Morrison case-stats oracle).
"""
