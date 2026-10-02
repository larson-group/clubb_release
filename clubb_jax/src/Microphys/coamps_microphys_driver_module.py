"""Disabled COAMPS interfaces from coamps_microphys_driver_module.F90.

The standalone initialization rejects this feature. These explicit dummy
interfaces preserve routine names and argument order and fail if called;
they do not produce a simulated enabled result. Fortran out arguments remain
in these dummy signatures because execution is never permitted.
"""


def coamps_microphys_driver(gr, runtype, timea_in, deltf_in, rtm, wm_zm, p_in_Pa, exner, rho, thlm, rim, rrm, rgm, rsm, rcm, Ncm, Nrm, Nim, saturation_formula, stats, icol, Nccnm, cond, Vrs, Vri, Vrr, VNr, Vrg, ritend, rrtend, rgtend, rstend, nrmtend, ncmtend, nimtend, rvm_mc, rcm_mc, thlm_mc):
    raise NotImplementedError("COAMPS is disabled in the JAX standalone")
