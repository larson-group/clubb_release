"""Disabled COAMPS interfaces from coamps_microphys_driver_module.F90.

The standalone initialization rejects this feature. These explicit dummy interfaces preserve
routine names and argument order and fail if called; they do not produce a simulated enabled
result. Fortran out arguments remain in these dummy signatures because execution is never
permitted.
"""


# -----------------------------------------------------------------------------
def coamps_microphys_driver(
    gr, runtype, timea_in, deltf_in,  # In
    rtm, wm_zm, p_in_Pa, exner, rho,  # In
    thlm, rim, rrm, rgm, rsm,         # In
    rcm, Ncm, Nrm, Nim,               # In
    saturation_formula,               # In
    stats,                            # InOut
    icol,                             # In
    Nccnm, cond,                      # InOut
    Vrs, Vri, Vrr, VNr, Vrg,          # Out
    ritend, rrtend, rgtend,           # Out
    rstend, nrmtend,                  # Out
    ncmtend, nimtend,                 # Out
    rvm_mc, rcm_mc, thlm_mc,          # Out
):
    """Retain the COAMPS interface as an explicitly disabled entry point.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        runtype: Benchmark case name.
        timea_in: Current model time [s]
        deltf_in: Timestep (i.e. dt_main in CLUBB) [s]
        rtm: Total water mixing ratio [kg/kg]
        wm_zm: Vertical wind [m/s]
        p_in_Pa: Pressure [Pa]
        exner: Mean exner function [-]
        rho: Mean density [kg/m^3]
        thlm: Liquid potential temperature [K]
        rim: Ice water mixing ratio [kg/kg]
        rrm: Rain water mixing ratio [kg/kg]
        rgm: Graupel water mixing ratio [kg/kg]
        rsm: Snow water mixing ratio [kg/kg]
        rcm: Cloud water mixing ratio [kg/kg]
        Ncm: Number of cloud droplets [count/kg]
        Nrm: Number of rain drops [count/kg]
        Nim: Number of ice crystals [count/kg]
        saturation_formula: Choice of liquid/ice saturation formula.
        stats: Immutable statistics state; return its updated value.
        icol: Source column index retained in the disabled interface.
        Nccnm: Number of cloud nuclei [count/kg]
        cond: condensation/evaporation of liquid water
        Vrs: Snow mixing ratio fall speed [m/s]
        Vri: Pristine ice mixing ratio fall speed [m/s]
        Vrr: Rain mixing ratio fall speed [m/s]
        VNr: Rain conc. fall speed [m/s]
        Vrg: Graupel fall speed [m/s]
        ritend: d(ri)/dt [kg/kg/s]
        rrtend: d(rr)/dt [kg/kg/s]
        rgtend: d(rg)/dt [kg/kg/s]
        rstend: d(rs)/dt [kg/kg/s]
        nrmtend: d(Nrm)/dt [count/kg/s]
        ncmtend: d(Ncm)/dt [count/kg/s]
        nimtend: d(Nim)/dt [count/kg/s]
        rvm_mc: d(rv)/dt [kg/kg/s]
        rcm_mc: d(rc)/dt [kg/kg/s]
        thlm_mc: d(thlm)/dt [K/s]
    """

    raise NotImplementedError("COAMPS is disabled in the JAX standalone")
