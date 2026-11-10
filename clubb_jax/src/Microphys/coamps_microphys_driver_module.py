"""Disabled COAMPS interfaces from coamps_microphys_driver_module.F90.

The standalone initialization rejects this feature. These explicit dummy interfaces preserve
routine names and argument order and fail if called; they do not produce a simulated enabled
result. Fortran out arguments remain in these dummy signatures because execution is never
permitted.
"""


# -----------------------------------------------------------------------------
def coamps_microphys_driver(
    gr, ngrdcol,
    runtype, timea_in, deltf_in,
    rtm, wm_zm, p_in_Pa, exner, rho,
    thlm, T_in_K,
    rim, rrm,
    rgm, rsm,
    rcm, Ncm,
    Nrm, Nim,
    saturation_formula,
    stats, Nccnm,
    Vrs, Vri,
    Vrr, VNr,
    Vrg,
    ritend, rrtend,
    rgtend, rstend,
    nrmtend, ncmtend, nimtend,
    rvm_mc, rcm_mc, thlm_mc,
):
    """Reject a COAMPS request when the imported scheme was not compiled.

    Fortran declarations (coamps_microphys_driver in src/Microphys/coamps_microphys_driver_module.F90):

        !--------------------------- Input Variables ---------------------------
        type(grid), intent(in) :: gr ! CLUBB vertical grids

        integer, intent(in) :: ngrdcol ! Number of grid columns [-]

        character(len=*), intent(in) :: runtype ! Case name

        real( kind = time_precision ), intent(in) :: timea_in ! Current model time [s]
        real( kind = core_rknd ), intent(in) :: deltf_in ! Main timestep [s]

        real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
          rtm ! Total water mixing ratio [kg/kg]

        real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
          wm_zm ! Vertical wind on momentum levels [m/s]

        real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
          p_in_Pa, & ! Pressure [Pa]
          exner,   & ! Exner function [-]
          rho,     & ! Density [kg/m^3]
          thlm,    & ! Liquid potential temperature [K]
          T_in_K,  & ! Absolute temperature computed for the full batch [K]
          rim,     & ! Ice water mixing ratio [kg/kg]
          rrm,     & ! Rain water mixing ratio [kg/kg]
          rgm,     & ! Graupel water mixing ratio [kg/kg]
          rsm,     & ! Snow water mixing ratio [kg/kg]
          rcm,     & ! Cloud water mixing ratio [kg/kg]
          Ncm,     & ! Number of cloud droplets [count/kg]
          Nrm,     & ! Number of rain drops [count/kg]
          Nim        ! Number of ice crystals [count/kg]

        integer, intent(in) :: saturation_formula ! Saturation formula selector [-]

        !------------------------ Input/Output Variables ------------------------
        type(stats_type), intent(inout) :: stats ! Statistics buffers and sample controls

        real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
          Nccnm ! Number of cloud nuclei [count/kg]

        !--------------------------- Output Variables ---------------------------
        ! All sedimentation velocities are on thermodynamic levels.
        real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
          Vrs,     & ! Snow mixing ratio fall speed [m/s]
          Vri,     & ! Pristine ice mixing ratio fall speed [m/s]
          Vrr,     & ! Rain mixing ratio fall speed [m/s]
          VNr,     & ! Rain concentration fall speed [m/s]
          Vrg,     & ! Graupel mixing ratio fall speed [m/s]
          ritend,  & ! Ice water tendency [kg/kg/s]
          rrtend,  & ! Rain water tendency [kg/kg/s]
          rgtend,  & ! Graupel water tendency [kg/kg/s]
          rstend,  & ! Snow water tendency [kg/kg/s]
          nrmtend, & ! Rain number tendency [count/kg/s]
          ncmtend, & ! Cloud droplet number tendency [count/kg/s]
          nimtend, & ! Ice number tendency [count/kg/s]
          rvm_mc,  & ! Vapor water tendency [kg/kg/s]
          rcm_mc,  & ! Cloud water tendency [kg/kg/s]
          thlm_mc    ! Liquid potential temperature tendency [K/s]
    """

    raise NotImplementedError("COAMPS is disabled in the JAX standalone")
