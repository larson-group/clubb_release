"""Disabled interfaces from gfdl_activation.F90.

Initialization rejects this feature. These source-signature dummies must never silently simulate
an enabled feature; replace them when its driver is ported.
"""


# -----------------------------------------------------------------------------
def aer_act_clubb_quadrature_Gauss(
    gr, ngrdcol, pdf_params, p_in_Pa,  # In
    aeromass_clubb,                    # InOut
    temp_clubb_act,                    # In
    Ndrop_max,                         # Out
):
    """The main subroutine used for the GFDL droplet activation.

    Arguments:
        gr: Grid coordinates, interpolation weights and vertical metrics.
        ngrdcol: Number of grid columns.
        pdf_params: CLUBB PDF component parameters.
        p_in_Pa: Pressure [Pa].
        aeromass_clubb: Aerosol mass fields (source input/output).
        temp_clubb_act: Temperature supplied to activation [K].
        Ndrop_max: Maximum activated droplet concentration (source output).
    """

    # Description:
    # The main subroutine used for the GFDL droplet activation.
    # =======================================================================

    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")


# -----------------------------------------------------------------------------
def erff(x):
    """Retain the GFDL error-function entry point until activation is ported.

    Arguments:
        x: Dimensionless error-function input.
    """
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")


# -----------------------------------------------------------------------------
def Loading(droplets, droplets2):
    """Retain the disabled host lookup-table loading interface.

    Arguments:
        droplets: Source activation lookup table.
        droplets2: Source secondary activation lookup table.
    """
    # Description:
    # Loads the lookup tables for droplet activation into memory from flat data files.
    # ---------------------------------------------------------------------------
    #
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")


# -----------------------------------------------------------------------------
def get_unit():
    """Retain the source file-unit allocator as a disabled host interface.
    """
    # Description:
    # Used to determine a free unit with which to open a file.
    # ---------------------------------------------------------------------------
    #
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")
