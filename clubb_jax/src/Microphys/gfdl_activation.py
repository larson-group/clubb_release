"""Disabled interfaces from gfdl_activation.F90.

Initialization rejects this feature. These source-signature dummies must never
silently simulate an enabled feature; replace them when its driver is ported.
"""


def aer_act_clubb_quadrature_Gauss(gr, ngrdcol, pdf_params, p_in_Pa, aeromass_clubb,
    temp_clubb_act, Ndrop_max):
    # Description:
    # The main subroutine used for the GFDL droplet activation.
    #=======================================================================
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")


def erff(x):
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")


def Loading(droplets, droplets2):
    # Description:
    # Loads the lookup tables for droplet activation into memory from flat data files.
    #---------------------------------------------------------------------------
    #
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")


def get_unit():
    # Description:
    # Used to determine a free unit with which to open a file.
    #---------------------------------------------------------------------------
    #
    raise NotImplementedError("GFDL aerosol activation is disabled in the JAX standalone")
