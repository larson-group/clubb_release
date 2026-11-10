"""Interpret the columns of an input file read by input_reader.py.

JAX port of input_interpret.F90. A case's sounding file can describe its
vertical coordinate, temperature and vertical velocity in more than one way.
Each routine here finds which form a file uses and returns that profile.

Usage:
    read_z_profile: return the altitudes of the sounding levels, converting
        them from pressure when the file gives pressure.
    read_theta_profile: return the temperature profile and its type.
    read_subs_profile: return the vertical velocity profile and its type.

    Pass the list from input_reader.read_one_dim_file, after
    input_reader.fill_blanks_one_dim_vars has filled its blank values.

Porting deviations:
    * Fortran returns its outputs through intent(out) arguments. These routines
      return them instead, in the same order.
    * Fortran leaves its outputs undefined when a file holds more than one form
      of the same quantity. These routines raise ValueError, because Python has
      no undefined value to return.

References:
    None.
"""

import jax.numpy as jnp
import numpy as np

from clubb_jax.src.CLUBB_core.calc_pressure import calculate_thvm
from clubb_jax.src.CLUBB_core.constants_clubb import (
    kappa,
    p0,
    Cp,
    Lv,
    ep2,
    one,
    zero_threshold,
)
from clubb_jax.src.CLUBB_core.saturation import (
    sat_mixrat_liq,
    rcm_sat_adj,
)
from clubb_jax.src.Input_fields.hydrostatic_module import inverse_hydrostatic
from clubb_jax.src.Input_fields.input_names import (
    z_name,
    pressure_name,
    rt_name,
    temperature_name,
    thetal_name,
    theta_name,
    wm_name,
    omega_name,
)
from clubb_jax.src.Input_fields.input_reader import read_x_profile


# -----------------------------------------------------------------------------
def read_z_profile(nvar, nsize, retVars, p_sfc, zm_init,
                   saturation_formula):
    """Searches for the variable specified by either 'z[m]' or 'Press[Pa]' in the
    collection of retVars.  If the subroutine finds the variable indicated by
    'z[m]', then it returns it.  If the subroutine finds 'Press[Pa]' then it
    converts it to values of altitude in meters.  If it does not find either
    or finds both, it raises ValueError.

    Args:
        nvar: number of elements in retVars
        nsize: number of values in each element of retVars
        retVars: collection being searched, from read_one_dim_file
        p_sfc: pressure at the surface [Pa]
        zm_init: height at zm(1) [m]
        saturation_formula: integer that stores the saturation formula to be used

    Returns:
        z: height sounding profile, shape (nsize,) [m]
        p_in_Pa: pressure sounding profile, shape (nsize,) [Pa]. Holds the fill
            value -999.9 when the file gives altitude.
        altitude_type: z_name or pressure_name, the column the file provides
    """
    names = [retVar.name for retVar in retVars]

    if z_name in names and pressure_name in names:
        raise ValueError("Sounding holds both " + z_name + " and " + pressure_name)

    if z_name in names:

        # The input sounding is given in terms of altitude.
        altitude_type = z_name

        # Obtain the value of altitude at each sounding level.
        z = read_x_profile(nvar, nsize, altitude_type, retVars)

        # Set the pressure at the sounding levels to the "fill value".
        p_in_Pa = np.full(nsize, -999.9)

    elif pressure_name in names:

        # The input sounding is given in terms of pressure.
        altitude_type = pressure_name

        # Obtain the value of total pressure at each pressure sounding level.
        p_in_Pa = read_x_profile(nvar, nsize, altitude_type, retVars)

        nlevels = retVars[0].values.size

        # Obtain the value of temperature_type (theta, theta_l, or temperature) at
        # each pressure sounding level.
        temperature_type, theta = read_theta_profile(nvar, nsize, retVars)

        # Obtain the value of total water mixing ratio at each pressure sounding
        # level.
        rtm = read_x_profile(nvar, nsize, rt_name, retVars)

        # Calculate exner from pressure.
        exner = (p_in_Pa / p0)**kappa

        if temperature_type == temperature_name:

            # The variable "theta" actually contains temperature (in Kelvin) at
            # this point.

            # Determine initial cloud water mixing ratio at sounding levels
            # based on temperature and rtm.  Again, "theta(k)" is actually
            # "temperature(k)" at this point in the code.
            rcm = np.maximum(
                rtm - np.asarray(sat_mixrat_liq(p_in_Pa, theta, saturation_formula)),
                zero_threshold,
            )

            # Convert temperature to potential temperature, theta.
            theta = theta / exner

            # Calculate theta_l from theta and cloud water mixing ratio, such
            # that:  theta_l = theta - [Lv/(Cp*exner)]*rcm.
            thlm = theta - Lv / (Cp * exner) * rcm

        elif temperature_type == theta_name:

            # This variable "theta" does indeed contain potential temperature.

            # Determine initial cloud water mixing ratio at sounding levels
            # based on potential temperature, exner, and rtm.
            rcm = np.maximum(
                rtm - np.asarray(sat_mixrat_liq(p_in_Pa, theta * exner, saturation_formula)),
                zero_threshold,
            )

            # Calculate theta_l from theta and cloud water mixing ratio, such
            # that:  theta_l = theta - [Lv/(Cp*exner)]*rcm.
            thlm = theta - Lv / (Cp * exner) * rcm

        elif temperature_type == thetal_name:

            # The variable "theta" actually contains liquid water potential
            # temperature, theta_l, at this point.
            thlm = theta

            # Determine initial cloud water mixing ratio.  If the profile is
            # unsaturated, then theta = theta_l.  If the profile is saturated at
            # any level, then cloud water mixing ratio must be determined using
            # an iterative method involving theta_l, total water mixing ratio,
            # pressure, and exner.
            rcm = np.asarray(rcm_sat_adj(thlm, rtm, p_in_Pa, exner, saturation_formula))

            # Calculate theta from theta_l and cloud water mixing ratio, such
            # that:  theta = theta_l + [Lv/(Cp*exner)]*rcm.
            theta = thlm + Lv / (Cp * exner) * rcm

        else:

            raise ValueError("Invalid temperature_type: " + temperature_type)

        # Now, the variable "theta" contains potential temperature, theta.
        # Compute initial theta_v based on theta, total water mixing ratio, and
        # cloud water mixing ratio.
        #
        # theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
        #
        # which can be rearranged as:
        #
        # theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
        #
        # The exact form of the equation for theta_v (including water vapor and
        # cloud water) is:
        #
        # theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
        #
        # which can be rearranged as:
        #
        # theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
        #                         * { r_v / ( 1 + r_v + r_c ) }
        #                     - { r_c / ( 1 + r_v + r_c ) } ].
        #
        # This version is written with r_t = r_v + r_c, such that:
        #
        # theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
        #                         * { ( r_t - r_c ) / ( 1 + r_t ) }
        #                     - { r_c / ( 1 + r_t ) } ].
        #
        # To use theta_l instead of theta, simply substitute the following for
        # theta in the above expression:  theta_l + {L_v/(C_p*exner)} * r_c.
        #
        # The CLUBB code uses a linearized version of the above equation (in
        # order to calculate thv' terms -- such as w'thv', etc.) for theta_v
        # throughout the model code, such that:
        #
        # theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
        #                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
        #
        # where thv_ds is used as a reference value to approximate theta_l.
        thvm = np.asarray(calculate_thvm(
            nlevels, 1,
            jnp.asarray(thlm), jnp.asarray(rtm), jnp.asarray(rcm), jnp.asarray(exner),
            jnp.asarray(theta * (one + ep2 * (rtm - rcm))**kappa),
        ))

        # Find the altitudes, z, of the pressure sounding levels.
        z = np.asarray(inverse_hydrostatic(float(p_sfc), float(zm_init), thvm, exner))

    else:

        raise ValueError("Could not read sounding vertical-coordinate variable")

    return z, p_in_Pa, altitude_type


# -----------------------------------------------------------------------------
def read_theta_profile(nvar, nsize, retVars):
    """Searches for the variable specified by either 'thetal[K]' or 'theta[K]' in the
    collection of retVars. If the function finds the variable then it returns
    it. If it does not it raises ValueError.

    Args:
        nvar: number of elements in retVars
        nsize: number of values in each element of retVars
        retVars: collection being searched, from read_one_dim_file

    Returns:
        temperature_type: theta_name, thetal_name or temperature_name
        theta: temperature profile of that type, shape (nsize,) [K]
    """
    names = [retVar.name for retVar in retVars]

    if sum(name in names for name in (theta_name, thetal_name, temperature_name)) > 1:
        raise ValueError("Sounding holds more than one theta compatible variable")

    if theta_name in names:
        temperature_type = theta_name
    elif thetal_name in names:
        temperature_type = thetal_name
    elif temperature_name in names:
        temperature_type = temperature_name
    else:
        raise ValueError("Could not read theta compatable variable")
    theta = read_x_profile(nvar, nsize, temperature_type, retVars)

    return temperature_type, theta


# -----------------------------------------------------------------------------
def read_subs_profile(nvar, nsize, retVars):
    """Searches for the variable specified by either 'w[m\\s]' or 'omega[Pa\\s]' in the
    collection of retVars. If the function finds the variable then it returns
    it. If it does not it raises ValueError.

    Args:
        nvar: number of elements in retVars
        nsize: number of values in each element of retVars
        retVars: collection being searched, from read_one_dim_file

    Returns:
        subs_type: wm_name or omega_name, the type of subsidence measurement
        subs: subsidence profile, shape (nsize,) [m/s or Pa/s]
    """
    names = [retVar.name for retVar in retVars]

    if wm_name in names and omega_name in names:
        raise ValueError("Sounding holds both " + wm_name + " and " + omega_name)

    if wm_name in names:
        subs_type = wm_name
    elif omega_name in names:
        subs_type = omega_name
    else:
        raise ValueError("Could not read vertical velocity compatable variable")
    subs = read_x_profile(nvar, nsize, subs_type, retVars)

    return subs_type, subs
