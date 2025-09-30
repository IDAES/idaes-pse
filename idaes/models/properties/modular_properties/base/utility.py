#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Common methods used by generic framework

Author: A Lee
"""
# TODO: Missing doc strings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from enum import Enum

from pyomo.environ import log, units as pyunits, value

from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
)
import idaes.logger as idaeslog
from idaes.core.scaling import CustomScalerBase

# Set up logger
_log = idaeslog.getLogger(__name__)


class StateIndex(Enum):
    true = 1
    apparent = 2


class GenericPropertyPackageError(PropertyPackageError):
    """
    Error message for when a property is called for but no option provided
    """

    def __init__(self, block, prop):
        super().__init__()
        self.prop = prop
        self.block = block

    def __str__(self):
        return (
            f"Generic Property Package instance {self.block} called for "
            f"{self.prop}, but was not provided with a method "
            f"for this property. Please add a method for this property "
            f"in the property parameter configuration."
        )


def get_method(self, config_arg, comp=None, phase=None, log_expression=False):
    """
    Method to inspect configuration argument and return the user-defined
    construction method associated with it.

    This method checks whether the value provided is a method or a
    module. If the value is a module, it looks in the module for a method
    with the same name as the config argument and returns this. If the
    value is a method, the method is returned. If the value is neither a
    module or a method, an ConfigurationError is raised.

    Args:
        config_arg : the configuration argument to look up
        comp : component name for which argument is to be retrieved
        phase : phase name indexing argument

    Returns:
        A callable method or a ConfigurationError
    """
    if comp is None:
        source_block = self.params.config
    else:
        source_block = self.params.get_component(comp).config

    try:
        c_arg = getattr(source_block, config_arg)
    except AttributeError:
        raise AttributeError(
            "{} Generic Property Package called for invalid "
            "configuration option {}. Please contact the "
            "developer of the property package.".format(self.name, config_arg)
        )

    if c_arg is None:
        raise GenericPropertyPackageError(self, config_arg)

    # Check to see if c_arg has an attribute with the name of the config_arg
    # If so, assume c_arg is a class or module holding property subclasses
    if hasattr(c_arg, config_arg):
        c_arg = getattr(c_arg, config_arg)
    if phase is not None:
        c_arg = c_arg[phase]

    # Try to get the return_expression method from c_arg
    # Otherwise assume c_arg is the return_expression method
    if not log_expression:
        try:
            mthd = c_arg.return_expression
        except AttributeError:
            mthd = c_arg
    else:
        if hasattr(c_arg, "return_log_expression"):
            mthd = c_arg.return_log_expression
        else:
            _log.warning(
                f"Failed to find log expression for {config_arg}."
                "Reverting to use of the log function of the "
                "non-log expression."
            )

            def mthd(*args, **kwargs):
                exp_mthd = get_method(
                    self, config_arg, comp=comp, phase=phase, log_expression=False
                )
                return log(exp_mthd(*args, **kwargs))

    # Call the return_expression method
    if callable(mthd):
        return mthd
    else:
        raise ConfigurationError(
            "{} Generic Property Package received invalid value "
            "for argument {}. Value must be a method, a class with a "
            "method named expression or a module containing one of the "
            "previous.".format(self.name, config_arg)
        )


def get_phase_method(self, config_arg, phase):
    """
    General method for finding and returning phase-specific configuration
    arguments.

    Args:
        config_arg : argument to find in Config block
        phase : phase in which to search for config_arg

    Returns:
        Pointer to method in Config block
    """
    p_config = self.params.get_phase(phase).config

    try:
        c_arg = getattr(p_config, config_arg)
    except AttributeError:
        raise AttributeError(
            "{} Generic Property Package called for invalid "
            "configuration option {}. Please contact the "
            "developer of the property package.".format(self.name, config_arg)
        )

    if c_arg is None:
        raise GenericPropertyPackageError(self, config_arg)

    # Check to see if c_arg has an attribute with the name of the config_arg
    # If so, assume c_arg is a class or module holding property subclasses
    if hasattr(c_arg, config_arg):
        c_arg = getattr(c_arg, config_arg)

    # Try to get the return_expression method from c_arg
    # Otherwise assume c_arg is the return_expression method
    try:
        mthd = c_arg.return_expression
    except AttributeError:
        mthd = c_arg

    # Check if method is callable
    if callable(mthd):
        return mthd
    else:
        raise ConfigurationError(
            "{} Generic Property Package received invalid value "
            "for argument {}. Value must be a method, a class with a "
            "method named expression or a module containing one of the "
            "previous.".format(self.name, config_arg)
        )

    return mthd


def get_component_object(self, comp):
    """
    Utility method to get a component object from the property parameter block.
    This code is used frequently throughout the generic property package
    libraries.

    Args:
        comp: name of the component object to be returned.

    Returns:
        Component: Component object with name comp.

    """
    return self.params.get_component(comp)


def get_bounds_from_config(b, state, base_units):
    """
    Method to take a 3- or 4-tuple state definition config argument and return
    tuples for the bounds and default value of the Var object.

    Expects the form (lower, default, upper, units) where units is optional

    Args:
        b - StateBlock on which the state vars are to be constructed
        state - name of state var as a string (to be matched with config dict)
        base_units - base units of state var to be used if conversion required

    Returns:
        bounds - 2-tuple of state var bounds in base units
        default_val - default value of state var in base units
    """
    try:
        var_config = b.params.config.state_bounds[state]
    except (KeyError, TypeError):
        # State definition missing
        return (None, None), None

    if len(var_config) == 4:
        # Units provided, need to convert values
        bounds = (
            pyunits.convert_value(
                var_config[0], from_units=var_config[3], to_units=base_units
            ),
            pyunits.convert_value(
                var_config[2], from_units=var_config[3], to_units=base_units
            ),
        )
        default_val = pyunits.convert_value(
            var_config[1], from_units=var_config[3], to_units=base_units
        )
    else:
        bounds = (var_config[0], var_config[2])
        default_val = var_config[1]

    return bounds, default_val


# Enumerate concentration form options
class ConcentrationForm(Enum):
    molarity = 1
    activity = 2
    molality = 3
    moleFraction = 4
    massFraction = 5
    partialPressure = 6


def get_concentration_term(blk, r_idx, log=False):
    """
    Get necessary concentration terms for reactions from property package, allowing for
    different bases.

    Args:
        blk: StateBlock of interest
        r_idx: index of reaction
        log: whether to use log concentration of not

    Returns:
        Var or Expression representing concentration

    """
    cfg = blk.params.config
    if "rate_reactions" in cfg:
        try:
            conc_form = cfg.rate_reactions[r_idx].concentration_form
        except KeyError:
            conc_form = cfg.equilibrium_reactions[r_idx].concentration_form
        state = blk.state_ref
    else:
        conc_form = cfg.inherent_reactions[r_idx].concentration_form
        state = blk

    # pylint: disable-next=protected-access
    if hasattr(state.params, "_electrolyte") and state.params._electrolyte:
        sub = "_true"
    else:
        sub = ""

    if log:
        pre = "log_"
    else:
        pre = ""

    if conc_form is None:
        raise ConfigurationError(
            "{} concentration_form configuration argument was not set. "
            "Please ensure that this argument is included in your "
            "configuration dict.".format(blk.name)
        )
    elif conc_form == ConcentrationForm.molarity:
        conc_term = getattr(state, pre + "conc_mol_phase_comp" + sub)
    elif conc_form == ConcentrationForm.activity:
        conc_term = getattr(state, pre + "act_phase_comp" + sub)
    elif conc_form == ConcentrationForm.molality:
        conc_term = getattr(state, pre + "molality_phase_comp" + sub)
    elif conc_form == ConcentrationForm.moleFraction:
        conc_term = getattr(state, pre + "mole_frac_phase_comp" + sub)
    elif conc_form == ConcentrationForm.massFraction:
        conc_term = getattr(state, pre + "mass_frac_phase_comp" + sub)
    elif conc_form == ConcentrationForm.partialPressure:
        conc_term = getattr(state, pre + "pressure_phase_comp" + sub)
    else:
        raise BurntToast(
            "{} get_concentration_term received unrecognised "
            "ConcentrationForm ({}). This should not happen - please contact "
            "the IDAES developers with this bug.".format(blk.name, conc_form)
        )

    return conc_term


def identify_VL_component_list(blk, phase_pair):
    """
    Identify liquid and vapor phases and which components are in VL equilibrium

    Args:
        blk: StateBlock of interest
        phase_pair: 2-tuple of Phases in equilibrium

    Returns:
        Lists of component names for:
            * liquid Phase object
            * vapor Phase object
            * components using Raoult's Law
            * components using Henry's Law
            * liquid only components,
            * vapor only components

    """
    vl_comps = []
    henry_comps = []
    l_only_comps = []
    v_only_comps = []

    pparams = blk.params

    if pparams.get_phase(phase_pair[0]).is_liquid_phase():
        l_phase = phase_pair[0]
        if pparams.get_phase(phase_pair[1]).is_vapor_phase():
            v_phase = phase_pair[1]
        else:
            raise PropertyPackageError(
                f"Phase pair {phase_pair[0]}-{phase_pair[1]} was identified as "
                f"being a VLE pair, however {phase_pair[0]} is liquid but "
                f"{phase_pair[1]} is not vapor."
            )
    elif pparams.get_phase(phase_pair[1]).is_liquid_phase():
        l_phase = phase_pair[1]
        if pparams.get_phase(phase_pair[0]).is_vapor_phase():
            v_phase = phase_pair[0]
        else:
            raise PropertyPackageError(
                f"Phase pair {phase_pair[0]}-{phase_pair[1]} was identified as "
                f"being a VLE pair, however {phase_pair[1]} is liquid but "
                f"{phase_pair[0]} is not vapor."
            )
    else:
        raise PropertyPackageError(
            f"Phase pair {phase_pair[0]}-{phase_pair[1]} was identified as "
            f"being a VLE pair, however neither {phase_pair[0]} nor "
            f"{phase_pair[1]} is liquid."
        )

    for j in blk.params.component_list:
        if (l_phase, j) in blk.phase_component_set and (
            v_phase,
            j,
        ) in blk.phase_component_set:
            cobj = pparams.get_component(j)
            if cobj.config.henry_component is not None and (
                phase_pair[0] in cobj.config.henry_component
                or phase_pair[1] in cobj.config.henry_component
            ):
                henry_comps.append(j)
            else:
                vl_comps.append(j)
        elif (l_phase, j) in blk.phase_component_set:
            l_only_comps.append(j)
        elif (v_phase, j) in blk.phase_component_set:
            v_only_comps.append(j)

    if len(vl_comps) == 0 and len(henry_comps) == 0:
        raise PropertyPackageError(
            f"Phase pair {phase_pair[0]}-{phase_pair[1]} was identified as "
            "being a VLE pair, however there are no components present in "
            "both the vapor and liquid phases simultaneously."
        )

    return l_phase, v_phase, vl_comps, henry_comps, l_only_comps, v_only_comps


TOL = 1e-1
MAX_ITER = 30


def estimate_Tbub(blk, T_units, raoult_comps, henry_comps, liquid_phase):
    """
    Function to estimate bubble point temperature

    Args:
        blk: StateBlock to use
        T_units: units of temperature
        raoult_comps: list of components that follow Raoult's Law
        henry_comps: list of components that follow Henry's Law
        liquid_phase: name of liquid phase

    Returns:
        Estimated bubble point temperature as a float.

    """
    # Use lowest component temperature_crit as starting point
    # Starting high and moving down generally works better,
    # as it under-predicts next step due to exponential form of
    # Psat.
    # Subtract 1 to avoid potential singularities at Tcrit
    Tbub0 = (
        min(blk.params.get_component(j).temperature_crit.value for j in raoult_comps)
        - 1
    )

    err = 1
    counter = 0

    # Newton solver with step limiter to prevent overshoot
    # Tolerance only needs to be ~1e-1
    # Iteration limit of 30
    while err > TOL and counter < MAX_ITER:
        f = value(
            sum(
                get_method(blk, "pressure_sat_comp", j)(
                    blk, blk.params.get_component(j), Tbub0 * T_units
                )
                * blk.mole_frac_comp[j]
                for j in raoult_comps
            )
            + sum(
                blk.mole_frac_comp[j]
                * blk.params.get_component(j)
                .config.henry_component[liquid_phase]["method"]
                .return_expression(blk, liquid_phase, j, Tbub0 * T_units)
                for j in henry_comps
            )
            - blk.pressure
        )
        df = value(
            sum(
                get_method(blk, "pressure_sat_comp", j)(
                    blk, blk.params.get_component(j), Tbub0 * T_units, dT=True
                )
                * blk.mole_frac_comp[j]
                for j in raoult_comps
            )
            + sum(
                blk.mole_frac_comp[j]
                * blk.params.get_component(j)
                .config.henry_component[liquid_phase]["method"]
                .dT_expression(blk, liquid_phase, j, Tbub0 * T_units)
                for j in henry_comps
            )
        )

        # Limit temperature step to avoid excessive overshoot
        if f / df < -50:
            Tbub1 = Tbub0 + 50
        elif f / df > 50:
            Tbub1 = Tbub0 - 50
        else:
            Tbub1 = Tbub0 - f / df

        err = abs(Tbub1 - Tbub0)
        Tbub0 = Tbub1
        counter += 1

    return Tbub0


def estimate_Tdew(blk, T_units, raoult_comps, henry_comps, liquid_phase):
    """
    Function to estimate dew point temperature

    Args:
        blk: StateBlock to use
        T_units: units of temperature
        raoult_comps: list of components that follow Raoult's Law
        henry_comps: list of components that follow Henry's Law
        liquid_phase: name of liquid phase

    Returns:
        Estimated dew point temperature as a float.

    """
    # Use lowest component critical temperature
    # as starting point
    # Subtract 1 to avoid potential singularities at Tcrit
    Tdew0 = (
        min(blk.params.get_component(j).temperature_crit.value for j in raoult_comps)
        - 1
    )

    err = 1
    counter = 0

    # Newton solver with step limiter to prevent overshoot
    # Tolerance only needs to be ~1e-1
    # Iteration limit of 30
    while err > TOL and counter < MAX_ITER:
        f = value(
            blk.pressure
            * (
                sum(
                    blk.mole_frac_comp[j]
                    / get_method(blk, "pressure_sat_comp", j)(
                        blk, blk.params.get_component(j), Tdew0 * T_units
                    )
                    for j in raoult_comps
                )
                + sum(
                    blk.mole_frac_comp[j]
                    / blk.params.get_component(j)
                    .config.henry_component[liquid_phase]["method"]
                    .return_expression(blk, liquid_phase, j, Tdew0 * T_units)
                    for j in henry_comps
                )
            )
            - 1
        )
        df = -value(
            blk.pressure
            * (
                sum(
                    blk.mole_frac_comp[j]
                    / get_method(blk, "pressure_sat_comp", j)(
                        blk, blk.params.get_component(j), Tdew0 * T_units
                    )
                    ** 2
                    * get_method(blk, "pressure_sat_comp", j)(
                        blk,
                        blk.params.get_component(j),
                        Tdew0 * T_units,
                        dT=True,
                    )
                    for j in raoult_comps
                )
                + sum(
                    blk.mole_frac_comp[j]
                    / blk.params.get_component(j)
                    .config.henry_component[liquid_phase]["method"]
                    .return_expression(blk, liquid_phase, j, Tdew0 * T_units)
                    ** 2
                    * blk.params.get_component(j)
                    .config.henry_component[liquid_phase]["method"]
                    .dT_expression(blk, liquid_phase, j, Tdew0 * T_units)
                    for j in henry_comps
                )
            )
        )

        # Limit temperature step to avoid excessive overshoot
        if f / df < -50:
            Tdew1 = Tdew0 + 50
        elif f / df > 50:
            Tdew1 = Tdew0 - 50
        else:
            Tdew1 = Tdew0 - f / df

        err = abs(Tdew1 - Tdew0)
        Tdew0 = Tdew1
        counter += 1

    return Tdew0


def estimate_Pbub(blk, raoult_comps, henry_comps, liquid_phase):
    """
    Function to estimate bubble point pressure

    Args:
        blk: StateBlock to use
        raoult_comps: list of components that follow Raoult's Law
        henry_comps: list of components that follow Henry's Law
        liquid_phase: name of liquid phase

    Returns:
        Estimated bubble point pressure as a float.

    """
    return value(
        sum(blk.mole_frac_comp[j] * blk.pressure_sat_comp[j] for j in raoult_comps)
        + sum(blk.mole_frac_comp[j] * blk.henry[liquid_phase, j] for j in henry_comps)
    )


def estimate_Pdew(blk, raoult_comps, henry_comps, liquid_phase):
    """
    Function to estimate dew point pressure

    Args:
        blk: StateBlock to use
        raoult_comps: list of components that follow Raoult's Law
        henry_comps: list of components that follow Henry's Law
        liquid_phase: name of liquid phase

    Returns:
        Estimated dew point pressure as a float.

    """
    # Safety catch for cases where Psat or Henry's constant might be 0
    # Not sure if this is meaningful, but if this is true then mathematically Pdew = 0
    if any(value(blk.pressure_sat_comp[j]) == 0 for j in raoult_comps) or any(
        value(blk.henry[liquid_phase, j]) == 0 for j in henry_comps
    ):
        return 0
    return value(
        1
        / (
            sum(blk.mole_frac_comp[j] / blk.pressure_sat_comp[j] for j in raoult_comps)
            + sum(
                blk.mole_frac_comp[j] / blk.henry[liquid_phase, j] for j in henry_comps
            )
        )
    )


class ModularPropertiesScalerBase(CustomScalerBase):
    """
    Base class for ModularPropertiesScaler and ModularReactionsScaler.
    Handles the logic for calling scaling methods for each individual module.
    """

    def call_module_scaling_method(
        self, model, module, index, method, overwrite: bool = False
    ):
        try:
            scaler_class = module.default_scaler
        # TODO create interface where the user can provide custom scalers for individual modules
        except AttributeError:
            _log.debug(
                f"No default Scaler set for module {module}. Cannot call {method}."
            )
            return
        scaler_obj = scaler_class(**self.CONFIG)
        try:
            method_func = getattr(scaler_obj, method)
        except AttributeError as err:
            raise AttributeError(
                f"Could not find {method} method on scaler for module {module}."
            ) from err
        if index is None:
            method_func(model, overwrite=overwrite)
        else:
            method_func(model, index, overwrite=overwrite)
