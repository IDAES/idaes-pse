#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Common methods used by generic framework

Author: A Lee
"""
from enum import Enum

from pyomo.environ import units as pyunits

from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
)
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


class StateIndex(Enum):
    true = 1
    apparent = 2


class GenericPropertyPackageError(PropertyPackageError):
    # Error message for when a property is called for but no option provided
    def __init__(self, block, prop):
        self.prop = prop
        self.block = block

    def __str__(self):
        return (
            f"Generic Property Package instance {self.block} called for "
            f"{self.prop}, but was not provided with a method "
            f"for this property. Please add a method for this property "
            f"in the property parameter configuration."
        )


def get_method(self, config_arg, comp=None, phase=None):
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
    try:
        mthd = c_arg.return_expression
    except AttributeError:
        mthd = c_arg

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
    This code is used frequently throughout the generic property pacakge
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
