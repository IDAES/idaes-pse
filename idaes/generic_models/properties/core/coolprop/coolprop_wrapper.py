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
This module contains methods for looking up and using CoolProp parameters
within the IDAES generic proeprties framework.
"""
import json

try:
    from CoolProp import CoolProp
except ModuleNotFoundError:
    raise ModuleNotFoundError(
        "The optional dependency CoolProp is required to use the IDAES "
        "CoolProp wrapper. Please consult the installation instructions for "
        "how to install CoolProp.")

from pyomo.environ import units as pyunits

import idaes.generic_models.properties.core.coolprop.coolprop_forms as cforms
from idaes.core.util.exceptions import BurntToast


# Map IDAES names to Coolprop names and where to look for them
name_map = {"dens_mol_crit": ("rhomolar", "CRITICAL"),
            "enth_mol_crit": ("hmolar", "CRITICAL"),
            "entr_mol_crit": ("smolar", "CRITICAL"),
            "mw": ("molar_mass", "EOS"),
            "omega": ("acentric", "EOS"),
            "pressure_crit": ("p", "CRITICAL"),
            "temperature_crit": ("T", "CRITICAL")}


class CoolPropExpressionError(ValueError):
    # Error message for when an unexpected expression form is used
    def __init__(self, prop, comp):
        self.prop = prop
        self.comp = comp

    def __str__(self):
        return (f"Found unsupported expression form for {self.prop} "
                f"of component {self.comp}. This likely occured due to "
                f"changes in CoolProp and the interface should be "
                f"updated.")


class CoolPropPropertyError(KeyError):
    # Error message for when CoolProp is missing entry for a property
    def __init__(self, prop, comp):
        self.prop = prop
        self.comp = comp

    def __str__(self):
        return (f"Could not retrieve parameters for {self.prop} of "
                f"component {self.comp} from CoolProp. This likely indicates "
                f"that CoolProp does not have values for the necessary "
                f"parameters.")


class CoolPropWrapper:
    cached_components = {}

    @staticmethod
    def get_parameter_value(comp_name, param):
        try:
            map_tuple = name_map[param]
        except KeyError:
            raise BurntToast(
                "Unrecognized property name in CoolProp wrapper. Please "
                "contact the IDAES developers with this bug.")

        if map_tuple[1] == "CRITICAL":
            ptuple = CoolPropWrapper._get_critical_property(
                comp_name, map_tuple[0])
        elif map_tuple[1] == "EOS":
            ptuple = CoolPropWrapper._get_eos_property(
                comp_name, map_tuple[0])
        else:
            raise BurntToast(
                "Unrecognized argument in CoolProp wrapper. Please contact "
                "the IDAES developers with this bug.")

        return ptuple

    @staticmethod
    def flush_cached_components():
        CoolPropWrapper.cached_components = {}

    # -------------------------------------------------------------------------
    # Pure component property methods
    class pressure_sat_comp():

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            try:
                # First, check to make sure the listed expression form is
                # supported.
                # 8-Dec-21: CoolProp only has two "types" for pressure_sat
                # which appear to be equivalent.
                if (cdict["ANCILLARIES"]["pS"]["type"] not in ["pL", "pV"] or
                        not cdict["ANCILLARIES"]["pS"]["using_tau_r"]):
                    # If not one of the forms we recognise, raise an exception
                    raise CoolPropExpressionError("pressure_sat", cname)

                ndict = cdict["ANCILLARIES"]["pS"]["n"]
                tdict = cdict["ANCILLARIES"]["pS"]["t"]
            except KeyError:
                raise CoolPropPropertyError("pressure_sat", cname)

            cforms.parameters_exponential(cobj, "pressure_sat", ndict, tdict)

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            if dT:
                return CoolPropWrapper.pressure_sat_comp.dT_expression(
                    b, cobj, T)

            return cforms.expression_exponential_tau(
                cobj, "pressure_sat", T, cobj.pressure_crit)

        @staticmethod
        def dT_expression(b, cobj, T):
            return cforms.dT_expression_exponential_tau(
                cobj, "pressure_sat", T, cobj.pressure_crit)

    # -------------------------------------------------------------------------
    # Internal methods

    @staticmethod
    def _get_component_data(comp_name):
        if comp_name in CoolPropWrapper.cached_components:
            # First check to see if component present by comp_name
            return CoolPropWrapper.cached_components[comp_name]
        else:
            # Check to see if comp_name is an alias for a cached component
            for v in CoolPropWrapper.cached_components.values():
                if (comp_name in v["INFO"]["ALIASES"] or
                        comp_name in v["INFO"]["NAME"]):
                    return v

        # If we haven't returned yet, then we need to load the component
        return CoolPropWrapper._load_component_data(comp_name)

    @staticmethod
    def _load_component_data(comp_name):
        try:
            prop_str = CoolProp.get_fluid_param_string(comp_name, "JSON")
        except RuntimeError:
            raise RuntimeError(
                f"Failed to find component {comp_name} in CoolProp JSON "
                f"database.")
        comp_prop = json.loads(prop_str)[0]

        CoolPropWrapper.cached_components[comp_name] = comp_prop

        return comp_prop

    @staticmethod
    def _get_critical_property(comp_name, prop_name):
        cdict = CoolPropWrapper._get_component_data(comp_name)

        pc = cdict["STATES"]["critical"][prop_name]
        punits = getattr(pyunits, cdict["STATES"]["critical"][
            prop_name+"_units"])

        return (pc, punits)

    @staticmethod
    def _get_eos_property(comp_name, prop_name):
        cdict = CoolPropWrapper._get_component_data(comp_name)

        pc = cdict["EOS"][0][prop_name]
        punits = cdict["EOS"][0][prop_name+"_units"]
        if punits == "-":
            punits = pyunits.dimensionless
        else:
            punits = getattr(pyunits, punits)

        return (pc, punits)
