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


class CoolPropWrapper:
    cached_components = {}

    @staticmethod
    def flush_cached_components():
        CoolPropWrapper.cached_components = {}

    @staticmethod
    def get_component_data(comp_name):
        if comp_name in CoolPropWrapper.cached_components:
            # First check to see if component present by comp_name
            return CoolPropWrapper.cached_components[comp_name]
        else:
            # Check to see if comp_name is an alias for a cached component
            for v in CoolPropWrapper.cached_components.values():
                if comp_name in v["INFO"]["ALIASES"]:
                    return v

        # If we haven't returned yet, then we need to load the component
        return CoolPropWrapper._load_component_data(comp_name)

    @staticmethod
    def _load_component_data(comp_name):
        # TODO : Handle alternative names for components
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
    def get_critical_property(comp_name, prop_name):
        cdict = CoolPropWrapper.get_component_data(comp_name)

        # Map IDAES names to Coolprop names
        name_map = {"dens_mol_crit": "rhomolar",
                    "enth_mol_crit": "hmolar",
                    "entr_mol_crit": "smolar",
                    "pressure_crit": "p",
                    "temperature_crit": "T"}

        pc = cdict["STATES"]["critical"][name_map[prop_name]]
        punits = getattr(pyunits, cdict["STATES"]["critical"][
            name_map[prop_name]+"_units"])

        return (pc, punits)

    @staticmethod
    def get_eos_property(comp_name, prop_name):
        cdict = CoolPropWrapper.get_component_data(comp_name)

        # Map IDAES names to Coolprop names
        name_map = {"omega": "acentric",
                    "mw": "molar_mass"}

        pc = cdict["EOS"][0][name_map[prop_name]]
        punits = cdict["EOS"][0][name_map[prop_name]+"_units"]
        if punits == "-":
            punits = pyunits.dimensionless
        else:
            punits = getattr(pyunits, punits)

        return (pc, punits)
