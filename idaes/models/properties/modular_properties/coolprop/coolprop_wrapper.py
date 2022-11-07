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
within the IDAES generic properties framework.

The CoolPropWrapper class contains a set of sub-classes for the supported
thermophysical properties required by the generic property framework, along
with some helper functions for common calls to the CoolProp database and
chaching data to avoid repeated database lookups.
"""
import json

from pyomo.common.dependencies import attempt_import
from pyomo.environ import units as pyunits, Var

import idaes.models.properties.modular_properties.coolprop.coolprop_forms as cforms
from idaes.core.util.exceptions import BurntToast

CoolProp, coolprop_available = attempt_import("CoolProp.CoolProp")


# Map IDAES names to Coolprop names and where to look for them
name_map = {
    "dens_mol_crit": ("rhomolar", "CRITICAL"),
    "enth_mol_crit": ("hmolar", "CRITICAL"),
    "entr_mol_crit": ("smolar", "CRITICAL"),
    "mw": ("molar_mass", "EOS"),
    "omega": ("acentric", "EOS"),
    "pressure_crit": ("p", "CRITICAL"),
    "temperature_crit": ("T", "CRITICAL"),
}


class CoolPropExpressionError(ValueError):
    """
    Error message for when an unexpected expression form is called for.

    This is mostly to future proof the code against changes in CoolProp where
    the expression form is changed to something we don't recognize.
    """

    def __init__(self, prop, comp):
        msg = (
            f"Found unsupported expression form for {prop} "
            f"of component {comp}. This likely occured due to "
            f"changes in CoolProp and the interface should be "
            f"updated."
        )

        super().__init__(self, msg)


class CoolPropPropertyError(KeyError):
    """
    Error message for when a parameter is called for that is not in CoolProp's
    database.
    """

    def __init__(self, prop, comp):
        msg = (
            f"Could not retrieve parameters for {prop} of "
            f"component {comp} from CoolProp. This likely indicates "
            f"that CoolProp does not have values for the necessary "
            f"parameters."
        )

        super().__init__(self, msg)


class CoolPropWrapper:
    """
    Interface wrapper for calling CoolProp parameter database from within the
    IDAES Generic Property Package Framework.

    This class is intended to be used in theromphysical property definition
    dicts and directs the properties framework to use forms and parameter data
    from the CoolProp libraries (if available). This requires that the user
    have CoolProp installed locally.
    """

    _cached_components = {}

    @staticmethod
    def get_parameter_value(comp_name, param):
        """
        Retrieve tuples of (value, units) for the specified property and
        component for the critical state section of the CoolProp json format.

        Args:
            comp_name: name of componet for which to retirve parameters
            param: name of parameter to return value and units

        Returns:
            tuple of parameter (value, units)
        """
        map_tuple = name_map.get(param, None)
        if map_tuple is None:
            raise BurntToast(
                "Unrecognized property name in CoolProp wrapper. Please "
                "contact the IDAES developers with this bug."
            )

        if map_tuple[1] == "CRITICAL":
            ptuple = CoolPropWrapper._get_critical_property(comp_name, map_tuple[0])
        elif map_tuple[1] == "EOS":
            ptuple = CoolPropWrapper._get_eos_property(comp_name, map_tuple[0])
        else:
            raise BurntToast(
                "Unrecognized argument in CoolProp wrapper. Please contact "
                "the IDAES developers with this bug."
            )

        return ptuple

    @staticmethod
    def flush_cached_components():
        """
        Clear cached component parameter data. This sets _cached_components to
        an empty dict.

        Args:
            None

        Returns:
            None
        """
        CoolPropWrapper._cached_components = {}

    # -------------------------------------------------------------------------
    # Pure component property sub-classes
    class dens_mol_liq_comp:
        """
        Calculate liquid molar density using CoolProp forms and parameters.
        """

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            # 8-Dec-21: CoolProp only uses rhoLnoexp for liquid density
            nlist, tlist = CoolPropWrapper._get_param_dicts(
                cname, cdict, "dens_mol_liq_comp", "rhoL", ["rhoLnoexp"]
            )

            cforms.parameters_nt_sum(cobj, "dens_mol_liq_comp", nlist, tlist)

        @staticmethod
        def return_expression(b, cobj, T):
            return cforms.expression_nonexponential(
                cobj, "dens_mol_liq_comp", T, cobj.dens_mol_crit
            )

    class enth_mol_liq_comp:
        """
        Calculate liquid molar enthalpy using CoolProp forms and parameters.
        """

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            # 29-Dec-21: CoolProp only uses rational_polynomial.
            Alist, Blist = CoolPropWrapper._get_param_dicts(
                cname, cdict, "enth_mol_liq_comp", "hL", ["rational_polynomial"]
            )
            href = cdict["EOS"][0]["STATES"]["hs_anchor"]["hmolar"]

            cforms.parameters_polynomial(
                cobj, "enth_mol_liq_comp", pyunits.J / pyunits.mol, Alist, Blist
            )

            href_var = Var(
                doc="Reference heat of formation", units=pyunits.J / pyunits.mol
            )
            cobj.add_component("enth_mol_liq_comp_anchor", href_var)
            href_var.fix(href)

        @staticmethod
        def return_expression(b, cobj, T):
            h = (
                cforms.expression_polynomial(cobj, "enth_mol_liq_comp", T)
                + cobj.enth_mol_liq_comp_anchor
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(h, units.ENERGY_MOLE)

    class enth_mol_ig_comp:
        """
        Calculate ideal gas molar enthalpy using CoolProp forms and parameters.
        """

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            # 29-Dec-21: CoolProp only uses rational_polynomial.
            Alist, Blist = CoolPropWrapper._get_param_dicts(
                cname, cdict, "enth_mol_ig_comp", "hLV", ["rational_polynomial"]
            )

            cforms.parameters_polynomial(
                cobj, "enth_mol_ig_comp", pyunits.J / pyunits.mol, Alist, Blist
            )

            # Next, build parameters for enth_mol_liq_comp if necessary
            if not hasattr(cobj, "enth_mol_liq_comp_coeff_A0"):
                CoolPropWrapper.enth_mol_liq_comp.build_parameters(cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            h = cforms.expression_polynomial(
                cobj, "enth_mol_ig_comp", T
            ) + CoolPropWrapper.enth_mol_liq_comp.return_expression(b, cobj, T)

            units = b.params.get_metadata().derived_units
            return pyunits.convert(h, units.ENERGY_MOLE)

    class entr_mol_liq_comp:
        """
        Calculate liquid molar entropy using CoolProp forms and parameters.
        """

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            # 29-Dec-21: CoolProp only uses rational_polynomial.
            Alist, Blist = CoolPropWrapper._get_param_dicts(
                cname, cdict, "entr_mol_liq_comp", "sL", ["rational_polynomial"]
            )
            sref = cdict["EOS"][0]["STATES"]["hs_anchor"]["smolar"]

            cforms.parameters_polynomial(
                cobj,
                "entr_mol_liq_comp",
                pyunits.J / pyunits.mol / pyunits.K,
                Alist,
                Blist,
            )

            sref_var = Var(
                doc="Reference heat of formation",
                units=pyunits.J / pyunits.mol / pyunits.K,
            )
            cobj.add_component("entr_mol_liq_comp_anchor", sref_var)
            sref_var.fix(sref)

        @staticmethod
        def return_expression(b, cobj, T):
            s = (
                cforms.expression_polynomial(cobj, "entr_mol_liq_comp", T)
                + cobj.entr_mol_liq_comp_anchor
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(s, units.ENTROPY_MOLE)

    class entr_mol_ig_comp:
        """
        Calculate ideal gas molar entropy using CoolProp forms and parameters.
        """

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            # 29-Dec-21: CoolProp only uses rational_polynomial.
            Alist, Blist = CoolPropWrapper._get_param_dicts(
                cname, cdict, "entr_mol_ig_comp", "sLV", ["rational_polynomial"]
            )

            cforms.parameters_polynomial(
                cobj,
                "entr_mol_ig_comp",
                pyunits.J / pyunits.mol / pyunits.K,
                Alist,
                Blist,
            )

            # Next, build parameters for entr_mol_liq_comp if necessary
            if not hasattr(cobj, "entr_mol_liq_comp_coeff_A0"):
                CoolPropWrapper.entr_mol_liq_comp.build_parameters(cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            s = cforms.expression_polynomial(
                cobj, "entr_mol_ig_comp", T
            ) + CoolPropWrapper.entr_mol_liq_comp.return_expression(b, cobj, T)

            units = b.params.get_metadata().derived_units
            return pyunits.convert(s, units.ENTROPY_MOLE)

    class pressure_sat_comp:
        """
        Calculate pure component saturation pressure using CoolProp forms and
        parameters.
        """

        @staticmethod
        def build_parameters(cobj):
            cname = cobj.local_name
            cdict = CoolPropWrapper._get_component_data(cname)

            # 8-Dec-21: CoolProp only has two "types" for pressure_sat
            # which appear to be equivalent.
            nlist, tlist = CoolPropWrapper._get_param_dicts(
                cname, cdict, "pressure_sat", "pS", ["pL", "pV"], using_tau_r=True
            )

            cforms.parameters_nt_sum(cobj, "pressure_sat", nlist, tlist)

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            if dT:
                return CoolPropWrapper.pressure_sat_comp.dT_expression(b, cobj, T)

            return cforms.expression_exponential(
                cobj, "pressure_sat", T, cobj.pressure_crit, tau=True
            )

        @staticmethod
        def dT_expression(b, cobj, T):
            return cforms.dT_expression_exponential(
                cobj, "pressure_sat", T, cobj.pressure_crit, tau=True
            )

    # -------------------------------------------------------------------------
    # Internal methods

    @staticmethod
    def _get_component_data(comp_name):
        """
        Method to get parameter data from cached components, and to call
        _load_component_data if it is not present. This method also includes
        checks to handle aliases of components.

        Args:
            comp_name: name of component to get data for

        Returns:
            dict constructed from json string retrieved from CoolProp database.
        """
        if comp_name in CoolPropWrapper._cached_components:
            # First check to see if component present by comp_name
            return CoolPropWrapper._cached_components[comp_name]
        else:
            # Check to see if comp_name is an alias for a cached component
            for v in CoolPropWrapper._cached_components.values():
                if comp_name in v["INFO"]["ALIASES"] or comp_name in v["INFO"]["NAME"]:
                    CoolPropWrapper._cached_components[comp_name] = v
                    return v

        # If we haven't returned yet, then we need to load the component
        return CoolPropWrapper._load_component_data(comp_name)

    @staticmethod
    def _load_component_data(comp_name):
        """
        Method to load parameter data for specified component from the CoolProp
        database in json format. Loaded data is in dict form and is stored in
        _cached_components to avoid need for repeated calls to CoolProp.

        Args:
            comp_name: name of component ot retrieve parameters for.

        Returns:
            dict constructed from json string retrieved from CoolProp database.

        Raises:
            RuntimeError is component is not found in database
        """
        try:
            prop_str = CoolProp.get_fluid_param_string(comp_name, "JSON")
        except RuntimeError:
            raise RuntimeError(
                f"Failed to find component {comp_name} in CoolProp JSON " f"database."
            )
        comp_prop = json.loads(prop_str)[0]

        CoolPropWrapper._cached_components[comp_name] = comp_prop

        return comp_prop

    @staticmethod
    def _get_critical_property(comp_name, prop_name):
        """
        Method to retrieve tuples of (value, units) for the specified property
        and component for the critical state section of the CoolProp json
        format.
        """
        cdict = CoolPropWrapper._get_component_data(comp_name)

        pc = cdict["STATES"]["critical"][prop_name]
        punits = getattr(pyunits, cdict["STATES"]["critical"][prop_name + "_units"])

        return (pc, punits)

    @staticmethod
    def _get_eos_property(comp_name, prop_name):
        """
        Method to retrieve tuples of (value, units) for the specified property
        and component for the EoS section of the CoolProp json format.
        """
        cdict = CoolPropWrapper._get_component_data(comp_name)

        pc = cdict["EOS"][0][prop_name]
        punits = cdict["EOS"][0][prop_name + "_units"]
        if punits == "-":
            punits = pyunits.dimensionless
        else:
            punits = getattr(pyunits, punits)

        return (pc, punits)

    @staticmethod
    def _get_param_dicts(
        comp_name,
        comp_data,
        prop_name,
        coolprop_name,
        expected_forms,
        using_tau_r=False,
    ):
        """
        Method to get parameter sets for expression forms. Also includes check
        to verify the expression form listed by CoolProp mathces the expected
        for in IDAES.

        Args:
            comp_name: name of the component for which parameters are required
            comp_data: dict of parameter values for component (from CoolProp)
            prop_name: IDAES name of property for which parameters are required
            coolprop_name: name used by CoolProp for property of interest
            expected_forms: list of expression forms supported by wrapper for
                              for property
            using_tau_r: (optional) flag indicating whether to check for use
                          of tau_r in the epxression or not (default = False,
                          do not check).

        Returns:
            tuple of lists of parameters required by expression
        """
        try:
            # First, check to make sure the listed expression form is
            # supported.
            if comp_data["ANCILLARIES"][coolprop_name]["type"] not in expected_forms:
                # If not one of the forms we recognise, raise an exception
                raise CoolPropExpressionError(prop_name, comp_name)
            elif (
                using_tau_r
                and not comp_data["ANCILLARIES"][coolprop_name]["using_tau_r"]
            ):
                # If not one of the forms we recognise, raise an exception
                raise CoolPropExpressionError(prop_name, comp_name)

            if expected_forms == ["rational_polynomial"]:
                list1 = comp_data["ANCILLARIES"][coolprop_name]["A"]
                list2 = comp_data["ANCILLARIES"][coolprop_name]["B"]
            else:
                list1 = comp_data["ANCILLARIES"][coolprop_name]["n"]
                list2 = comp_data["ANCILLARIES"][coolprop_name]["t"]
        except KeyError:
            raise CoolPropPropertyError(prop_name, comp_name)

        return (list1, list2)
