##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 3rd Edition
Reid, Prausnitz and Polling, 1977, McGraw-Hill

All parameter indicies based on conventions used by the source
"""

from pyomo.environ import exp, log, Var, units as pyunits

from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class cp_mol_ig_comp():

    @staticmethod
    def build_parameters(cobj):
        cobj.cp_mol_ig_comp_coeff_A = Var(
                doc="Parameter A for ideal gas molar heat capacity",
                units=pyunits.cal/pyunits.mol/pyunits.K)
        set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="A")

        cobj.cp_mol_ig_comp_coeff_B = Var(
                doc="Parameter B for ideal gas molar heat capacity",
                units=pyunits.cal/pyunits.mol/pyunits.K**2)
        set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="B")

        cobj.cp_mol_ig_comp_coeff_C = Var(
                doc="Parameter C for ideal gas molar heat capacity",
                units=pyunits.cal/pyunits.mol/pyunits.K**3)
        set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="C")

        cobj.cp_mol_ig_comp_coeff_D = Var(
                doc="Parameter D for ideal gas molar heat capacity",
                units=pyunits.cal/pyunits.mol/pyunits.K**4)
        set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="D")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific heat capacity
        T = pyunits.convert(T, to_units=pyunits.K)

        cp = (cobj.cp_mol_ig_comp_coeff_D*T**3 +
              cobj.cp_mol_ig_comp_coeff_C*T**2 +
              cobj.cp_mol_ig_comp_coeff_B*T +
              cobj.cp_mol_ig_comp_coeff_A)

        base_units = b.params.get_metadata().default_units
        cp_units = (base_units["mass"] *
                    base_units["length"]**2 *
                    base_units["time"]**-2 *
                    base_units["amount"]**-1 *
                    base_units["temperature"]**-1)
        return pyunits.convert(cp, cp_units)


class enth_mol_ig_comp():

    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
            cp_mol_ig_comp.build_parameters(cobj)

        base_units = cobj.parent_block().get_metadata().default_units
        h_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        cobj.enth_mol_form_vap_comp_ref = Var(
                doc="Vapor phase molar heat of formation @ Tref",
                units=h_units)
        set_param_from_config(cobj, param="enth_mol_form_vap_comp_ref")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific enthalpy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        base_units = b.params.get_metadata().default_units
        h_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        h = (pyunits.convert(
                (cobj.cp_mol_ig_comp_coeff_D/4)*(T**4-Tr**4) +
                (cobj.cp_mol_ig_comp_coeff_C/3)*(T**3-Tr**3) +
                (cobj.cp_mol_ig_comp_coeff_B/2)*(T**2-Tr**2) +
                cobj.cp_mol_ig_comp_coeff_A*(T-Tr), h_units) +
             cobj.enth_mol_form_vap_comp_ref)

        return h


class entr_mol_ig_comp():

    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
            cp_mol_ig_comp.build_parameters(cobj)

        base_units = cobj.parent_block().get_metadata().default_units
        s_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1 *
                   base_units["temperature"]**-1)

        cobj.entr_mol_form_vap_comp_ref = Var(
                doc="Vapor phase molar entropy of formation @ Tref",
                units=s_units)
        set_param_from_config(cobj, param="entr_mol_form_vap_comp_ref")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific entropy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        base_units = b.params.get_metadata().default_units
        s_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1 *
                   base_units["temperature"]**-1)

        s = (pyunits.convert(
                (cobj.cp_mol_ig_comp_coeff_D/3)*(T**3-Tr**3) +
                (cobj.cp_mol_ig_comp_coeff_C/2)*(T**2-Tr**2) +
                cobj.cp_mol_ig_comp_coeff_B*(T-Tr) +
                cobj.cp_mol_ig_comp_coeff_A*log(T/Tr), s_units) +
             cobj.entr_mol_form_vap_comp_ref)

        return s


# -----------------------------------------------------------------------------
# Saturation pressure
# Note that this equation in not valid beyond the critical temperature
class pressure_sat_comp():

    @staticmethod
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff_A = Var(
                doc="Antoine A coefficient for calculating Psat",
                units=pyunits.dimensionless)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="A")

        cobj.pressure_sat_comp_coeff_B = Var(
                doc="Antoine B coefficient for calculating Psat",
                units=pyunits.K)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="B")

        cobj.pressure_sat_comp_coeff_C = Var(
                doc="Antoine C coefficient for calculating Psat",
                units=pyunits.K)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C")

    @staticmethod
    def return_expression(b, cobj, T, dT=False):
        if dT:
            return pressure_sat_comp.dT_expression(b, cobj, T)

        psat = (exp(cobj.pressure_sat_comp_coeff_A -
                    cobj.pressure_sat_comp_coeff_B /
                    (pyunits.convert(T, to_units=pyunits.K) +
                     cobj.pressure_sat_comp_coeff_C)))*pyunits.mmHg

        base_units = b.params.get_metadata().default_units
        p_units = (base_units["mass"] *
                   base_units["length"]**-1 *
                   base_units["time"]**-2)
        return pyunits.convert(psat, to_units=p_units)

    @staticmethod
    def dT_expression(b, cobj, T):
        p_sat_dT = (pressure_sat_comp.return_expression(b, cobj, T) *
                    cobj.pressure_sat_comp_coeff_B /
                    (pyunits.convert(T, to_units=pyunits.K) +
                     cobj.pressure_sat_comp_coeff_C)**2)

        base_units = b.params.get_metadata().default_units
        dp_units = (base_units["mass"] *
                    base_units["length"]**-1 *
                    base_units["time"]**-2 *
                    base_units["temperature"]**-1)
        return pyunits.convert(p_sat_dT, to_units=dp_units)
