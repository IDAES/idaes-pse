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

The Properties of Gases & Liquids, 4th Edition
Reid, Prausnitz and Polling, 1987, McGraw-Hill

All parameter indicies and units based on conventions used by the source
"""
from pyomo.environ import exp, log, Var, units as pyunits


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class cp_mol_ig_comp():
    def build_parameters(cobj):
        cobj.cp_mol_ig_comp_coeff_A = Var(
                initialize=cobj.config.parameter_data["cp_mol_ig_comp_coeff"]["A"],
                doc="Parameter A for ideal gas molar heat capacity",
                units=pyunits.J/pyunits.mol/pyunits.K)
        cobj.cp_mol_ig_comp_coeff_B = Var(
                initialize=cobj.config.parameter_data["cp_mol_ig_comp_coeff"]["B"],
                doc="Parameter B for ideal gas molar heat capacity",
                units=pyunits.J/pyunits.mol/pyunits.K**2)
        cobj.cp_mol_ig_comp_coeff_C = Var(
                initialize=cobj.config.parameter_data["cp_mol_ig_comp_coeff"]["C"],
                doc="Parameter C for ideal gas molar heat capacity",
                units=pyunits.J/pyunits.mol/pyunits.K**3)
        cobj.cp_mol_ig_comp_coeff_D = Var(
                initialize=cobj.config.parameter_data["cp_mol_ig_comp_coeff"]["D"],
                doc="Parameter D for ideal gas molar heat capacity",
                units=pyunits.J/pyunits.mol/pyunits.K**4)

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
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
            cp_mol_ig_comp.build_parameters(cobj)

        base_units = cobj.parent_block().get_metadata().default_units
        h_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        cobj.enth_mol_form_vap_comp_ref = Var(
                initialize=cobj.config.parameter_data[
                    "enth_mol_form_vap_comp_ref"],
                doc="Vapor phase molar heat of formation @ Tref",
                units=h_units)

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
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
            cp_mol_ig_comp.build_parameters(cobj)

        base_units = cobj.parent_block().get_metadata().default_units
        s_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1 *
                   base_units["temperature"]**-1)

        cobj.entr_mol_form_vap_comp_ref = Var(
                initialize=cobj.config.parameter_data[
                    "entr_mol_form_vap_comp_ref"],
                doc="Vapor phase molar entropy of formation @ Tref",
                units=s_units)

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
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff_A = Var(
                initialize=cobj.config.parameter_data[
                    "pressure_sat_comp_coeff"]["A"],
                doc="Coefficient A for calculating Psat",
                units=None)
        cobj.pressure_sat_comp_coeff_B = Var(
                initialize=cobj.config.parameter_data[
                    "pressure_sat_comp_coeff"]["B"],
                doc="Coefficient B for calculating Psat",
                units=None)
        cobj.pressure_sat_comp_coeff_C = Var(
                initialize=cobj.config.parameter_data[
                    "pressure_sat_comp_coeff"]["C"],
                doc="Coefficient C for calculating Psat",
                units=None)
        cobj.pressure_sat_comp_coeff_D = Var(
                initialize=cobj.config.parameter_data[
                    "pressure_sat_comp_coeff"]["D"],
                doc="Coefficient D for calculating Psat",
                units=None)

    def return_expression(b, cobj, T, dT=False):
        if dT:
            return pressure_sat_comp.dT_expression(b, cobj, T)

        x = 1 - T/cobj.temperature_crit

        return (exp((1-x)**-1 * (cobj.pressure_sat_comp_coeff_A*x +
                                 cobj.pressure_sat_comp_coeff_B*x**1.5 +
                                 cobj.pressure_sat_comp_coeff_C*x**3 +
                                 cobj.pressure_sat_comp_coeff_D*x**6)) *
                cobj.pressure_crit)

    def dT_expression(b, cobj, T):
        x = 1 - T/cobj.temperature_crit

        return (-pressure_sat_comp.return_expression(b, cobj, T) *
                ((cobj.pressure_sat_comp_coeff_A +
                  1.5*cobj.pressure_sat_comp_coeff_B*x**0.5 +
                  3*cobj.pressure_sat_comp_coeff_C*x**2 +
                  6*cobj.pressure_sat_comp_coeff_D*x**5)/T +
                 (cobj.temperature_crit/T**2) *
                 (cobj.pressure_sat_comp_coeff_A*x +
                  cobj.pressure_sat_comp_coeff_B*x**1.5 +
                  cobj.pressure_sat_comp_coeff_C*x**3 +
                  cobj.pressure_sat_comp_coeff_D*x**6)))
