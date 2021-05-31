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
Custom Methods for calculating pure component properties of MEA, H2O:
"""
from pyomo.environ import log, exp, Var, units as pyunits

from idaes.core.util.misc import set_param_from_config 

# Specific Heat capacity Correlation for MEA (J/mol.K)
# Data Source [6]
# Eqn 6, Table 1
class cp_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        cobj.cp_mol_liq_comp_coeff_a0 = Var(
            doc="Parameter a0 for liquid phase heat capacity",
            units=pyunits.J*pyunits.mol**-1*pyunits.K**-1)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="a0")

        cobj.cp_mol_liq_comp_coeff_a1 = Var(
            doc="Parameter a1 for liquid phase heat capacity",
            units=pyunits.J*pyunits.mol**-1*pyunits.K**-2)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="a1")


    @staticmethod
    def return_expression(b, cobj, T):
        T = pyunits.convert(T, to_units=pyunits.K)
        cp = cobj.cp_mol_liq_comp_coeff_a0 + (cobj.cp_mol_liq_comp_coeff_a1*T)
        units = b.params.get_metadata().derived_units
        return pyunits.convert(cp, units["heat_capacity_mole"])

# Specific Enthalpy Calculation for MEA (J/mol)
# (based on Specific Heat capacity Correlation)
class enth_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff_a0"):
            cp_mol_liq_comp.build_parameters(cobj)

        units = cobj.parent_block().get_metadata().derived_units

        cobj.enth_mol_form_liq_comp_ref = Var(
                doc="Liquid phase molar heat of formation @ Tref",
                units=units["energy_mole"])
        set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

    @staticmethod
    def return_expression(b, cobj, T):
        T = pyunits.convert(T, to_units=pyunits.K)
        Tref = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)
        specific_enthalpy = (cobj.cp_mol_liq_comp_coeff_a0*(T-Tref)) + \
            (cobj.cp_mol_liq_comp_coeff_a1*(T-Tref)**2/2) + \
                cobj.enth_mol_form_liq_comp_ref
        return specific_enthalpy

# Molar Volume Correlation for MEA, H2O (cm3/mol)
# Data Source [5] 
# MEA - Eqn 16, Table 3
# H2O - Eqn 25
class vol_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        cobj.mol_vol_liq_comp_coeff_1 = Var(
                doc="Parameter 1 for liquid phase molar volume",
                units=pyunits.g/pyunits.cm**3/pyunits.K**2)
        set_param_from_config(cobj, param="mol_vol_liq_comp_coeff", index="1")

        cobj.mol_vol_liq_comp_coeff_2 = Var(
                doc="Parameter 2 for liquid phase molar volume",
                units=pyunits.g/pyunits.cm**3/pyunits.K)
        set_param_from_config(cobj, param="mol_vol_liq_comp_coeff", index="2")

        cobj.mol_vol_liq_comp_coeff_3 = Var(
                doc="Parameter 3 for liquid phase molar volume",
                units=pyunits.g/pyunits.cm**3)
        set_param_from_config(cobj, param="mol_vol_liq_comp_coeff", index="3")
        

    @staticmethod
    def return_expression(b, cobj, T):
        T = pyunits.convert(T, to_units=pyunits.K)        
        mw = pyunits.convert(cobj.mw, to_units=pyunits.g/pyunits.mol)

        vol = mw/((cobj.mol_vol_liq_comp_coeff_1*T*T)+\
                  (cobj.mol_vol_liq_comp_coeff_2*T)+\
                      cobj.mol_vol_liq_comp_coeff_3)

        return vol
    
# Molar Density Correlation for MEA, H2O (based on molar volume)
class dens_mol_liq_comp():
    
    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "mol_vol_liq_comp_coeff_1"):
            vol_mol_liq_comp.build_parameters(cobj)
    
    @staticmethod
    def return_expression(b, cobj, T):
        vol = vol_mol_liq_comp.return_expression(b, cobj, T)
        units = b.params.get_metadata().derived_units
        return pyunits.convert(1/vol, units["density_mole"])
    
# Extended Antoines Correlation for H2O (optional: for better accuracy)
# Data source [4]
# Eqn 7-4, Table 7.4-2
class pressure_sat_comp_option_1():

    @staticmethod
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff_C1 = Var(
                doc="Extended Antoine coefficient C1 for calculating Psat",
                units=None)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C1")

        cobj.pressure_sat_comp_coeff_C2 = Var(
                doc="Extended Antoine coefficient C2 for calculating Psat",
                units=pyunits.K)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C2")

        cobj.pressure_sat_comp_coeff_C3 = Var(
                doc="Extended Antoine coefficient C3 for calculating Psat",
                units=pyunits.K)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C3")
        cobj.pressure_sat_comp_coeff_C4 = Var(
                doc="Extended Antoine coefficient C4 for calculating Psat",
                units=1/(pyunits.K))
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C4")
        cobj.pressure_sat_comp_coeff_C5 = Var(
                doc="Extended Antoine coefficient C5 for calculating Psat",
                units=None)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C5")
        cobj.pressure_sat_comp_coeff_C6 = Var(
                doc="Extended Antoine coefficient C6 for calculating Psat",
                units=1/(pyunits.K)**2)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C6")
        cobj.pressure_sat_comp_coeff_C7 = Var(
                doc="Extended Antoine coefficient C7 for calculating Psat",
                units=None)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C7")

    @staticmethod
    def return_expression(b, cobj, T, dT=False):
        T = pyunits.convert(T, to_units=pyunits.K)
        if dT:
            return pressure_sat_comp.dT_expression(b, cobj, T)
        psat = (exp(cobj.pressure_sat_comp_coeff_C1 + \
                    (cobj.pressure_sat_comp_coeff_C2/(T + \
                    cobj.pressure_sat_comp_coeff_C3)) \
                   + (cobj.pressure_sat_comp_coeff_C4*T) + \
                     (cobj.pressure_sat_comp_coeff_C5*log(T*1/pyunits.K)) + \
                     (cobj.pressure_sat_comp_coeff_C6*(T)**\
                      cobj.pressure_sat_comp_coeff_C7)))*1000*pyunits.Pa

        units = b.params.get_metadata().derived_units
        return pyunits.convert(psat, to_units=pyunits.Pa)
    
    @staticmethod
    def dT_expression(b, cobj, T):            
        p_sat_dT = pressure_sat_comp.return_expression(b, cobj, T)*\
            ((-cobj.pressure_sat_comp_coeff_C2/\
              (T+cobj.pressure_sat_comp_coeff_C3)**2) +\
             cobj.pressure_sat_comp_coeff_C4 - \
                 (cobj.pressure_sat_comp_coeff_C5/T) +\
                 (cobj.pressure_sat_comp_coeff_C6*\
                  cobj.pressure_sat_comp_coeff_C7*\
                  T**(cobj.pressure_sat_comp_coeff_C7-1)))
        units = b.params.get_metadata().derived_units
        dp_units = units["pressure"]/units["temperature"]
        return pyunits.convert(p_sat_dT, to_units=dp_units)